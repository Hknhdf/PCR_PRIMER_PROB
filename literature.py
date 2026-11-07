"""
literature.py
Literatur / yayinlanmis assay karsilastirma ve assay kutuphanesi.

Amac: kullanicinin yapistirdigi yayin primer/prob setini referans genom uzerinde
konumlandirmak, uygulamanin aktif assay'i ile örtüsme/benzerlik/uyum skorunu
hesaplamak. Kutuphane baslangicta bos gelir; yanlis/dogrulanmamis sekans
uydurmamak icin ornek sekanslar gomulmez.
"""

from __future__ import annotations

import io
import json
import re
from datetime import datetime

import pandas as pd
import requests
from Bio import Align

from pcr_engine import calculate_gc_content, calculate_tm, revcomp
from specificity import in_silico_pcr


def global_identity(a: str, b: str):
    """Global hizalama ile yuzde identite."""
    a = str(a).upper()
    b = str(b).upper()
    if not a or not b:
        return 0.0
    if a == b:
        return 100.0
    if len(a) == len(b):
        matches = sum(1 for x, y in zip(a, b) if x == y)
        return round(100.0 * matches / len(a), 1)
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    try:
        aln = aligner.align(a, b)[0]
        ref_aln, sub_aln = str(aln[0]), str(aln[1])
    except Exception:
        return 0.0
    matches = compared = 0
    for x, y in zip(ref_aln, sub_aln):
        if x == "-" or y == "-":
            continue
        compared += 1
        if x == y:
            matches += 1
    return round(100.0 * matches / compared, 1) if compared else 0.0


def search_published_papers(query: str, limit: int = 25, min_citations: int = 0,
                            year_from: int | None = None, year_to: int | None = None):
    """Yayınlanmış makaleleri arar ve atıf sayısına göre sıralar.

    Önemli: bu fonksiyon yalnızca bibliyografik metadata döndürür. Makalelerdeki
    primer/prob dizileri ayrı bir veritabanından uydurulamaz; diziler kullanıcı
    tarafından makale tablosundan kopyalanmalıdır. Atıf sayısı Semantic Scholar
    üzerinden, yedek olarak Crossref üzerinden alınır.
    """
    query = str(query or "").strip()
    if not query:
        return {"ok": False, "papers": [], "error": "Arama terimi boş.", "source": ""}
    limit = max(1, min(int(limit or 25), 100))
    headers = {"User-Agent": "PCR-Primer-Probe/1.0 (in-silico assay literatur arama)"}

    # 1) Semantic Scholar: citationCount + DOI/PMID + açık erişim bağlantısı.
    try:
        params = {
            "query": query,
            "limit": limit,
            "fields": "title,authors,year,venue,externalIds,url,citationCount,openAccessPdf,abstract",
        }
        resp = requests.get(
            "https://api.semanticscholar.org/graph/v1/paper/search",
            params=params, headers=headers, timeout=30,
        )
        resp.raise_for_status()
        payload = resp.json()
        papers = []
        for item in payload.get("data", []) or []:
            cites = int(item.get("citationCount") or 0)
            if cites < int(min_citations or 0):
                continue
            try:
                year = int(item.get("year")) if item.get("year") else None
            except Exception:
                year = None
            if year_from and isinstance(year, int) and year < int(year_from):
                continue
            if year_to and isinstance(year, int) and year > int(year_to):
                continue
            external = item.get("externalIds") or {}
            authors = ", ".join(
                str(a.get("name", "")).strip()
                for a in (item.get("authors") or [])
                if a.get("name")
            )[:240]
            pdf = item.get("openAccessPdf") or {}
            papers.append({
                "title": str(item.get("title") or "").strip(),
                "authors": authors,
                "year": year,
                "venue": str(item.get("venue") or "").strip(),
                "doi": str(external.get("DOI") or "").strip(),
                "pmid": str(external.get("PubMed") or "").strip(),
                "url": str(item.get("url") or "").strip(),
                "pdf_url": str(pdf.get("url") or "").strip(),
                "citation_count": cites,
                "source": "Semantic Scholar",
                "abstract": str(item.get("abstract") or "").strip(),
            })
        papers.sort(key=lambda p: (-int(p.get("citation_count") or 0), -(p.get("year") or 0)))
        if papers:
            return {"ok": True, "papers": papers[:limit], "source": "Semantic Scholar", "error": ""}
        return {"ok": False, "papers": [], "error": "Sonuç bulunamadı.", "source": "Semantic Scholar"}
    except Exception as semantic_exc:
        semantic_error = str(semantic_exc)

    # 2) Crossref yedeği: is-referenced-by-count atıf sayımı verir.
    try:
        params = {
            "query.bibliographic": query,
            "rows": max(limit * 2, limit),
            "select": "title,author,container-title,published,DOI,URL,is-referenced-by-count,abstract,publisher,link,resource",
        }
        resp = requests.get(
            "https://api.crossref.org/works",
            params=params, headers=headers, timeout=30,
        )
        resp.raise_for_status()
        payload = resp.json()
        papers = []
        for item in (payload.get("message") or {}).get("items", []) or []:
            cites = int(item.get("is-referenced-by-count") or 0)
            if cites < int(min_citations or 0):
                continue
            published = item.get("published") or item.get("issued") or {}
            year = None
            try:
                year = int((published.get("date-parts") or [[None]])[0][0])
            except Exception:
                year = None
            if year_from and isinstance(year, int) and year < int(year_from):
                continue
            if year_to and isinstance(year, int) and year > int(year_to):
                continue
            title = " ".join(str(x) for x in (item.get("title") or []))
            authors = ", ".join(
                " ".join(str(a.get(k, "")).strip() for k in ("given", "family") if a.get(k))
                for a in (item.get("author") or [])
            ).strip()[:240]
            pdf_url = ""
            for link in item.get("link") or []:
                try:
                    url = str((link or {}).get("URL") or "").strip()
                    if url and (url.lower().endswith(".pdf") or "/pdf" in url.lower()):
                        pdf_url = url
                        break
                except Exception:
                    continue
            papers.append({
                "title": title.strip(),
                "authors": authors,
                "year": year,
                "venue": " ".join(str(x) for x in (item.get("container-title") or [])),
                "doi": str(item.get("DOI") or "").strip(),
                "pmid": "",
                "url": str(item.get("URL") or "").strip(),
                "pdf_url": pdf_url,
                "citation_count": cites,
                "source": "Crossref",
                "abstract": str(item.get("abstract") or "").strip(),
            })
        papers = [p for p in papers if p.get("title")]
        papers.sort(key=lambda p: (-int(p.get("citation_count") or 0), -(p.get("year") or 0)))
        if papers:
            return {"ok": True, "papers": papers[:limit], "source": "Crossref", "error": ""}
        return {
            "ok": False,
            "papers": [],
            "error": f"Semantic Scholar: {semantic_error}; Crossref: sonuç yok.",
            "source": "Crossref",
        }
    except Exception as crossref_exc:
        return {
            "ok": False,
            "papers": [],
            "error": f"Semantic Scholar: {semantic_error}; Crossref: {crossref_exc}",
            "source": "",
        }


def _extract_sequence_candidates_from_text(text: str):
    """Tam metinden primer/prob adayı dizileri bulur (pozisyon bağlamı ile)."""
    if not text:
        return []
    # PDF satır kırılmaları/boşluklar oligo'yu bölebilir; önce metni normalize et.
    compact = re.sub(r"([ACGTacgtryskmbdhvn])\s+([ACGTacgtryskmbdhvn])",
                     lambda m: m.group(1) + m.group(2), text)
    pattern = re.compile(r"\b([ACGTRYSWKMBDHVN]{15,35})\b", re.IGNORECASE)
    seen = []
    out = []
    for m in pattern.finditer(compact):
        seq = m.group(1).upper()
        # İngilizce metinde ACGT-benzeri sahte pozitifleri ele.
        if set(seq) <= set("ACGT") and len(seq) < 16:
            continue
        if len(set(seq)) < 3:
            continue
        ctx0 = max(0, m.start() - 240)
        ctx1 = min(len(compact), m.end() + 240)
        context = re.sub(r"\s+", " ", compact[ctx0:ctx1]).strip()
        low = context.lower()
        score = 0
        if "forward" in low or "forward primer" in low or "f primer" in low or "5'" in low or "5’" in low:
            score += 3
        if "reverse" in low or "reverse primer" in low or "r primer" in low or "3'" in low or "3’" in low:
            score += 3
        if "probe" in low or "taqman" in low or "hybridization" in low:
            score += 4
        if "primer" in low or "oligonucleotide" in low:
            score += 2
        if "table" in low or "sequence" in low:
            score += 1
        if score <= 0:
            continue
        role = "aday"
        if "probe" in low or "taqman" in low:
            role = "prob"
        elif "reverse" in low:
            role = "reverse"
        elif "forward" in low:
            role = "forward"
        if (seq, context) not in seen:
            seen.append((seq, context))
            out.append({
                "sequence": seq,
                "role": role,
                "score": score,
                "context": context[:420],
            })
    return sorted(out, key=lambda x: (-x["score"], len(x["sequence"])))


def extract_published_assay_candidates(pdf_url: str = "", doi: str = "", pmid: str = "",
                                       max_candidates: int = 30):
    """Açık erişim tam metinden primer/prob adayı dizileri çıkarır.

    Yalnızca erişilebilir PDF veya Europe PMC tam metni kullanılmalıdır. Dizi uydurulmaz;
    benzergah/rol sınıflandırması zayıfsa 'aday' olarak döner ve kullanıcı onaylar.
    """
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)",
        "Accept": "application/pdf,*/*;q=0.8",
    }
    text = ""
    source = ""
    notes = []

    if pdf_url:
        try:
            resp = requests.get(str(pdf_url), headers=headers, timeout=45)
            resp.raise_for_status()
            ctype = (resp.headers.get("content-type") or "").lower()
            if "pdf" in ctype or resp.content[:4] == b"%PDF":
                from pypdf import PdfReader

                reader = PdfReader(io.BytesIO(resp.content))
                parts = []
                for page in reader.pages:
                    try:
                        parts.append(page.extract_text() or "")
                    except Exception:
                        continue
                text = "\n".join(parts)
                source = str(pdf_url)
            else:
                notes.append("PDF bağlantısı HTML döndürdü.")
        except Exception as exc:
            notes.append(f"PDF okunamadı: {exc}")

    if not text and (doi or pmid):
        try:
            q = f'DOI:"{doi}"' if doi else f'EXT_ID:{pmid}'
            resp = requests.get(
                "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                params={"query": q, "format": "json", "resultType": "core", "pageSize": 1},
                headers={"User-Agent": "PCR-Primer-Probe/1.0"},
                timeout=30,
            )
            resp.raise_for_status()
            results = (resp.json().get("resultList") or {}).get("result") or []
            pmcid = results[0].get("pmcid") if results else ""
            if pmcid:
                full = requests.get(
                    f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML",
                    headers={"User-Agent": "PCR-Primer-Probe/1.0"},
                    timeout=45,
                )
                full.raise_for_status()
                text = full.text
                source = f"Europe PMC {pmcid}"
            else:
                notes.append("Europe PMC'de tam metin (PMCID) bulunamadı.")
        except Exception as exc:
            notes.append(f"Europe PMC tam metni alınamadı: {exc}")

    candidates = _extract_sequence_candidates_from_text(text) if text else []
    return {
        "ok": bool(candidates),
        "candidates": candidates[:max_candidates],
        "source": source,
        "n_chars": len(text),
        "notes": notes,
    }


def locate_assay(ref_seq: str, forward: str, reverse: str, probe: str | None = None,
                 max_mismatches: int = 4, min_amp: int = 30, max_amp: int = 2000):
    """Yayinlanmis assay'i referans genom uzerinde konumlandir."""
    forward = str(forward).upper().strip()
    reverse = str(reverse).upper().strip()
    probe = str(probe).upper().strip() if probe else None
    hits = in_silico_pcr(
        forward, reverse, ref_seq, probe,
        max_mismatches=max_mismatches, min_amp=min_amp, max_amp=max_amp,
    )
    if not hits:
        return {"found": False, "message": "Bu primer çifti referansta amplikon oluşturmadı.",
                "hits": [], "best": None}
    best = hits[0]
    best_info = best.as_dict()
    best_info.update({
        "forward_seq": forward,
        "reverse_seq": reverse,
        "probe_seq": probe,
        "forward_tm": calculate_tm(forward),
        "reverse_tm": calculate_tm(reverse),
        "probe_tm": calculate_tm(probe) if probe else None,
        "forward_gc": calculate_gc_content(forward),
        "reverse_gc": calculate_gc_content(reverse),
        "probe_gc": calculate_gc_content(probe) if probe else None,
        "amplicon": best.amplicon,
    })
    return {"found": True, "message": "", "best": best_info, "hits": [h.as_dict() for h in hits]}


def _overlap(a0, a1, b0, b1):
    lo = max(int(a0), int(b0))
    hi = min(int(a1), int(b1))
    return max(0, hi - lo)


def compare_assays(active: dict, published: dict):
    """
    active: {'forward','reverse','probe','amplicon_start','amplicon_end','amplicon'}
    published: locate_assay(...)['best']
    """
    if not active or not published:
        return {"ok": False, "message": "Karşılaştırma için iki assay gerekir."}
    a_f = str(active.get("forward", "")).upper()
    a_r = str(active.get("reverse", "")).upper()
    a_p = str(active.get("probe", "") or "").upper()
    p_f = str(published.get("forward_seq", "")).upper()
    p_r = str(published.get("reverse_seq", "")).upper()
    p_p = str(published.get("probe_seq", "") or "").upper()

    a0, a1 = int(active.get("amplicon_start", 0)), int(active.get("amplicon_end", 0))
    p0 = int(published.get("amplicon_start", 0))
    p1 = int(published.get("amplicon_end", 0))
    ov = _overlap(a0, a1, p0, p1)
    union = max(a1, p1) - min(a0, p0) if max(a1, p1) > min(a0, p0) else 0
    jaccard = ov / union if union else 0.0

    a_amp = str(active.get("amplicon", "")).upper()
    p_amp = str(published.get("amplicon", "")).upper()
    amp_identity = global_identity(a_amp, p_amp) if a_amp and p_amp else 0.0

    f_ident = global_identity(a_f, p_f) if a_f and p_f else None
    r_ident = global_identity(a_r, p_r) if a_r and p_r else None
    p_ident = global_identity(a_p, p_p) if a_p and p_p else None
    oligo_idents = [x for x in (f_ident, r_ident, p_ident) if x is not None]
    oligo_mean = sum(oligo_idents) / len(oligo_idents) if oligo_idents else 0.0

    f_tm_diff = abs((calculate_tm(a_f) or 0) - (published.get("forward_tm") or 0)) if a_f else None
    r_tm_diff = abs((calculate_tm(a_r) or 0) - (published.get("reverse_tm") or 0)) if a_r else None
    p_tm_diff = abs((calculate_tm(a_p) or 0) - (published.get("probe_tm") or 0)) if a_p and published.get("probe_tm") else None

    score = 0.30 * (100 * jaccard) + 0.30 * amp_identity + 0.40 * oligo_mean
    if jaccard < 0.1:
        score *= 0.5
    return {
        "ok": True,
        "amplicon_overlap_bp": ov,
        "amplicon_jaccard": round(jaccard, 3),
        "amplicon_identity": amp_identity,
        "forward_identity": f_ident,
        "reverse_identity": r_ident,
        "probe_identity": p_ident,
        "oligo_mean_identity": round(oligo_mean, 1),
        "forward_tm_diff": round(f_tm_diff, 1) if f_tm_diff is not None else None,
        "reverse_tm_diff": round(r_tm_diff, 1) if r_tm_diff is not None else None,
        "probe_tm_diff": round(p_tm_diff, 1) if p_tm_diff is not None else None,
        "concordance_score": round(score, 1),
    }


def comparison_table(comparison: dict):
    if not comparison.get("ok"):
        return pd.DataFrame()
    rows = [
        ("Amplikon örtüşmesi (bp)", comparison["amplicon_overlap_bp"]),
        ("Amplikon Jaccard", comparison["amplicon_jaccard"]),
        ("Amplikon identitesi (%)", comparison["amplicon_identity"]),
        ("Forward identitesi (%)", comparison["forward_identity"]),
        ("Reverse identitesi (%)", comparison["reverse_identity"]),
        ("Prob identitesi (%)", comparison["probe_identity"]),
        ("Forward Tm farkı (°C)", comparison["forward_tm_diff"]),
        ("Reverse Tm farkı (°C)", comparison["reverse_tm_diff"]),
        ("Prob Tm farkı (°C)", comparison["probe_tm_diff"]),
        ("Uyum skoru (0-100)", comparison["concordance_score"]),
    ]
    return pd.DataFrame(rows, columns=["Metrik", "Değer"])


def library_entry(name: str, reference: str, forward: str, reverse: str, probe: str = "",
                  notes: str = "", source: str = "", tags=()):
    return {
        "name": str(name).strip(),
        "reference": str(reference).strip(),
        "source": str(source).strip(),
        "forward": str(forward).upper().strip(),
        "reverse": str(reverse).upper().strip(),
        "probe": str(probe).upper().strip(),
        "notes": str(notes).strip(),
        "tags": list(tags),
        "added": datetime.now().isoformat(timespec="seconds"),
    }


def validate_entry(entry: dict):
    missing = [k for k in ("name", "forward", "reverse") if not str(entry.get(k, "")).strip()]
    if missing:
        return False, "Eksik alan: " + ", ".join(missing)
    for key in ("forward", "reverse", "probe"):
        seq = str(entry.get(key, "") or "").upper()
        if seq and (len(seq) < 15 or len(seq) > 40):
            return False, f"{key} uzunluğu 15-40 nt olmalı."
        if seq and set(seq) - set("ACGTN"):
            return False, f"{key} içinde ACGTN dışı karakter var."
    return True, ""


def library_to_json(entries):
    return json.dumps(list(entries), ensure_ascii=False, indent=2)


def library_from_json(text: str):
    try:
        data = json.loads(text)
    except Exception as exc:
        return [], f"JSON okunamadı: {exc}"
    if not isinstance(data, list):
        return [], "JSON bir liste olmalı."
    out = []
    for item in data:
        if not isinstance(item, dict):
            continue
        ok, msg = validate_entry(item)
        if ok:
            out.append(item)
    return out, ""


def library_dataframe(entries):
    if not entries:
        return pd.DataFrame(columns=["Ad", "Kaynak", "Forward", "Reverse", "Prob", "Not"])
    return pd.DataFrame([
        {
            "Ad": e.get("name", ""),
            "Kaynak": e.get("source", e.get("reference", "")),
            "Forward": e.get("forward", ""),
            "Reverse": e.get("reverse", ""),
            "Prob": e.get("probe", ""),
            "Not": e.get("notes", ""),
        }
        for e in entries
    ])
