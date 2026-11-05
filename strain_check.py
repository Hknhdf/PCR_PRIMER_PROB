"""
strain_check.py
Türkiye suşu uyum kontrolü ve sipariş öncesi QC.

Amaç:
  - NCBI'da ülke filtresiyle (varsayılan Turkey) bulunan dizileri veya kullanıcının
    yüklediği FASTA dizilerini aktif primer/prob setiyle karşılaştırmak.
  - Primer/prob bağlanma bölgelerindeki mismatch, 3' uç mismatch'i ve amplikon
    varlığını raporlamak.
  - Sipariş öncesi termodinamik + özgüllük + Türkiye kapsaması + manuel teyit
    kontrollerini tek bir karar tablosunda toplamak.

Bu modül in-silico ön kontroldür; ıslak laboratuvar validasyonunun yerine geçmez.
"""

from __future__ import annotations

import copy
import io
import re
from datetime import datetime

from Bio import SeqIO

from genotype import analyze_genotypes, fetch_subjects
from thermo import analyze_oligo, dimer_matrix

PASS = "PASS"
WARN = "WARN"
FAIL = "FAIL"
BLOCKED = "BLOCKED"


def _clean_sequence(text: str) -> str:
    seq = str(text or "").upper().replace("U", "T")
    return re.sub(r"[^ACGTN]", "N", seq)


def parse_fasta_text(text: str) -> list[dict]:
    """FASTA veya başlıksız ham diziyi [{'id','description','seq'}] biçimine çevirir."""
    text = str(text or "").strip()
    if not text:
        return []
    records: list[dict] = []
    try:
        for rec in SeqIO.parse(io.StringIO(text), "fasta"):
            seq = _clean_sequence(str(rec.seq))
            if len(seq) >= 20:
                records.append({
                    "id": str(rec.id) or f"MANUEL-{len(records) + 1}",
                    "description": str(rec.description or rec.id),
                    "seq": seq,
                })
    except Exception:
        records = []
    if records:
        return records

    # FASTA başlığı yoksa: boş satırlarla ayrılmış ham diziler
    for i, block in enumerate(re.split(r"\n\s*\n", text), start=1):
        seq = _clean_sequence(block)
        if len(seq) >= 20:
            records.append({
                "id": f"MANUEL-{i}",
                "description": f"Manuel dizi {i}",
                "seq": seq,
            })
    return records


def fetch_turkey_subjects(organism: str, reference_amp: str, n: int = 20,
                          keywords=(), accession: str = "", country: str = "Turkey"):
    """NCBI'da organizma + ülke filtresiyle dizileri çeker.

    Ülke metadata'sı boş olan kayıtlar bu yolla gelmez; kullanıcı FASTA yüklemelidir.
    """
    if not organism:
        return []
    return fetch_subjects(
        organism,
        exclude_accessions=[accession] if accession else [],
        reference_seq=reference_amp,
        n=int(n),
        keywords=tuple(keywords or ()),
        min_subject_len=60,
        country=country or None,
    )


def classify_record(rec: dict, expected_size: int, has_probe: bool = True) -> tuple[str, str]:
    """Tek bir dizinin primer/prob uyumunu PASS/WARN/FAIL olarak sınıflandırır."""
    if rec.get("error") and not rec.get("found") and not rec.get("partial"):
        return FAIL, str(rec["error"])[:160]
    if not rec.get("found") and not rec.get("partial"):
        return FAIL, "Amplikon bulunamadı."

    f = rec.get("F") or {}
    r = rec.get("R") or {}
    p = rec.get("P") or {}
    f_mm = f.get("n_mismatches")
    r_mm = r.get("n_mismatches")
    p_mm = p.get("n_mismatches")
    three = int(rec.get("three_prime_mm") or 0)
    total = rec.get("total_mm")
    identity = rec.get("identity")
    notes: list[str] = []

    if rec.get("partial"):
        notes.append("Tam amplikon yok; yalnızca bağımsız bağlanma bölgeleri bulundu.")
    if three > 0:
        notes.append("Primer 3' ucunda mismatch var (uzama riski).")
    if any(v is not None and int(v) >= 3 for v in (f_mm, r_mm)):
        notes.append("Primer bağlanma bölgesinde ≥3 mismatch.")
    if has_probe and p_mm is not None and int(p_mm) >= 2:
        notes.append("Prob bölgesinde ≥2 mismatch (sinyal kaybı riski).")
    if rec.get("amplicon_size") and expected_size:
        diff = abs(int(rec["amplicon_size"]) - int(expected_size)) / max(1, int(expected_size))
        if diff > 0.25:
            notes.append("Amplikon boyu beklenenden >%25 farklı.")
    if identity is not None and identity < 0.90:
        notes.append("Amplikon identitesi %90'ın altında.")
    elif identity is not None and identity < 0.97:
        notes.append("Amplikon identitesi %97'nin altında.")

    critical = (
        three > 0
        or any(v is not None and int(v) >= 3 for v in (f_mm, r_mm))
        or (has_probe and p_mm is not None and int(p_mm) >= 2)
        or (identity is not None and identity < 0.90)
    )
    if critical:
        return FAIL, " ".join(notes) or "Kritik mismatch."

    if rec.get("partial"):
        return WARN, " ".join(notes) or "Kısmi eşleşme."

    if total == 0 and (identity is None or identity >= 0.99) and not notes:
        return PASS, "Tam eşleşme."
    if total is not None and int(total) <= 2:
        return WARN, " ".join(notes) or "1–2 mismatch; laboratuvarda doğrulanmalı."
    return FAIL, " ".join(notes) or "Toplam mismatch çok yüksek."


def evaluate_turkey_strains(active: dict, subjects: list[dict], max_mismatches: int = 4) -> dict:
    """Aktif assay'i verilen Türkiye/yerel dizilerle karşılaştırır."""
    if not active:
        return {"error": "Aktif assay yok.", "rows": [], "n_total": 0}
    if not subjects:
        return {"error": "Değerlendirilecek dizi yok.", "rows": [], "n_total": 0}

    result = analyze_genotypes(
        active["ref_seq"],
        active["amplicon_start"],
        active["amplicon_end"],
        active["forward"],
        active["reverse"],
        active.get("probe") or None,
        subjects,
        probe_offset=active.get("probe_offset"),
        max_mismatches=int(max_mismatches),
        include_reference=False,
    )
    if result.get("error"):
        return {**result, "rows": [], "n_total": 0}

    has_probe = bool(active.get("probe"))
    rows = []
    for rec in result.get("records", []):
        status, note = classify_record(
            rec, expected_size=int(active.get("amplicon_size") or 0), has_probe=has_probe
        )
        rows.append({
            "accession": rec.get("id", "?"),
            "description": rec.get("description", ""),
            "status": status,
            "note": note,
            "amplicon": "var" if rec.get("found") else ("kısmi" if rec.get("partial") else "yok"),
            "size": rec.get("amplicon_size"),
            "identity": rec.get("identity"),
            "f_mm": (rec.get("F") or {}).get("n_mismatches"),
            "r_mm": (rec.get("R") or {}).get("n_mismatches"),
            "p_mm": (rec.get("P") or {}).get("n_mismatches"),
            "three_mm": rec.get("three_prime_mm"),
        })

    result["rows"] = rows
    result["n_total"] = len(rows)
    result["n_pass"] = sum(1 for r in rows if r["status"] == PASS)
    result["n_warn"] = sum(1 for r in rows if r["status"] == WARN)
    result["n_fail"] = sum(1 for r in rows if r["status"] == FAIL)
    result["n_three_prime_risk"] = sum(1 for r in rows if (r.get("three_mm") or 0) > 0)
    if result["n_fail"]:
        result["overall"] = FAIL
    elif result["n_warn"]:
        result["overall"] = WARN
    else:
        result["overall"] = PASS
    return result


def turkey_rows(result: dict) -> list[dict]:
    """UI/CSV için sade Türkiye sonuç tablosu."""
    return [{
        "Accession": r.get("accession"),
        "Açıklama": r.get("description"),
        "Durum": r.get("status"),
        "Amplikon": r.get("amplicon"),
        "Boyut (bp)": r.get("size"),
        "Identite": r.get("identity"),
        "F mm": r.get("f_mm"),
        "R mm": r.get("r_mm"),
        "P mm": r.get("p_mm"),
        "3' mm": r.get("three_mm"),
        "Not": r.get("note"),
    } for r in result.get("rows", [])]


def consensus_to_active(active: dict, candidate: dict) -> dict:
    """Consensus adayını, diğer sekmelerin kullanacağı aktif assay biçimine çevirir.

    candidate koordinatları, aktif assay'ın referans amplikonu içinde 1-based'tir.
    """
    ref_seq = str(active["ref_seq"]).upper()
    amp_start = int(active["amplicon_start"])
    f_local0 = amp_start + int(candidate["forward_start"]) - 1
    r_local1 = amp_start + int(candidate["reverse_end"])
    p_local0 = amp_start + int(candidate["probe_start"]) - 1
    amp_seq = ref_seq[f_local0:r_local1]
    if not amp_seq:
        raise ValueError("Consensus adayı için amplikon oluşturulamadı.")

    forward = str(candidate["forward"]).upper()
    reverse = str(candidate["reverse"]).upper()
    probe = str(candidate.get("probe") or "").upper()
    ref_offset = int(active.get("ref_offset", 0))
    new = copy.deepcopy(active)
    new.update({
        "ref_seq": ref_seq,
        "forward": forward,
        "reverse": reverse,
        "probe": probe,
        "forward_tm": candidate.get("forward_tm"),
        "reverse_tm": candidate.get("reverse_tm"),
        "probe_tm": candidate.get("probe_tm"),
        "forward_gc": candidate.get("forward_gc"),
        "reverse_gc": candidate.get("reverse_gc"),
        "probe_gc": candidate.get("probe_gc"),
        "amplicon_start": f_local0,
        "amplicon_end": r_local1,
        "amplicon": amp_seq,
        "amplicon_size": len(amp_seq),
        "start": f_local0,
        "end": r_local1,
        "region_start_local": f_local0,
        "region_end_local": r_local1,
        "forward_start": ref_offset + f_local0 + 1,
        "forward_end": ref_offset + f_local0 + len(forward),
        "reverse_start": ref_offset + r_local1 - len(reverse) + 1,
        "reverse_end": ref_offset + r_local1,
        "probe_start": ref_offset + p_local0 + 1 if probe else None,
        "probe_end": ref_offset + p_local0 + len(probe) if probe else None,
        "probe_offset": p_local0 - f_local0 if probe else None,
        "source": "turkey_consensus",
        "created": datetime.now().isoformat(timespec="seconds"),
    })
    new["best"] = {
        "forward": forward,
        "reverse": reverse,
        "forward_start": new["forward_start"],
        "reverse_end": new["reverse_end"],
        "forward_len": len(forward),
        "reverse_len": len(reverse),
        "amplicon_size": len(amp_seq),
        "amplicon_seq": amp_seq,
        "forward_tm": candidate.get("forward_tm"),
        "reverse_tm": candidate.get("reverse_tm"),
        "forward_gc": candidate.get("forward_gc"),
        "reverse_gc": candidate.get("reverse_gc"),
        "score": candidate.get("score"),
    }
    if probe:
        new["probe_info"] = {
            "sequence": probe,
            "position": p_local0 - f_local0 + 1,
            "length": len(probe),
            "tm": candidate.get("probe_tm"),
            "gc": candidate.get("probe_gc"),
            "dimer_risk": False,
        }
    else:
        new["probe_info"] = None
    new["strict"] = bool(candidate.get("strict", True))
    new["cons_info"] = None
    new["cross_risk"] = None
    return new


def _check(key: str, label: str, status: str, detail: str) -> dict:
    return {"key": key, "label": label, "status": status, "detail": detail}


def build_order_gate(active: dict, turkey_result: dict | None = None,
                     specificity: dict | None = None,
                     primer_blast_done: bool = False,
                     synthesis_confirmed: bool = False,
                     wetlab_planned: bool = False) -> dict:
    """Sipariş öncesi kontrolleri tek karar tablosunda toplar.

    Not: PASS yalnızca in-silico kontrollerin geçtiğini gösterir; klinik validasyon değildir.
    """
    checks: list[dict] = []
    if not active:
        return {
            "overall": BLOCKED,
            "checks": [_check("active", "Aktif assay", BLOCKED, "Önce primer/prob tasarlayın.")],
            "blocking": ["Aktif assay yok."],
        }

    forward = str(active.get("forward") or "").upper()
    reverse = str(active.get("reverse") or "").upper()
    probe = str(active.get("probe") or "").upper()
    sequences = [forward, reverse] + ([probe] if probe else [])

    if all(s and re.fullmatch(r"[ACGT]+", s) for s in sequences):
        checks.append(_check("sequence", "Dizi bütünlüğü", PASS,
                             "Yalnızca A/C/G/T; N veya belirsiz baz yok."))
    else:
        checks.append(_check("sequence", "Dizi bütünlüğü", FAIL,
                             "Dizide N veya geçersiz karakter var."))

    size = int(active.get("amplicon_size") or 0)
    if 70 <= size <= 200:
        checks.append(_check("amplicon", "Amplikon boyu", PASS, f"{size} bp (TaqMan için uygun)."))
    elif 60 <= size <= 250:
        checks.append(_check("amplicon", "Amplikon boyu", WARN, f"{size} bp; hedef aralık 70–200 bp."))
    else:
        checks.append(_check("amplicon", "Amplikon boyu", FAIL, f"{size} bp; hedef aralık 70–200 bp."))

    try:
        f_info = analyze_oligo(forward, role="primer")
        r_info = analyze_oligo(reverse, role="primer")
        p_info = analyze_oligo(probe, role="probe") if probe else None
        _, dimer_flags = dimer_matrix({
            "F": forward, "R": reverse, **({"P": probe} if probe else {}),
        })
    except Exception as exc:
        f_info = r_info = p_info = None
        dimer_flags = []
        checks.append(_check("thermo", "Termodinamik", FAIL, f"Hesaplanamadı: {exc}"))

    if f_info and r_info:
        tms = [f_info.get("tm"), r_info.get("tm")]
        if all(t is not None and 58 <= t <= 65 for t in tms) and abs(tms[0] - tms[1]) <= 2:
            checks.append(_check("primer_tm", "Primer Tm", PASS,
                                 f"F {tms[0]} °C · R {tms[1]} °C · fark {abs(tms[0]-tms[1]):.1f} °C."))
        elif all(t is not None and 56 <= t <= 67 for t in tms) and abs(tms[0] - tms[1]) <= 3:
            checks.append(_check("primer_tm", "Primer Tm", WARN,
                                 f"F {tms[0]} °C · R {tms[1]} °C; 58–65 °C / ≤2 °C hedefi dışında."))
        else:
            checks.append(_check("primer_tm", "Primer Tm", FAIL,
                                 f"F {tms[0]} °C · R {tms[1]} °C; ciddi Tm uyumsuzluğu."))

        if p_info and p_info.get("tm") is not None:
            gap = p_info["tm"] - max(tms)
            if 3 <= gap <= 12 and 66 <= p_info["tm"] <= 73:
                checks.append(_check("probe_tm", "Prob Tm / fark", PASS,
                                     f"Prob {p_info['tm']} °C · fark {gap:+.1f} °C."))
            elif 0 <= gap <= 15:
                checks.append(_check("probe_tm", "Prob Tm / fark", WARN,
                                     f"Prob {p_info['tm']} °C · fark {gap:+.1f} °C; hedef 66–73 °C ve +3…+12 °C."))
            else:
                checks.append(_check("probe_tm", "Prob Tm / fark", FAIL,
                                     f"Prob {p_info['tm']} °C · fark {gap:+.1f} °C."))
        else:
            checks.append(_check("probe_tm", "Prob Tm / fark", WARN,
                                 "Prob yok; bu bir TaqMan assay değil."))

        structure_warn = []
        for info, name in ((f_info, "F"), (r_info, "R"), (p_info, "P")):
            if not info:
                continue
            if info.get("hairpin_dg", 0) <= -3.0:
                structure_warn.append(f"{name} hairpin {info['hairpin_dg']} kcal/mol")
            if info.get("three_prime_dg", 0) <= -6.0:
                structure_warn.append(f"{name} 3' uç stabilitesi {info['three_prime_dg']} kcal/mol")
            if info.get("max_run", 0) >= 5:
                structure_warn.append(f"{name} homopolymer {info['max_run']}")
        if dimer_flags:
            checks.append(_check("dimer", "Dimer / yapı", FAIL,
                                 "; ".join(f"{x['oligo_1']}-{x['oligo_2']} ΔG {x['dg']}" for x in dimer_flags[:3])))
        elif structure_warn:
            checks.append(_check("dimer", "Dimer / yapı", WARN, "; ".join(structure_warn[:4])))
        else:
            checks.append(_check("dimer", "Dimer / yapı", PASS, "Güçlü dimer/hairpin uyarısı yok."))

    if not probe:
        checks.append(_check("probe_mod", "Prob modifikasyonu", WARN, "Prob yok; raportör/quencher seçilemez."))
    elif probe.startswith("G"):
        checks.append(_check("probe_mod", "Prob modifikasyonu", WARN, "Prob 5' ucu G; raportör sondurmesi riski."))
    else:
        checks.append(_check("probe_mod", "Prob modifikasyonu", PASS, "Prob 5' ucu G değil."))

    if specificity:
        n_off = int(specificity.get("n_offtargets") or 0)
        score = float(specificity.get("score") or 0)
        if n_off == 0 and score >= 90:
            checks.append(_check("specificity", "Özgüllük taraması", PASS,
                                 f"Skor {score:.0f}/100 · off-target yok."))
        elif n_off == 0 and score >= 75:
            checks.append(_check("specificity", "Özgüllük taraması", WARN,
                                 f"Skor {score:.0f}/100; düşük riskli off-target yok."))
        else:
            checks.append(_check("specificity", "Özgüllük taraması", FAIL,
                                 f"Skor {score:.0f}/100 · {n_off} off-target."))
    else:
        checks.append(_check("specificity", "Özgüllük taraması", WARN,
                             "Henüz çalıştırılmadı; Analiz → Özgüllük sekmesini kullanın."))

    if turkey_result and turkey_result.get("n_total"):
        n_fail = int(turkey_result.get("n_fail") or 0)
        n_warn = int(turkey_result.get("n_warn") or 0)
        n_pass = int(turkey_result.get("n_pass") or 0)
        detail = f"{n_pass} PASS · {n_warn} WARN · {n_fail} FAIL / {turkey_result['n_total']} dizi."
        if n_fail:
            checks.append(_check("turkey", "Türkiye suşları", FAIL, detail))
        elif n_warn:
            checks.append(_check("turkey", "Türkiye suşları", WARN, detail))
        else:
            checks.append(_check("turkey", "Türkiye suşları", PASS, detail))
    else:
        checks.append(_check("turkey", "Türkiye suşları", WARN,
                             "Türkiye dizisi değerlendirilmedi; NCBI veya FASTA ile kontrol edin."))

    checks.append(_check("blast", "Primer-BLAST / BLASTn", PASS if primer_blast_done else WARN,
                         "Yapıldı olarak işaretlendi." if primer_blast_done else "Yapılmadı; önerilir."))
    checks.append(_check("synthesis", "Sentez/modifikasyon teyidi", PASS if synthesis_confirmed else WARN,
                         "Teyit edildi." if synthesis_confirmed else "Raportör, quencher, ölçek ve saflaştırma teyit edilmedi."))
    checks.append(_check("wetlab", "Islak laboratuvar planı", PASS if wetlab_planned else WARN,
                         "Planlandı." if wetlab_planned else "LOD/LOQ ve validasyon planı teyit edilmedi."))

    statuses = [c["status"] for c in checks]
    if FAIL in statuses:
        overall = FAIL
    elif WARN in statuses:
        overall = WARN
    else:
        overall = PASS
    blocking = [f"{c['label']}: {c['detail']}" for c in checks if c["status"] == FAIL]
    return {
        "overall": overall,
        "checks": checks,
        "blocking": blocking,
        "n_pass": statuses.count(PASS),
        "n_warn": statuses.count(WARN),
        "n_fail": statuses.count(FAIL),
    }


def order_sheet_rows(active: dict, reporter: str = "FAM", quencher: str = "BHQ1",
                     scale: str = "0.1 µmol", purification: str = "HPLC (prob) / standart (primer)") -> list[dict]:
    """Sipariş listesi satırları."""
    rows = []
    for name, key, role in (("Forward", "forward", "primer"),
                            ("Reverse", "reverse", "primer"),
                            ("Probe", "probe", "probe")):
        seq = str(active.get(key) or "").upper()
        if not seq:
            continue
        info = analyze_oligo(seq, role=role)
        rows.append({
            "Oligo": name,
            "Dizi (5'→3')": seq,
            "Uzunluk": len(seq),
            "Tm (°C)": info.get("tm"),
            "GC (%)": info.get("gc"),
            "5' modifikasyon": reporter if role == "probe" else "—",
            "3' modifikasyon": quencher if role == "probe" else "—",
            "Sentez ölçeği": scale,
            "Saflaştırma": purification,
            "Not": "Prob için HPLC önerilir." if role == "probe" else "Standart desalt yeterli olabilir.",
        })
    return rows
