"""
literature.py
Literatur / yayinlanmis assay karsilastirma ve assay kutuphanesi.

Amac: kullanicinin yapistirdigi yayin primer/prob setini referans genom uzerinde
konumlandirmak, uygulamanin aktif assay'i ile örtüsme/benzerlik/uyum skorunu
hesaplamak. Kutuphane baslangicta bos gelir; yanlis/dogrulanmamis sekans
uydurmamak icin ornek sekanslar gomulmez.
"""

from __future__ import annotations

import json
from datetime import datetime

import pandas as pd
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
