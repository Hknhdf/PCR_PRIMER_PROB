"""
advanced_design.py
Ileri tasarim modulleri:
  - consensus / dejenere primer tasarimi (coklu hizalama uzerinden)
  - SNP / alel ayrimci primer tasarimi (3' uc SNP)
  - RPA primer + prob onerisi (beta)
  - LAMP F3/B3/FIP/BIP tasarimi (beta) ve LAMP primer kontrolu

Not: RPA ve LAMP modulleri laboratuvar optimizasyonu gerektiren sezgisel
tasarimlardir; sonuclar in-silico kontrol edilir, wet-lab dogrulama yerine gecmez.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import numpy as np
import pandas as pd
from Bio import Align

from pcr_engine import (
    calculate_gc_content,
    calculate_tm,
    design_primers_best,
    design_probe,
    dimer_penalty,
    dimer_risk,
    has_bad_run,
    revcomp,
)
from thermo import heterodimer, homodimer as thermo_homodimer


IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "AC": "M", "AG": "R", "AT": "W", "CG": "S", "CT": "Y", "GT": "K",
    "ACG": "V", "ACT": "H", "AGT": "D", "CGT": "B", "ACGT": "N",
}
IUPAC_COUNT = {code: len(bases) for bases, code in IUPAC.items()}
IUPAC_COUNT.update({"N": 4})


def iupac_code(bases: str):
    key = "".join(sorted(set(b for b in bases.upper() if b in "ACGT")))
    return IUPAC.get(key, "N")


def degeneracy(seq: str):
    n = 1
    for b in str(seq).upper():
        n *= IUPAC_COUNT.get(b, 4)
    return n


def _anchor_offset(ref: str, subject: str, k: int = 20):
    L = len(ref)
    if L < k:
        return None
    for frac in (0.5, 0.35, 0.65, 0.2, 0.8, 0.45, 0.6, 0.3, 0.1, 0.9):
        p = max(0, min(L - k, int(L * frac) - k // 2))
        probe = ref[p:p + k]
        idx = subject.find(probe)
        if idx >= 0:
            return idx - p
    return None


def align_to_reference(ref_region: str, subject: str):
    """Subject'i referans bolgeye global hizalar; referans pozisyonlari icin baz listesi."""
    ref_region = str(ref_region).upper()
    subject = str(subject).upper()
    if not ref_region or not subject:
        return None, None
    if ref_region in subject:
        sub = subject[subject.find(ref_region):subject.find(ref_region) + len(ref_region)]
    else:
        off = _anchor_offset(ref_region, subject)
        if off is None:
            return None, None
        lo = max(0, off - 50)
        hi = min(len(subject), off + len(ref_region) + 50)
        sub = subject[lo:hi]
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    try:
        aln = aligner.align(ref_region, sub)[0]
        ref_aln, sub_aln = str(aln[0]), str(aln[1])
    except Exception:
        return None, None
    bases = [None] * len(ref_region)
    ri = matches = compared = 0
    for a, b in zip(ref_aln, sub_aln):
        if a == "-":
            continue
        if ri >= len(bases):
            break
        bases[ri] = b
        if b != "-":
            compared += 1
            if a == b:
                matches += 1
        ri += 1
    identity = matches / compared if compared else None
    return bases, identity


def build_consensus(ref_region: str, subjects):
    """
    subjects: [{'id','description','seq'}, ...] veya [(id, desc, seq)]
    Donus: majority (ACGT), degenerate (IUPAC), depth, identity, rows
    """
    ref_region = str(ref_region).upper()
    rows = []
    meta = []
    for item in subjects:
        if isinstance(item, dict):
            sid, desc, seq = item.get("id", "?"), item.get("description", ""), item.get("seq", "")
        else:
            sid, desc, seq = item[0], item[1], item[2]
        bases, identity = align_to_reference(ref_region, seq)
        if bases is None:
            continue
        rows.append(bases)
        meta.append({"id": sid, "description": desc, "identity": identity})
    if not rows:
        return {"ok": False, "message": "Diziler referans bölgeye hizalanamadı."}

    majority = []
    degenerate = []
    depth = []
    frequencies = []
    for i in range(len(ref_region)):
        counts = {b: 0 for b in "ACGT"}
        for row in rows:
            b = row[i] if i < len(row) else None
            if b in counts:
                counts[b] += 1
        n = sum(counts.values())
        depth.append(n)
        frequencies.append(counts)
        if n == 0:
            majority.append("N")
            degenerate.append("N")
            continue
        top = max(counts, key=counts.get)
        majority.append(top)
        present = [b for b in "ACGT" if counts[b] > 0]
        if counts[top] / n >= 0.85:
            degenerate.append(top)
        else:
            # en az %10 frekansa sahip bazlari IUPAC'a dahil et
            incl = [b for b in present if counts[b] / n >= 0.10]
            degenerate.append(iupac_code("".join(incl)) if incl else top)
    return {
        "ok": True,
        "majority": "".join(majority),
        "degenerate": "".join(degenerate),
        "depth": depth,
        "frequencies": frequencies,
        "rows": rows,
        "meta": meta,
        "n": len(rows),
        "mean_identity": round(float(np.mean([m["identity"] for m in meta if m["identity"] is not None])), 4)
        if any(m["identity"] is not None for m in meta) else None,
    }


def _coverage_for_segment(rows, start, end, majority):
    """Satirlarin segment boyunca majority'ye <=1 mismatch orani."""
    if not rows:
        return 0.0
    ok = 0
    for row in rows:
        mm = 0
        for i in range(start, min(end, len(majority))):
            b = row[i] if i < len(row) else None
            if b is None or b == "-" or b != majority[i]:
                mm += 1
        if mm <= 1:
            ok += 1
    return ok / len(rows)


def design_consensus_primers(ref_region: str, subjects, region_start: int = 0,
                             target_amp: int = 150, top_n: int = 10,
                             min_coverage: float = 0.80, max_degeneracy: int = 256):
    """Consensus/dejenere primer + prob onerisi."""
    cons = build_consensus(ref_region, subjects)
    if not cons.get("ok"):
        return cons
    majority = cons["majority"]
    degenerate = cons["degenerate"]
    rows = cons["rows"]
    if "N" in majority:
        # tasarim icin N'leri en sik bazla degistir
        majority = "".join(
            max(cons["frequencies"][i], key=cons["frequencies"][i].get) if majority[i] == "N" else majority[i]
            for i in range(len(majority))
        )

    sets, strict = design_primers_best(
        majority, 0, len(majority), target_amp, top_n=max(20, top_n * 2), conservation=None
    )
    out = []
    for cand in sets:
        f_start = cand["forward_start"] - 1
        f_end = f_start + cand["forward_len"]
        r_start = cand["reverse_end"] - cand["reverse_len"]
        r_end = cand["reverse_end"]
        cov_f = _coverage_for_segment(rows, f_start, f_end, majority)
        cov_r = _coverage_for_segment(rows, r_start, r_end, majority)
        if cov_f < min_coverage or cov_r < min_coverage:
            continue
        deg_f = degeneracy(degenerate[f_start:f_end])
        deg_r = degeneracy(degenerate[r_start:r_end])
        if deg_f > max_degeneracy or deg_r > max_degeneracy:
            continue
        probe = design_probe(
            cand["amplicon_seq"], cand["forward_len"], cand["reverse_len"],
            avoid_sequences=(cand["forward"], cand["reverse"]),
        )
        if not probe:
            continue
        p_start = f_start + probe["position"] - 1
        p_end = p_start + probe["length"]
        cov_p = _coverage_for_segment(rows, p_start, p_end, majority)
        deg_p = degeneracy(degenerate[p_start:p_end])
        if cov_p < min_coverage or deg_p > max_degeneracy:
            continue
        score = cand["score"] + 25 * (cov_f + cov_r + cov_p) / 3 - 3 * math.log2(max(1, deg_f * deg_r * deg_p))
        out.append({
            "forward": cand["forward"],
            "reverse": cand["reverse"],
            "probe": probe["sequence"],
            "forward_degenerate": degenerate[f_start:f_end],
            "reverse_degenerate": revcomp(degenerate[r_start:r_end]),
            "reverse_site_degenerate": degenerate[r_start:r_end],
            "probe_degenerate": degenerate[p_start:p_end],
            "forward_deg": deg_f,
            "reverse_deg": deg_r,
            "probe_deg": deg_p,
            "forward_tm": cand["forward_tm"],
            "reverse_tm": cand["reverse_tm"],
            "probe_tm": probe["tm"],
            "forward_gc": cand["forward_gc"],
            "reverse_gc": cand["reverse_gc"],
            "probe_gc": probe["gc"],
            "amplicon_size": cand["amplicon_size"],
            "forward_start": region_start + f_start + 1,
            "reverse_end": region_start + r_end,
            "probe_start": region_start + p_start + 1,
            "coverage_f": round(100 * cov_f, 1),
            "coverage_r": round(100 * cov_r, 1),
            "coverage_p": round(100 * cov_p, 1),
            "score": round(score, 1),
        })
    out.sort(key=lambda x: x["score"], reverse=True)
    return {
        "ok": True,
        "consensus": cons,
        "majority": majority,
        "degenerate": degenerate,
        "strict": strict,
        "candidates": out[:top_n],
    }


def design_snp_primers(ref_seq: str, region_start: int, region_end: int,
                       snp_pos: int, alt_base: str, reverse: str | None = None,
                       target_amp: int = 120, top_n: int = 12):
    """
    SNP pozisyonu 0-based genomik koordinat; primerin 3' ucu bu pozisyona denk gelir.
    Mutant (alt) primer ile WT kontrol primerini birlikte dondurur.
    """
    ref_seq = str(ref_seq).upper()
    alt = str(alt_base).upper().strip()
    if alt not in "ACGT":
        return {"ok": False, "message": "Alt baz A/C/G/T olmalı."}
    region_start = max(0, int(region_start))
    region_end = min(len(ref_seq), int(region_end))
    if not (region_start <= snp_pos < region_end):
        return {"ok": False, "message": "SNP pozisyonu hedef bölge içinde değil."}
    ref_base = ref_seq[snp_pos]
    if alt == ref_base:
        return {"ok": False, "message": "Alt baz referans baz ile aynı."}
    region = ref_seq[region_start:region_end]
    snp_rel = snp_pos - region_start

    reverse_candidates = []
    if reverse:
        reverse = str(reverse).upper().strip()
        r_pos = region.find(reverse)
        if r_pos < 0:
            r_pos = region.find(revcomp(reverse))
        if r_pos >= 0 and r_pos > snp_rel:
            reverse_candidates.append({"reverse": reverse, "start": r_pos, "end": r_pos + len(reverse)})
    if not reverse_candidates:
        for end in range(max(snp_rel + 60, snp_rel + target_amp - 45),
                         min(len(region), snp_rel + target_amp + 45) + 1):
            for r_len in range(18, 26):
                start = end - r_len
                if start <= snp_rel + 20:
                    continue
                site = region[start:end]
                if "N" in site or has_bad_run(site, 5):
                    continue
                r_seq = revcomp(site)
                r_tm = calculate_tm(r_seq)
                r_gc = calculate_gc_content(r_seq)
                if r_tm is None or not (58 <= r_tm <= 66) or not (35 <= r_gc <= 65):
                    continue
                reverse_candidates.append({"reverse": r_seq, "start": start, "end": end,
                                           "tm": r_tm, "gc": r_gc})
        reverse_candidates.sort(key=lambda x: (abs(x.get("tm", 61) - 61), abs(x["end"] - snp_rel - target_amp)))
        reverse_candidates = reverse_candidates[:8]

    out = []
    for f_len in range(18, 26):
        start = snp_rel - f_len + 1
        if start < 0:
            continue
        wt_forward = region[start:snp_rel] + ref_base
        mut_forward = region[start:snp_rel] + alt
        if "N" in mut_forward or has_bad_run(mut_forward, 5):
            continue
        f_tm = calculate_tm(mut_forward)
        f_gc = calculate_gc_content(mut_forward)
        if f_tm is None or not (56 <= f_tm <= 68) or not (30 <= f_gc <= 70):
            continue
        target_wt = region[start:snp_rel + 1]
        target_alt = region[start:snp_rel] + alt
        wt_anti = revcomp(target_wt)
        alt_anti = revcomp(target_alt)
        perfect = heterodimer(mut_forward, alt_anti)
        mismatch = heterodimer(mut_forward, wt_anti)
        if not perfect["structure_found"] and not mismatch["structure_found"]:
            perfect = heterodimer(mut_forward[-8:], alt_anti[-8:])
            mismatch = heterodimer(mut_forward[-8:], wt_anti[-8:])
        delta_dg = round(mismatch["dg"] - perfect["dg"], 2)
        # Forward primer sense ile ayni yonde; 3' uc karsilastirmasi sense hedefle yapilir.
        mm_wt = sum(1 for a, b in zip(mut_forward[-5:], target_wt[-5:]) if a != b)
        mm_alt = sum(1 for a, b in zip(mut_forward[-5:], target_alt[-5:]) if a != b)
        terminal_wt = mut_forward[-1] != target_wt[-1]
        terminal_alt = mut_forward[-1] != target_alt[-1]
        wt_tm = calculate_tm(wt_forward)
        for rc in reverse_candidates:
            r_seq = rc["reverse"]
            amp_size = rc["end"] - start
            if not (60 <= amp_size <= 300):
                continue
            if thermo_homodimer(mut_forward)["dg"] <= -6.0:
                continue
            if thermo_homodimer(r_seq)["dg"] <= -6.0:
                continue
            if heterodimer(mut_forward, r_seq)["dg"] <= -6.0:
                continue
            if dimer_risk(mut_forward[-6:], r_seq[-6:], min_overlap=5)[0]:
                continue
            r_tm = rc.get("tm") or calculate_tm(r_seq)
            r_gc = rc.get("gc") or calculate_gc_content(r_seq)
            if terminal_wt and not terminal_alt:
                disc = "güçlü (3' uç mismatch)"
            elif mm_wt > mm_alt:
                disc = "orta"
            else:
                disc = "zayıf"
            out.append({
                "forward": mut_forward,
                "forward_wt": wt_forward,
                "reverse": r_seq,
                "forward_tm": f_tm,
                "reverse_tm": r_tm,
                "forward_gc": f_gc,
                "reverse_gc": r_gc,
                "wt_forward_tm": wt_tm,
                "tm_delta_wt": round((f_tm or 0) - (wt_tm or 0), 1),
                "delta_dg": delta_dg,
                "wt_last5_mismatch": mm_wt,
                "alt_last5_mismatch": mm_alt,
                "terminal_mismatch_wt": bool(terminal_wt),
                "terminal_mismatch_alt": bool(terminal_alt),
                "amplicon_size": amp_size,
                "snp_pos": snp_pos + 1,
                "ref_base": ref_base,
                "alt_base": alt,
                "forward_start": region_start + start + 1,
                "reverse_end": region_start + rc["end"],
                "discrimination": disc,
            })
    # her forward icin en iyi reverse'i tut
    best = {}
    for row in out:
        key = row["forward"]
        if key not in best or (row["wt_last5_mismatch"], -row["alt_last5_mismatch"]) > (
                best[key]["wt_last5_mismatch"], -best[key]["alt_last5_mismatch"]):
            best[key] = row
    rank = {"güçlü (3' uç mismatch)": 0, "orta": 1, "zayıf": 2}
    rows = sorted(best.values(), key=lambda x: (
        rank.get(x["discrimination"], 3), -x["wt_last5_mismatch"], abs((x["forward_tm"] or 0) - 61)
    ))
    return {"ok": True, "ref_base": ref_base, "alt_base": alt, "snp_pos": snp_pos + 1,
            "candidates": rows[:top_n]}


# ---------------------------------------------------------------------------
# RPA (beta)
# ---------------------------------------------------------------------------

def design_rpa(ref_region: str, region_start: int = 0, target_amp: int = 200, top_n: int = 6):
    """30-35 nt RPA primer cifti + 46-52 nt prob onerisi (beta)."""
    region = str(ref_region).upper()
    if len(region) < 120:
        return {"ok": False, "message": "RPA için bölge en az 120 bp olmalı."}
    out = []
    for amp_len in range(max(100, target_amp - 50), min(len(region), target_amp + 50) + 1, 5):
        for offset in range(0, len(region) - amp_len + 1, 8):
            window = region[offset:offset + amp_len]
            for f_len in range(30, 36):
                forward = window[:f_len]
                if "N" in forward or has_bad_run(forward, 5):
                    continue
                f_gc = calculate_gc_content(forward)
                if not (30 <= f_gc <= 70):
                    continue
                f_tm = calculate_tm(forward)
                for r_len in range(30, 36):
                    reverse_site = window[-r_len:]
                    if "N" in reverse_site or has_bad_run(reverse_site, 5):
                        continue
                    reverse = revcomp(reverse_site)
                    r_gc = calculate_gc_content(reverse)
                    if not (30 <= r_gc <= 70):
                        continue
                    if dimer_risk(forward, reverse)[0] or dimer_risk(forward)[0] or dimer_risk(reverse)[0]:
                        continue
                    r_tm = calculate_tm(reverse)
                    score = 100 - abs(f_gc - 50) - abs(r_gc - 50) - abs(amp_len - target_amp) * 0.5
                    score -= dimer_penalty(forward, reverse) * 0.3
                    out.append({
                        "forward": forward,
                        "reverse": reverse,
                        "forward_tm": f_tm,
                        "reverse_tm": r_tm,
                        "forward_gc": f_gc,
                        "reverse_gc": r_gc,
                        "amplicon_size": amp_len,
                        "forward_start": region_start + offset + 1,
                        "reverse_end": region_start + offset + amp_len,
                        "score": round(score, 1),
                    })
    out.sort(key=lambda x: x["score"], reverse=True)
    unique, seen = [], set()
    for row in out:
        key = (row["forward"], row["reverse"])
        if key in seen:
            continue
        seen.add(key)
        # prob onerisi
        amp = region[row["forward_start"] - region_start - 1: row["reverse_end"] - region_start]
        probe = _rpa_probe(amp, len(row["forward"]), len(row["reverse"]))
        row["probe"] = probe
        unique.append(row)
        if len(unique) >= top_n:
            break
    return {"ok": True, "candidates": unique}


def _rpa_probe(amplicon: str, f_len: int, r_len: int):
    inner = amplicon[f_len:len(amplicon) - r_len]
    best = None
    for p_len in range(46, 53):
        for i in range(0, max(0, len(inner) - p_len + 1), 1):
            probe = inner[i:i + p_len]
            if "N" in probe or has_bad_run(probe, 5) or probe[0] == "G":
                continue
            gc = calculate_gc_content(probe)
            if not (30 <= gc <= 70):
                continue
            thf = 28
            if probe[thf] in "ACGT" and not has_bad_run(probe[thf - 2:thf + 3], 4):
                score = -abs(gc - 50) - abs(len(probe) - 49) * 0.5
                if best is None or score > best["score"]:
                    best = {
                        "sequence": probe,
                        "length": len(probe),
                        "gc": gc,
                        "tm": calculate_tm(probe),
                        "thf_position": thf + 1,
                        "label_5": "FAM / Texas Red",
                        "quencher_3": "BHQ",
                        "score": score,
                    }
    return best


# ---------------------------------------------------------------------------
# LAMP (beta)
# ---------------------------------------------------------------------------

def _candidate_primers(region: str, start: int, end: int, lengths=range(18, 26), top=4,
                       max_end: int | None = None):
    """Belirli bir pencerede Tm/GC'ye gore en iyi oligolari dondurur."""
    out = []
    for L in lengths:
        for pos in range(start, min(end, len(region) - L) + 1):
            if max_end is not None and pos + L > max_end:
                continue
            seq = region[pos:pos + L]
            if "N" in seq or has_bad_run(seq, 5):
                continue
            tm = calculate_tm(seq)
            gc = calculate_gc_content(seq)
            if tm is None or not (58 <= tm <= 68) or not (35 <= gc <= 65):
                continue
            score = -abs(tm - 62) * 2 - abs(gc - 50)
            out.append({"seq": seq, "pos": pos, "len": L, "tm": tm, "gc": gc, "score": score})
    out.sort(key=lambda x: x["score"], reverse=True)
    # ayni pozisyondan tek aday
    unique, seen = [], set()
    for item in out:
        if item["pos"] in seen:
            continue
        seen.add(item["pos"])
        unique.append(item)
        if len(unique) >= top:
            break
    return unique


def _lamp_dimer_ok(a: str, b: str):
    """LAMP icin 3' uc komplementerligi + 60 °C heterodimer DG kontrolu."""
    if dimer_risk(a[-6:], b[-6:], min_overlap=5)[0]:
        return False
    return heterodimer(a, b)["dg"] > -8.0


def design_lamp(ref_region: str, region_start: int = 0, target_amp: int = 220, top_n: int = 5):
    """Basit F3/B3/FIP/BIP tasarimi (beta; loop primerleri yok)."""
    region = str(ref_region).upper()
    if len(region) < 160:
        return {"ok": False, "message": "LAMP için bölge en az 160 bp olmalı."}
    sets = []
    for win_len in range(max(160, target_amp - 40), min(len(region), target_amp + 40) + 1, 20):
        for w0 in range(0, len(region) - win_len + 1, 10):
            w = region[w0:w0 + win_len]
            f3_c = _candidate_primers(w, 0, 12, top=3)
            if not f3_c:
                continue
            f2_c, f1_c, b3_c, b2_c, b1_c = [], [], [], [], []
            for f3 in f3_c:
                f2_c += _candidate_primers(w, f3["pos"] + f3["len"], f3["pos"] + f3["len"] + 15, top=2)
            for f2 in f2_c:
                f1_c += _candidate_primers(w, f2["pos"] + f2["len"] + 20, f2["pos"] + f2["len"] + 55, top=2)
            b3_c = _candidate_primers(w, win_len - 25, win_len - 15, top=3)
            for b3 in b3_c:
                b2_c += _candidate_primers(w, max(0, b3["pos"] - 60), max(0, b3["pos"] - 18),
                                           top=3, max_end=b3["pos"])
            for b2 in b2_c:
                b1_c += _candidate_primers(w, max(0, b2["pos"] - 60), max(0, b2["pos"] - 18),
                                           top=3, max_end=b2["pos"])

            # benzersizleri sinirla
            def dedup(items):
                out, seen = [], set()
                for x in sorted(items, key=lambda z: z["score"], reverse=True):
                    if x["pos"] in seen:
                        continue
                    seen.add(x["pos"])
                    out.append(x)
                return out[:2]

            for f3 in f3_c:
                for f2 in dedup([x for x in f2_c if x["pos"] >= f3["pos"] + f3["len"]]):
                    for f1 in dedup([x for x in f1_c if x["pos"] >= f2["pos"] + f2["len"] + 20]):
                        for b1 in dedup([x for x in b1_c if x["pos"] > f1["pos"] + f1["len"]]):
                            for b2 in dedup([x for x in b2_c if x["pos"] >= b1["pos"] + b1["len"] + 20]):
                                for b3 in dedup([x for x in b3_c if x["pos"] >= b2["pos"] + b2["len"]]):
                                    if b1["pos"] < f1["pos"] + f1["len"]:
                                        continue
                                    fip = f1["seq"] + revcomp(f2["seq"])
                                    bip = b1["seq"] + revcomp(b2["seq"])
                                    f3_seq = f3["seq"]
                                    b3_seq = revcomp(b3["seq"])
                                    oligos = [f3_seq, b3_seq, fip, bip]
                                    if any(has_bad_run(x, 5) for x in oligos):
                                        continue
                                    if any(thermo_homodimer(x)["dg"] <= -8.0 for x in oligos):
                                        continue
                                    if any(not _lamp_dimer_ok(a, b) for a, b in itertools.combinations(oligos, 2)):
                                        continue
                                    score = (f3["score"] + f2["score"] + f1["score"] + b1["score"] + b2["score"] + b3["score"])
                                    score -= abs(win_len - target_amp) * 0.4
                                    sets.append({
                                        "f3": f3_seq,
                                        "b3": b3_seq,
                                        "fip": fip,
                                        "bip": bip,
                                        "f3_tm": f3["tm"], "b3_tm": b3["tm"],
                                        "fip_tm": calculate_tm(fip), "bip_tm": calculate_tm(bip),
                                        "f3_gc": f3["gc"], "b3_gc": b3["gc"],
                                        "fip_gc": calculate_gc_content(fip), "bip_gc": calculate_gc_content(bip),
                                        "f2": f2["seq"], "f1c": f1["seq"],
                                        "b1c": b1["seq"], "b2": revcomp(b2["seq"]),
                                        "window_start": region_start + w0 + 1,
                                        "window_end": region_start + w0 + win_len,
                                        "f3_start": region_start + w0 + f3["pos"] + 1,
                                        "b3_start": region_start + w0 + b3["pos"] + 1,
                                        "score": round(score, 1),
                                    })
    sets.sort(key=lambda x: x["score"], reverse=True)
    unique, seen = [], set()
    for s in sets:
        key = (s["f3"], s["b3"], s["fip"], s["bip"])
        if key in seen:
            continue
        seen.add(key)
        unique.append(s)
        if len(unique) >= top_n:
            break
    return {"ok": True, "candidates": unique, "n_total": len(sets)}


def check_lamp_set(f3: str, b3: str, fip: str, bip: str, lf: str = "", lb: str = ""):
    """Yayinlanmis LAMP primer setini termodinamik/dimer acisindan kontrol eder."""
    seqs = {
        "F3": str(f3).upper().strip(),
        "B3": str(b3).upper().strip(),
        "FIP": str(fip).upper().strip(),
        "BIP": str(bip).upper().strip(),
    }
    if lf:
        seqs["LF"] = str(lf).upper().strip()
    if lb:
        seqs["LB"] = str(lb).upper().strip()
    rows = []
    for name, seq in seqs.items():
        rows.append({
            "Primer": name,
            "Dizi": seq,
            "Uzunluk": len(seq),
            "Tm (°C)": calculate_tm(seq),
            "GC (%)": calculate_gc_content(seq),
            "3' baz": seq[-1:] if seq else "",
            "Homopolymer": "VAR" if has_bad_run(seq, 5) else "yok",
        })
    conflicts = []
    for (n1, s1), (n2, s2) in itertools.combinations(seqs.items(), 2):
        risk, kmer = dimer_risk(s1, s2)
        if risk:
            conflicts.append({"Primer 1": n1, "Primer 2": n2, "5-mer": kmer, "Risk": "dimer"})
    return {"table": pd.DataFrame(rows), "conflicts": pd.DataFrame(conflicts)}
