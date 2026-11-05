"""
thermo.py
Primer3 tabanli termodinamik analiz: Tm, hairpin, homodimer/heterodimer ΔG,
3' uc stabilitesi ve Tm belirsizligi (Monte Carlo).

DG degerleri kcal/mol'e cevrilir (primer3 cal/mol dondurur).
"""

from __future__ import annotations

import random
import statistics

import primer3

# qPCR kosullari (pcr_engine.QPCR_PARAMS ile uyumlu)
QPCR = dict(mv_conc=50.0, dv_conc=1.5, dntp_conc=0.2, dna_conc=250.0)
QPCR_STRUCT = {**QPCR, "temp_c": 60.0}   # annealing sicakliginda yapi analizi

DIMER_WARN_DG = -6.0          # kcal/mol, 60 °C'de heterodimer/homodimer
HAIRPIN_WARN_DG = -3.0        # kcal/mol, 60 °C'de hairpin
THREE_PRIME_WARN_DG = -6.0    # kcal/mol, 37 °C'de son 5 baz


def _kcal(dg):
    try:
        return dg / 1000.0
    except Exception:
        return 0.0


def calc_tm(seq, **overrides):
    """Primer3 SantaLucia Tm (qPCR tuz kosullari)."""
    params = {**QPCR, **overrides}
    try:
        return round(float(primer3.calc_tm(str(seq).upper(), **params)), 1)
    except Exception:
        return None


def _result(seq1, seq2=None, temp_c=60.0, is_hairpin=False):
    params = {**QPCR_STRUCT, "temp_c": temp_c}
    try:
        if is_hairpin:
            r = primer3.calc_hairpin(str(seq1).upper(), **params)
        elif seq2 is None:
            r = primer3.calc_homodimer(str(seq1).upper(), **params)
        else:
            r = primer3.calc_heterodimer(str(seq1).upper(), str(seq2).upper(), **params)
    except Exception:
        return {"tm": None, "dg": 0.0, "dh": None, "ds": None, "structure_found": False}
    if not r.structure_found or r.dg >= 0:
        return {"tm": None, "dg": 0.0, "dh": None, "ds": None, "structure_found": False}
    return {
        "tm": round(float(r.tm), 1),
        "dg": round(_kcal(r.dg), 2),
        "dh": round(_kcal(r.dh), 1),
        "ds": round(float(r.ds), 2),
        "structure_found": True,
    }


def hairpin(seq):
    return _result(seq, is_hairpin=True)


def homodimer(seq):
    return _result(seq)


def heterodimer(seq1, seq2):
    return _result(seq1, seq2)


def three_prime_stability(seq, n=5):
    """Son n bazin kendi tamamlayicisiyla dupleks DG'si (37 °C, goreli 3' stabilite)."""
    seq = str(seq).upper()
    tail = seq[-n:] if len(seq) >= n else seq
    return _result(tail, temp_c=37.0)


def gc_content(seq):
    seq = str(seq).upper()
    if not seq:
        return 0.0
    return round(100.0 * (seq.count("G") + seq.count("C")) / len(seq), 1)


def max_run(seq):
    seq = str(seq).upper()
    best = run = 1 if seq else 0
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def analyze_oligo(seq, role="primer"):
    """Tek oligo icin tum termodinamik metrikler + uyarilar."""
    seq = str(seq).upper().strip()
    tm = calc_tm(seq)
    gc = gc_content(seq)
    hp = hairpin(seq)
    hd = homodimer(seq)
    tp = three_prime_stability(seq)

    warnings = []
    if "N" in seq:
        warnings.append("Dizide N var; sentez/analiz icin uygun degil.")
    if hp["dg"] <= HAIRPIN_WARN_DG:
        warnings.append(f"Guclu hairpin (ΔG {hp['dg']} kcal/mol).")
    if hd["dg"] <= DIMER_WARN_DG:
        warnings.append(f"Guclu homodimer (ΔG {hd['dg']} kcal/mol).")
    if tp["dg"] <= THREE_PRIME_WARN_DG:
        warnings.append(f"3' uc stabilitesi yuksek (ΔG {tp['dg']} kcal/mol).")
    if max_run(seq) >= 5:
        warnings.append("5+ ayni baz dizisi (homopolymer).")
    if role == "probe" and seq[:1] == "G":
        warnings.append("Prob 5' ucu G; raportor sondurmesi riski.")
    if role == "primer":
        if tm is not None and not (58.0 <= tm <= 65.0):
            warnings.append("Primer Tm 58-65 °C araligi disinda.")
        if not (40.0 <= gc <= 60.0):
            warnings.append("Primer GC %40-60 araligi disinda.")
        if seq[-1:] not in ("G", "C"):
            warnings.append("3' uc G/C degil (uzama verimi dusebilir).")
    if role == "probe":
        if tm is not None and not (66.0 <= tm <= 73.0):
            warnings.append("Prob Tm 66-73 °C araligi disinda.")
        if not (40.0 <= gc <= 65.0):
            warnings.append("Prob GC %40-65 araligi disinda.")

    return {
        "sequence": seq,
        "length": len(seq),
        "tm": tm,
        "gc": gc,
        "hairpin_dg": hp["dg"],
        "hairpin_tm": hp["tm"],
        "homodimer_dg": hd["dg"],
        "homodimer_tm": hd["tm"],
        "three_prime_dg": tp["dg"],
        "three_prime_tm": tp["tm"],
        "gc_clamp": sum(1 for b in seq[-5:] if b in "GC"),
        "max_run": max_run(seq),
        "n_count": seq.count("N"),
        "five_prime_base": seq[:1],
        "three_prime_base": seq[-1:],
        "warnings": warnings,
    }


def dimer_matrix(oligos):
    """
    oligos: {"F": "...", "R": "...", "P": "..."}
    Donus: (pandas.DataFrame ΔG matrisi, uyari listesi)
    """
    import pandas as pd

    names = list(oligos)
    rows = []
    flags = []
    for a in names:
        row = {"Oligo": a}
        for b in names:
            if a == b:
                r = homodimer(oligos[a])
                kind = "homodimer"
            else:
                r = heterodimer(oligos[a], oligos[b])
                kind = "heterodimer"
            row[b] = r["dg"]
            if r["dg"] <= DIMER_WARN_DG:
                flags.append({
                    "oligo_1": a,
                    "oligo_2": b,
                    "tip": kind,
                    "dg": r["dg"],
                    "tm": r["tm"],
                })
        rows.append(row)
    df = pd.DataFrame(rows).set_index("Oligo")
    # simetrik matriste tekrar eden uyarilari tekille
    seen = set()
    unique = []
    for f in flags:
        key = tuple(sorted([f["oligo_1"], f["oligo_2"]])) + (f["tip"],)
        if key not in seen:
            seen.add(key)
            unique.append(f)
    return df, unique


def tm_uncertainty(seq, n=1000, seed=42):
    """Tuz/Mg/dNTP/oligo konsantrasyonu araliklarinda Tm dagilimi."""
    seq = str(seq).upper()
    rng = random.Random(seed)
    vals = []
    for _ in range(max(50, int(n))):
        try:
            vals.append(float(primer3.calc_tm(
                seq,
                mv_conc=rng.uniform(40.0, 60.0),
                dv_conc=rng.uniform(1.0, 3.0),
                dntp_conc=rng.uniform(0.1, 0.3),
                dna_conc=rng.uniform(200.0, 300.0),
            )))
        except Exception:
            continue
    if not vals:
        return None
    vals.sort()
    mean = statistics.fmean(vals)
    sd = statistics.pstdev(vals) if len(vals) > 1 else 0.0
    lo = vals[int(0.025 * (len(vals) - 1))]
    hi = vals[int(0.975 * (len(vals) - 1))]
    return {
        "n": len(vals),
        "mean": round(mean, 1),
        "sd": round(sd, 2),
        "p2.5": round(lo, 1),
        "p97.5": round(hi, 1),
        "min": round(vals[0], 1),
        "max": round(vals[-1], 1),
    }


def thermo_score(analysis, role="primer"):
    """0-100 arasi aciklanabilir termodinamik skor."""
    score = 100.0
    tm = analysis.get("tm")
    gc = analysis.get("gc", 0.0)
    if role == "primer":
        if tm is not None and not (58.0 <= tm <= 65.0):
            score -= 20
        if not (40.0 <= gc <= 60.0):
            score -= 15
        if analysis.get("three_prime_base") not in ("G", "C"):
            score -= 5
    else:
        if tm is not None and not (66.0 <= tm <= 73.0):
            score -= 20
        if not (40.0 <= gc <= 65.0):
            score -= 15
        if analysis.get("five_prime_base") == "G":
            score -= 10
    if analysis.get("hairpin_dg", 0) <= HAIRPIN_WARN_DG:
        score -= 15
    if analysis.get("homodimer_dg", 0) <= DIMER_WARN_DG:
        score -= 20
    if analysis.get("three_prime_dg", 0) <= THREE_PRIME_WARN_DG:
        score -= 10
    if analysis.get("max_run", 0) >= 5:
        score -= 10
    if analysis.get("n_count", 0):
        score -= 20
    return round(max(0.0, min(100.0, score)), 1)
