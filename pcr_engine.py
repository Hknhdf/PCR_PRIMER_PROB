"""
pcr_engine.py
Real-Time PCR primer / TaqMan prob tasarim motoru (Streamlit'ten bagimsiz, test edilebilir).

Tm hesabi: SantaLucia nearest-neighbour (Biopython MeltingTemp.Tm_NN)
qPCR kosullari: 50 mM Na+, 1.5 mM Mg2+, 0.2 mM dNTP, 250 nM oligo
"""

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

QPCR_PARAMS = dict(Na=50, Mg=1.5, dNTPs=0.2, dnac1=250, dnac2=250, saltcorr=5)

# Varsayilan kabul kriterleri
PRIMER_TM_MIN, PRIMER_TM_MAX = 58.0, 65.0
PRIMER_GC_MIN, PRIMER_GC_MAX = 40.0, 60.0
PRIMER_TM_MAX_DIFF = 2.0
PROBE_TM_MIN, PROBE_TM_MAX = 66.0, 73.0
PROBE_GC_MIN, PROBE_GC_MAX = 40.0, 65.0


def calculate_tm(dna_sequence):
    """Nearest-neighbour Tm; basarisiz olursa GC bazli yaklasik formul."""
    try:
        sequence = str(dna_sequence).upper()
        if len(sequence) < 8:
            return None
        try:
            return round(float(mt.Tm_NN(sequence, **QPCR_PARAMS)), 1)
        except Exception:
            gc = sequence.count("G") + sequence.count("C")
            return round(64.9 + 41 * (gc - 16.4) / len(sequence), 1)
    except Exception:
        return None


def calculate_gc_content(dna_sequence):
    try:
        sequence = str(dna_sequence).upper()
        if not sequence:
            return None
        gc = sequence.count("G") + sequence.count("C")
        return round((gc / len(sequence)) * 100, 1)
    except Exception:
        return None


def revcomp(dna_sequence):
    return str(Seq(str(dna_sequence)).reverse_complement()).upper()


def has_bad_run(sequence, max_run=5):
    """max_run veya daha uzun homopolymer (AAAAA / GGGGG) var mi?"""
    run = 1
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i - 1]:
            run += 1
            if run >= max_run:
                return True
        else:
            run = 1
    return False


def dimer_risk(seq1, seq2=None, min_overlap=5):
    """Antiparalel min_overlap baz tamamlayici eslesme (primer-dimer) kontrolu."""
    a = str(seq1).upper()
    b = revcomp(seq2) if seq2 else revcomp(seq1)
    for i in range(len(a) - min_overlap + 1):
        kmer = a[i:i + min_overlap]
        if kmer in b:
            return True, kmer
    return False, None


def dimer_penalty(seq1, seq2=None):
    """Dimer riski icin yumusak ceza puani (5-mer siddetli, 4-mer hafif)."""
    a = str(seq1).upper()
    b = revcomp(seq2) if seq2 else revcomp(seq1)
    penalty = 0
    for i in range(len(a) - 4 + 1):
        kmer = a[i:i + 4]
        if kmer in b:
            penalty += 8
    if dimer_risk(seq1, seq2, min_overlap=5)[0]:
        penalty += 20
    return penalty


def validate_primer_primer(forward, reverse):
    """(bool, mesaj) dondurur: F/R dimer ve self-dimer kontrolu (5+ baz)."""
    problems = []
    risk, kmer = dimer_risk(forward, reverse)
    if risk:
        problems.append(f"F/R dimer ({kmer})")
    risk, kmer = dimer_risk(forward)
    if risk:
        problems.append(f"F self-dimer ({kmer})")
    risk, kmer = dimer_risk(reverse)
    if risk:
        problems.append(f"R self-dimer ({kmer})")
    if problems:
        return False, "Dimer riski: " + ", ".join(problems)
    return True, "Primer-dimer / capraz reaksiyon riski yok (5+ baz)"


def _mean(values):
    """None (veri yok) degerlerini atlar; hic veri yoksa 1.0 (ceza yok)."""
    vals = [v for v in (values or []) if v is not None]
    return sum(vals) / len(vals) if vals else 1.0


def design_primers(dna_sequence, start_pos, end_pos, target_amplicon_size=150, top_n=8, strict=True,
                   conservation=None, min_cons=0.85):
    """
    Hedef bolge icinde gercek amplikon mantigina uygun primer ciftleri tasarlar:
      forward = sense dizinin basi
      reverse = amplikon sonunun TERS KOMPLEMENTI
    Puanlamaya gore siralanmis aday listesi dondurur.
    """
    sequence = str(dna_sequence).upper()
    region = sequence[start_pos:end_pos]
    if len(region) < 60:
        return []

    primer_lengths = (18, 19, 20, 21, 22, 23, 24, 25)

    tm_cache, gc_cache = {}, {}

    def tm_of(seq):
        if seq not in tm_cache:
            tm_cache[seq] = calculate_tm(seq)
        return tm_cache[seq]

    def gc_of(seq):
        if seq not in gc_cache:
            gc_cache[seq] = calculate_gc_content(seq)
        return gc_cache[seq]

    candidates = []
    amp_min = max(60, target_amplicon_size - 35)
    amp_max = min(len(region), target_amplicon_size + 35)

    # Yumusak kabul araliklari (puanlama ile en iyiler secilir)
    if strict:
        tm_lo, tm_hi = PRIMER_TM_MIN, PRIMER_TM_MAX      # 58-65 C
        gc_lo, gc_hi = PRIMER_GC_MIN, PRIMER_GC_MAX      # %40-60
    else:
        tm_lo, tm_hi = 55.0, 68.0
        gc_lo, gc_hi = 30.0, 70.0

    for amp_len in range(amp_min, amp_max + 1, 3):
        max_offset = len(region) - amp_len
        if max_offset < 0:
            continue
        step = max(1, max_offset // 150 + 1)
        for offset in range(0, max_offset + 1, step):
            window = region[offset:offset + amp_len]
            if "N" in window:
                continue

            for f_len in primer_lengths:
                forward = window[:f_len]
                if "N" in forward or has_bad_run(forward, 6):
                    continue
                f_tm, f_gc = tm_of(forward), gc_of(forward)
                if f_tm is None or f_gc is None:
                    continue
                if not (tm_lo <= f_tm <= tm_hi) or not (gc_lo <= f_gc <= gc_hi):
                    continue

                for r_len in primer_lengths:
                    reverse_site = window[-r_len:]
                    if "N" in reverse_site:
                        continue
                    reverse = revcomp(reverse_site)
                    if has_bad_run(reverse, 6):
                        continue
                    r_tm, r_gc = tm_of(reverse), gc_of(reverse)
                    if r_tm is None or r_gc is None:
                        continue
                    if not (tm_lo <= r_tm <= tm_hi) or not (gc_lo <= r_gc <= gc_hi):
                        continue

                    # --- Korunmusluk filtresi ---
                    cons_f = cons_r = cons_amp = 1.0
                    if conservation is not None and len(conservation) >= len(region):
                        cons_f = _mean(conservation[offset:offset + f_len])
                        cons_r = _mean(conservation[offset + amp_len - r_len:offset + amp_len])
                        cons_amp = _mean(conservation[offset:offset + amp_len])
                        if cons_f < min_cons or cons_r < min_cons:
                            continue

                    if strict:
                        if forward[-1] not in "GC" or reverse[-1] not in "GC":
                            continue
                        if has_bad_run(forward) or has_bad_run(reverse):
                            continue
                        if abs(f_tm - r_tm) > PRIMER_TM_MAX_DIFF:
                            continue
                        if not validate_primer_primer(forward, reverse)[0]:
                            continue

                    # --- Puanlama (yuksek = iyi) ---
                    score = 100.0
                    score -= abs(f_tm - 61.0) * 3.0
                    score -= abs(r_tm - 61.0) * 3.0
                    score -= abs(f_gc - 50.0) * 1.0
                    score -= abs(r_gc - 50.0) * 1.0
                    score -= abs(f_tm - r_tm) * 4.0
                    score -= abs(amp_len - target_amplicon_size) * 1.2
                    score -= abs((f_gc + r_gc) / 2.0 - 50.0) * 0.8
                    if forward[-1] in "GC":
                        score += 6
                    if reverse[-1] in "GC":
                        score += 6
                    if has_bad_run(forward):
                        score -= 8
                    if has_bad_run(reverse):
                        score -= 8
                    score -= dimer_penalty(forward, reverse)
                    score -= dimer_penalty(forward) * 0.5
                    score -= dimer_penalty(reverse) * 0.5
                    if conservation is not None:
                        score += (cons_f + cons_r) * 12.0 + cons_amp * 8.0

                    candidates.append({
                        "forward": forward,
                        "reverse": reverse,
                        "forward_tm": f_tm,
                        "reverse_tm": r_tm,
                        "forward_gc": f_gc,
                        "reverse_gc": r_gc,
                        "amplicon_size": amp_len,
                        "amplicon_seq": window,
                        "forward_start": start_pos + offset + 1,
                        "reverse_end": start_pos + offset + amp_len,
                        "forward_len": f_len,
                        "reverse_len": r_len,
                        "score": round(score, 1),
                    })

    if not candidates:
        return []

    candidates.sort(key=lambda c: c["score"], reverse=True)

    unique, seen = [], set()
    for c in candidates:
        key = (c["forward_start"], c["amplicon_size"])
        if key in seen:
            continue
        seen.add(key)
        unique.append(c)
        if len(unique) >= top_n:
            break
    return unique


def design_primers_best(dna_sequence, start_pos, end_pos, target_amplicon_size=150, top_n=8,
                        conservation=None, min_cons=0.85):
    """
    Once katı kriterlerle (Tm 58-65, GC %40-60, dimer yok, 3' G/C) arar.
    Bulamazsa gevşek kriterlere düşer ve puanlamaya güvenir.
    conservation: hedef bölge koordinatlarıyla aynı uzunlukta korunmuşluk dizisi.
    Dönen değer: (adaylar, strict_mi)
    """
    strict_hits = design_primers(
        dna_sequence, start_pos, end_pos, target_amplicon_size, top_n, strict=True,
        conservation=conservation, min_cons=min_cons,
    )
    if strict_hits:
        return strict_hits, True
    return design_primers(
        dna_sequence, start_pos, end_pos, target_amplicon_size, top_n, strict=False,
        conservation=conservation, min_cons=min_cons,
    ), False


def design_probe(amplicon_seq, forward_len, reverse_len, target_tm=69.0,
                 conservation=None, min_cons=0.85,
                 tm_min=None, tm_max=None, gc_min=None, gc_max=None,
                 avoid_sequences=(), avoid_overlap=6):
    """
    TaqMan probu: primerlerin arasinda, 5' ucunda G yok,
    Tm ~69 C, 18-30 nt, homopolymer siniri 4.
    conservation verilirse amplikon koordinatlariyla eslesen dizi beklenir.
    tm_min/tm_max/gc_min/gc_max: AT-zengin bolgeler icin gevsetilebilir.
    """
    tm_lo = PROBE_TM_MIN if tm_min is None else float(tm_min)
    tm_hi = PROBE_TM_MAX if tm_max is None else float(tm_max)
    gc_lo = PROBE_GC_MIN if gc_min is None else float(gc_min)
    gc_hi = PROBE_GC_MAX if gc_max is None else float(gc_max)
    amplicon = str(amplicon_seq).upper()
    inner_start = forward_len
    inner_end = len(amplicon) - reverse_len
    inner = amplicon[inner_start:inner_end]
    if len(inner) < 18:
        return None

    best_clean = None
    best_any = None
    for plen in range(18, min(31, len(inner)) + 1):
        for i in range(0, len(inner) - plen + 1):
            probe = inner[i:i + plen]
            if "N" in probe or probe[0] == "G" or has_bad_run(probe, 4):
                continue
            tm, gc = calculate_tm(probe), calculate_gc_content(probe)
            if tm is None or gc is None:
                continue
            if not (tm_lo <= tm <= tm_hi):
                continue
            if not (gc_lo <= gc <= gc_hi):
                continue

            cons = 1.0
            if conservation is not None:
                start = inner_start + i
                end = start + plen
                if len(conservation) >= end:
                    cons = _mean(conservation[start:end])
                    if cons < min_cons:
                        continue

            risky = bool(avoid_sequences) and any(
                dimer_risk(probe, other, min_overlap=avoid_overlap)[0]
                for other in avoid_sequences
            )

            score = -abs(tm - target_tm) * 5 - abs(gc - 55.0) * 0.8 - abs(plen - 25) * 0.5
            if conservation is not None:
                score += cons * 15.0
            entry = {
                "sequence": probe,
                "tm": tm,
                "gc": gc,
                "length": plen,
                "position": inner_start + i + 1,
                "conservation": round(cons * 100, 1) if conservation is not None else None,
                "dimer_risk": risky,
                "_score": score,
            }
            if best_any is None or score > best_any["_score"]:
                best_any = entry
            if not risky and (best_clean is None or score > best_clean["_score"]):
                best_clean = entry

    best = best_clean or best_any
    if best:
        best = dict(best)
        best.pop("_score", None)
    return best


def check_specificity(primer, sequences, min_identity_region=0.85):
    """
    Basit in silico ozgulluk: primerin diger dizilerde tam/neredeyse tam
    eslesip eslesmedigini kontrol eder (Primer-BLAST yerine hizli on kontrol).
    """
    primer = str(primer).upper()
    hits = []
    for name, seq in sequences:
        seq = str(seq).upper()
        if primer in seq or revcomp(primer) in seq:
            hits.append(name)
    return hits
