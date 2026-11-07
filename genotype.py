"""
genotype.py
Genotip kapsama analizi:
  - toleransli in-silico PCR ile her dizide amplikon bulma
  - referans amplikona global hizalama
  - primer/prob bazinda mismatch, 3' uc riski, Wilson GA
  - mismatch isi haritasi, entropi/korunmusluk izi, sequence logo
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from Bio import Align

from browser_track import navy_template
from conservation import fetch_related_sequences
from pcr_engine import revcomp
from specificity import BindingSite, in_silico_pcr, find_binding_sites
from stats_engine import wilson_ci


@dataclass
class OligoCoverage:
    name: str
    sequence: str
    n_total: int
    n_evaluable: int
    perfect: int = 0
    one: int = 0
    two: int = 0
    three_plus: int = 0
    mean_mm: float | None = None
    three_prime_risk: int = 0
    perfect_ci: tuple[float, float] = (0.0, 0.0)
    evaluable_ci: tuple[float, float] = (0.0, 0.0)
    perfect_pct: float = 0.0
    evaluable_pct: float = 0.0
    three_prime_risk_pct: float = 0.0

    def as_dict(self):
        return {
            "name": self.name,
            "sequence": self.sequence,
            "n_total": self.n_total,
            "n_evaluable": self.n_evaluable,
            "perfect": self.perfect,
            "one": self.one,
            "two": self.two,
            "three_plus": self.three_plus,
            "mean_mm": round(self.mean_mm, 2) if self.mean_mm is not None else None,
            "three_prime_risk": self.three_prime_risk,
            "perfect_pct": round(100 * self.perfect_pct, 1),
            "perfect_ci": [round(100 * x, 1) for x in self.perfect_ci],
            "evaluable_pct": round(100 * self.evaluable_pct, 1),
            "evaluable_ci": [round(100 * x, 1) for x in self.evaluable_ci],
            "three_prime_risk_pct": round(100 * self.three_prime_risk_pct, 1),
        }


def _align_to_reference(ref_amp: str, subject_amp: str):
    """Global hizalama; referans pozisyonlari icin subject bazlarini dondurur."""
    ref_amp = str(ref_amp).upper()
    subject_amp = str(subject_amp).upper()
    if not ref_amp or not subject_amp:
        return [None] * len(ref_amp), None, None
    if ref_amp == subject_amp:
        return list(subject_amp), 1.0, 0

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    try:
        aln = aligner.align(ref_amp, subject_amp)[0]
        ref_aln = str(aln[0])
        subj_aln = str(aln[1])
    except Exception:
        return [None] * len(ref_amp), None, None

    bases: list[str | None] = [None] * len(ref_amp)
    ri = 0
    matches = compared = gaps = 0
    for a, b in zip(ref_aln, subj_aln):
        if a == "-":
            continue
        if ri >= len(bases):
            break
        bases[ri] = b
        if b == "-":
            gaps += 1
        else:
            compared += 1
            if a == b:
                matches += 1
        ri += 1
    identity = (matches / compared) if compared else None
    return bases, identity, gaps


def _mm_from_aligned(ref_amp: str, aligned: list[str | None], start: int, end: int,
                     three_prime_side: str = "right"):
    """Hizalanmis subject bazlarindan oligo mismatch sayilari."""
    if aligned is None or start is None or end is None:
        return None, None, []
    start = max(0, int(start))
    end = min(len(ref_amp), int(end))
    if end <= start:
        return None, None, []
    mm_idx = []
    for i in range(start, end):
        base = aligned[i] if i < len(aligned) else None
        ref_base = ref_amp[i]
        if base is None or base == "-" or base != ref_base:
            mm_idx.append(i - start)
    if three_prime_side == "left":
        three = sum(1 for i in mm_idx if i < 5)
    elif three_prime_side == "right":
        three = sum(1 for i in mm_idx if i >= (end - start) - 5)
    else:
        three = 0
    return len(mm_idx), three, mm_idx


def _best_site(sites):
    if not sites:
        return None
    return sorted(sites, key=lambda s: (s.n_mismatches, s.three_prime_mm, s.position))[0]


def _independent_sites(seq: str, forward: str, reverse: str, probe: str | None,
                       max_mismatches: int = 5):
    f = _best_site(find_binding_sites(forward, seq, max_mismatches, "right"))
    r = _best_site(find_binding_sites(revcomp(reverse), seq, max_mismatches, "left"))
    p = _best_site(find_binding_sites(probe, seq, max_mismatches, "right")) if probe else None
    return f, r, p


def _record_for_subject(subject: dict, ref_amp: str, forward: str, reverse: str,
                        probe: str | None, probe_offset: int | None,
                        max_mismatches: int, expected_size: int):
    seq = str(subject.get("seq", "")).upper()
    sid = str(subject.get("id", "?"))
    desc = str(subject.get("description", ""))[:160]
    rec = {
        "id": sid,
        "description": desc,
        "found": False,
        "partial": False,
        "amplicon_start": None,
        "amplicon_end": None,
        "amplicon_size": None,
        "identity": None,
        "gaps": None,
        "aligned": None,
        "F": None,
        "R": None,
        "P": None,
        "total_mm": None,
        "three_prime_mm": None,
        "error": None,
    }
    if not seq:
        rec["error"] = "Bos dizi"
        return rec

    hits = in_silico_pcr(
        forward, reverse, seq, probe,
        max_mismatches=max_mismatches,
        min_amp=max(30, int(expected_size * 0.5)),
        max_amp=max(100, int(expected_size * 1.6)),
    )
    if hits:
        # probu olan, az mismatchli ve beklenen boya yakin hit'i sec
        def key(h):
            return (0 if h.probe is not None else 1, h.total_mismatches,
                    abs(h.size - expected_size), h.three_prime_mismatches)
        hit = sorted(hits, key=key)[0]
        aligned, identity, gaps = _align_to_reference(ref_amp, hit.amplicon)
        rec.update({
            "found": True,
            "partial": False,
            "amplicon_start": hit.amplicon_start,
            "amplicon_end": hit.amplicon_end,
            "amplicon_size": hit.size,
            "identity": round(identity, 4) if identity is not None else None,
            "gaps": gaps,
            "aligned": aligned,
        })
        f_mm, f_three, _ = _mm_from_aligned(ref_amp, aligned, 0, len(forward), "right")
        r_start = len(ref_amp) - len(reverse)
        r_mm, r_three, _ = _mm_from_aligned(ref_amp, aligned, r_start, len(ref_amp), "left")
        p_mm = p_three = None
        if probe and probe_offset is not None:
            p_mm, p_three, _ = _mm_from_aligned(
                ref_amp, aligned, probe_offset, probe_offset + len(probe), "none"
            )
        rec["F"] = {"n_mismatches": f_mm, "three_prime_mm": f_three}
        rec["R"] = {"n_mismatches": r_mm, "three_prime_mm": r_three}
        rec["P"] = {"n_mismatches": p_mm, "three_prime_mm": p_three}
        rec["total_mm"] = sum(x for x in (f_mm, r_mm, p_mm) if x is not None)
        rec["three_prime_mm"] = (f_three or 0) + (r_three or 0)
        return rec

    # Amplikon bulunamadi: oligolari bagimsiz ara (kismi diziler icin)
    f, r, p = _independent_sites(seq, forward, reverse, probe, max_mismatches=max_mismatches + 1)
    rec["partial"] = True
    if f:
        rec["F"] = {"n_mismatches": f.n_mismatches, "three_prime_mm": f.three_prime_mm}
    if r:
        rec["R"] = {"n_mismatches": r.n_mismatches, "three_prime_mm": r.three_prime_mm}
    if p:
        rec["P"] = {"n_mismatches": p.n_mismatches, "three_prime_mm": 0}
    if any(rec[k] is not None for k in ("F", "R", "P")):
        rec["total_mm"] = sum(
            rec[k]["n_mismatches"] for k in ("F", "R", "P") if rec[k] is not None
        )
        rec["three_prime_mm"] = sum(
            rec[k]["three_prime_mm"] for k in ("F", "R") if rec[k] is not None
        )
    else:
        rec["error"] = "Amplikon ve oligo baglanma bolgesi bulunamadi"
    return rec


def _coverage_for(name: str, sequence: str, records: list[dict], n_total: int) -> OligoCoverage:
    values = []
    three_risk = 0
    for rec in records:
        info = rec.get(name)
        if info is None or info.get("n_mismatches") is None:
            continue
        mm = int(info["n_mismatches"])
        values.append(mm)
        if name in ("F", "R") and int(info.get("three_prime_mm") or 0) > 0:
            three_risk += 1
    n_eval = len(values)
    cov = OligoCoverage(name=name, sequence=sequence, n_total=n_total, n_evaluable=n_eval)
    for mm in values:
        if mm == 0:
            cov.perfect += 1
        elif mm == 1:
            cov.one += 1
        elif mm == 2:
            cov.two += 1
        else:
            cov.three_plus += 1
    if values:
        cov.mean_mm = float(np.mean(values))
        cov.perfect_pct = cov.perfect / n_eval if n_eval else 0.0
        cov.perfect_ci = wilson_ci(cov.perfect, n_eval)
        cov.three_prime_risk_pct = three_risk / n_eval if n_eval else 0.0
    cov.three_prime_risk = three_risk
    cov.evaluable_pct = n_eval / n_total if n_total else 0.0
    cov.evaluable_ci = wilson_ci(n_eval, n_total)
    return cov


def analyze_genotypes(ref_seq: str, amplicon_start: int, amplicon_end: int,
                      forward: str, reverse: str, probe: str | None,
                      subjects: list[dict], probe_offset: int | None = None,
                      max_mismatches: int = 4, include_reference: bool = True):
    """
    Ana giris noktasi.

    subjects: [{'id','description','seq'}, ...]
    probe_offset: probun amplikon icindeki 0-based baslangici.
    """
    ref_seq = str(ref_seq).upper()
    amplicon_start = max(0, int(amplicon_start))
    amplicon_end = min(len(ref_seq), int(amplicon_end))
    ref_amp = ref_seq[amplicon_start:amplicon_end]
    if not ref_amp:
        return {"error": "Bos amplikon."}
    if probe_offset is None and probe:
        probe_offset = ref_amp.find(str(probe).upper())
        if probe_offset < 0:
            probe_offset = None

    records = []
    for subj in subjects:
        try:
            rec = _record_for_subject(
                subj, ref_amp, forward, reverse, probe, probe_offset,
                max_mismatches=max_mismatches, expected_size=len(ref_amp),
            )
        except Exception as exc:
            rec = {
                "id": str(subj.get("id", "?")), "description": str(subj.get("description", ""))[:160],
                "found": False, "partial": False, "aligned": None, "F": None, "R": None, "P": None,
                "total_mm": None, "three_prime_mm": None, "error": str(exc)[:120],
            }
        records.append(rec)

    n_total = len(records)
    cov_f = _coverage_for("F", forward, records, n_total)
    cov_r = _coverage_for("R", reverse, records, n_total)
    cov_p = _coverage_for("P", probe or "", records, n_total) if probe else None
    coverages = [c for c in (cov_f, cov_r, cov_p) if c is not None]

    found = [r for r in records if r["found"]]
    amplicon_pct = len(found) / n_total if n_total else 0.0
    amplicon_ci = wilson_ci(len(found), n_total)

    # Mismatch matrisi (referans satiri dahil)
    rows = []
    labels = []
    if include_reference:
        rows.append([0] * len(ref_amp))
        labels.append("REF")
    for rec in records:
        aligned = rec.get("aligned")
        if aligned is None:
            rows.append([None] * len(ref_amp))
        else:
            vals = []
            for i, base in enumerate(aligned):
                if base is None:
                    vals.append(None)
                elif base == "-":
                    vals.append(2)
                elif base == ref_amp[i]:
                    vals.append(0)
                else:
                    vals.append(1)
            rows.append(vals)
        labels.append(rec["id"])

    # Pozisyon bazli entropi / korunmusluk (referans haric)
    depth_list, match_list, gap_list = [], [], []
    base_counts = []
    for i in range(len(ref_amp)):
        counts = {"A": 0, "C": 0, "G": 0, "T": 0}
        gaps = 0
        for rec in records:
            aligned = rec.get("aligned")
            if aligned is None or i >= len(aligned):
                continue
            b = aligned[i]
            if b is None:
                continue
            if b == "-":
                gaps += 1
            elif b in counts:
                counts[b] += 1
        n = sum(counts.values())
        depth_list.append(n)
        gap_list.append(gaps)
        match_list.append(counts.get(ref_amp[i], 0))
        base_counts.append(counts)

    conservation = []
    entropy = []
    for i, n in enumerate(depth_list):
        if n <= 0:
            conservation.append(None)
            entropy.append(None)
            continue
        conservation.append(match_list[i] / n)
        probs = [base_counts[i][b] / n for b in "ACGT"]
        h = -sum(p * math.log2(p) for p in probs if p > 0)
        entropy.append(h)

    logo_rows = []
    logo_positions = []
    for i, counts in enumerate(base_counts):
        n = sum(counts.values())
        logo_positions.append(amplicon_start + i + 1)
        if n <= 0:
            logo_rows.append([0.25, 0.25, 0.25, 0.25])
        else:
            logo_rows.append([counts[b] / n for b in "ACGT"])
    logo_df = pd.DataFrame(logo_rows, columns=list("ACGT"), index=logo_positions)

    return {
        "ref_amp": ref_amp,
        "amplicon_start": amplicon_start,
        "amplicon_end": amplicon_end,
        "probe_offset": probe_offset,
        "forward_len": len(forward),
        "reverse_len": len(reverse),
        "probe_len": len(probe) if probe else 0,
        "records": records,
        "labels": labels,
        "matrix": rows,
        "matrix_labels": labels,
        "depth": depth_list,
        "gaps": gap_list,
        "conservation": conservation,
        "entropy": entropy,
        "base_counts": base_counts,
        "logo_df": logo_df,
        "coverage": {c.name: c.as_dict() for c in coverages},
        "coverages": coverages,
        "n_total": n_total,
        "n_found": len(found),
        "amplicon_found_pct": round(100 * amplicon_pct, 1),
        "amplicon_found_ci": [round(100 * x, 1) for x in amplicon_ci],
        "identity_mean": round(float(np.mean([r["identity"] for r in found if r["identity"] is not None])), 4)
        if any(r["identity"] is not None for r in found) else None,
        "conservation_mean": round(float(np.mean([c for c in conservation if c is not None])), 4)
        if any(c is not None for c in conservation) else None,
        "entropy_mean": round(float(np.mean([e for e in entropy if e is not None])), 3)
        if any(e is not None for e in entropy) else None,
    }


def fetch_subjects(organism: str, exclude_accessions, reference_seq: str, n: int = 12,
                   keywords=(), max_candidates=None, min_subject_len=None,
                   country: str | None = None, min_similarity=None):
    """NCBI'den ayni organizmaya ait dizileri ceker; genotype.analyze_genotypes icin hazirlar.

    country verilirse yalnizca o ulkeye ait kayitlar aranir (or. Turkey).
    min_similarity=0.0 ulke taramasinda farkli lineage'lari da getirir; hedef bolge
    yoksa analyze_genotypes bunu 'degerlendirilemez' olarak isaretler.
    """
    kwargs = {}
    if min_similarity is not None:
        kwargs["min_similarity"] = float(min_similarity)
    rows = fetch_related_sequences(
        organism,
        exclude_accessions=exclude_accessions,
        max_n=n,
        reference_seq=reference_seq,
        keywords=tuple(keywords),
        max_candidates=max_candidates,
        min_subject_len=min_subject_len,
        country=country,
        **kwargs,
    )
    return [{"id": acc, "description": desc, "seq": seq} for acc, desc, seq in rows]


# ---------------------------------------------------------------------------
# Gorseller
# ---------------------------------------------------------------------------

def heatmap_figure(result: dict, max_rows: int = 80):
    """Mismatch isi haritasi + entropi/korunmusluk izi."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    matrix = result["matrix"]
    labels = result["matrix_labels"]
    if not matrix:
        return None
    if len(matrix) > max_rows:
        matrix = matrix[:max_rows]
        labels = labels[:max_rows]
    x = list(range(result["amplicon_start"] + 1, result["amplicon_end"] + 1))
    n_rows = len(matrix)
    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
        row_heights=[0.28, 0.72], specs=[[{"secondary_y": True}], [{}]],
        subplot_titles=("Korunmuşluk ve entropi", "Dizi bazında mismatch haritası"),
    )
    cons = result["conservation"]
    ent = result["entropy"]
    fig.add_trace(go.Scatter(
        x=x, y=[None if c is None else 100 * c for c in cons],
        mode="lines", name="Korunmuşluk (%)", line=dict(color="#1b5e20", width=2),
        hovertemplate="%{x} bp<br>Korunmuşluk: %{y:.1f}%<extra></extra>",
    ), row=1, col=1, secondary_y=False)
    fig.add_trace(go.Scatter(
        x=x, y=ent, mode="lines", name="Entropi (bit)", line=dict(color="#c62828", width=1.5),
        fill="tozeroy", opacity=0.25,
        hovertemplate="%{x} bp<br>Entropi: %{y:.2f} bit<extra></extra>",
    ), row=1, col=1, secondary_y=True)

    fig.add_trace(go.Heatmap(
        z=matrix, x=x, y=labels, zmin=0, zmax=2,
        colorscale=[[0.0, "#c8e6c9"], [0.49, "#c8e6c9"], [0.5, "#e53935"],
                    [0.99, "#e53935"], [1.0, "#fb8c00"]],
        showscale=False, hoverongaps=False,
        hovertemplate="Pozisyon %{x} bp<br>Dizi: %{y}<br>Durum: %{z}<extra></extra>",
    ), row=2, col=1)

    # primer/prob bolgeleri
    regions = [
        ("F", result["amplicon_start"] + 1, result["amplicon_start"] + result["forward_len"], "#1565c0"),
        ("P", result["amplicon_start"] + (result.get("probe_offset") or 0) + 1,
         result["amplicon_start"] + (result.get("probe_offset") or 0) + result["probe_len"], "#f9a825")
        if result.get("probe_len") else None,
        ("R", result["amplicon_end"] - result["reverse_len"] + 1, result["amplicon_end"], "#6a1b9a"),
    ]
    for reg in regions:
        if reg is None:
            continue
        name, x0, x1, color = reg
        fig.add_vrect(x0=x0 - 0.5, x1=x1 + 0.5, fillcolor=color, opacity=0.10,
                      line_width=1, line_color=color, row=2, col=1)
        fig.add_annotation(x=(x0 + x1) / 2, y=1.02, yref="y2 domain", text=name,
                           showarrow=False, font=dict(color=color, size=11))

    fig.update_yaxes(title_text="Korunmuşluk (%)", range=[0, 105], row=1, col=1, secondary_y=False)
    fig.update_yaxes(title_text="Entropi", range=[0, 2.1], row=1, col=1, secondary_y=True)
    fig.update_yaxes(title_text="", row=2, col=1, tickfont=dict(size=9))
    fig.update_xaxes(title_text="Genomik pozisyon (bp)", row=2, col=1)
    height = max(430, 120 + 18 * n_rows)
    fig.update_layout(
        height=height, margin=dict(l=20, r=20, t=50, b=40),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0.01),
        template=navy_template(),
    )
    return fig


def coverage_figure(result: dict):
    """Primer/prob basina mismatch dagilimi (stacked bar)."""
    import plotly.graph_objects as go

    coverages = result.get("coverages") or []
    if not coverages:
        return None
    names = [c.name for c in coverages]
    cats = [("perfect", "0 mismatch", "#2e7d32"),
            ("one", "1 mismatch", "#f9a825"),
            ("two", "2 mismatch", "#ef6c00"),
            ("three_plus", "3+ mismatch", "#c62828")]
    fig = go.Figure()
    for key, label, color in cats:
        vals = [getattr(c, key) for c in coverages]
        fig.add_trace(go.Bar(
            y=names, x=vals, name=label, orientation="h", marker_color=color,
            hovertemplate="%{y}: %{x} dizi<extra>" + label + "</extra>",
        ))
    for c in coverages:
        fig.add_annotation(
            x=c.n_evaluable, y=c.name, text=f"değerlendirilebilir {c.n_evaluable}/{c.n_total}",
            showarrow=False, xanchor="left", xshift=5, font=dict(size=9, color="#E6EDF7"),
        )
    fig.update_layout(
        barmode="stack", height=250, margin=dict(l=10, r=10, t=30, b=30),
        xaxis_title="Dizi sayısı", yaxis_title="", legend=dict(orientation="h", y=1.12),
        template=navy_template(),
    )
    return fig


def logo_figure(result: dict, max_bases: int = 180):
    """Sequence logo (logomaker/matplotlib)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import logomaker
    except Exception:
        return None
    df = result.get("logo_df")
    if df is None or df.empty:
        return None
    if len(df) > max_bases:
        # en yuksek entropili pencereyi sec
        ent = result.get("entropy") or []
        if len(ent) == len(df):
            score = [(-1 if e is None else e) for e in ent]
            start = 0
            best = sum(score[:max_bases])
            cur = best
            for i in range(1, len(score) - max_bases + 1):
                cur += score[i + max_bases - 1] - score[i - 1]
                if cur > best:
                    best, start = cur, i
            df = df.iloc[start:start + max_bases]
    try:
        info = logomaker.transform_matrix(df, from_type="probability", to_type="information")
        fig, ax = plt.subplots(figsize=(max(6, len(df) / 12), 1.9))
        fig.patch.set_facecolor("#0B1526")
        ax.set_facecolor("#0B1526")
        logomaker.Logo(info, ax=ax, color_scheme="classic", vpad=0.02)
        ax.set_xlabel("Amplikon pozisyonu (bp)", fontsize=8, color="#E6EDF7")
        ax.set_ylabel("Bilgi (bit)", fontsize=8, color="#E6EDF7")
        ax.tick_params(labelsize=7, colors="#E6EDF7")
        fig.tight_layout()
        return fig
    except Exception:
        return None
