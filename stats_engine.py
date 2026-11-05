"""
stats_engine.py
Biyoistatistik ve validasyon motoru:
  - Wilson guven araligi
  - probit MLE + bootstrap ile LOD95
  - Monte Carlo qPCR simulasyonu (Poisson dagilim, verim degiskenligi, Ct gurultusu)
  - LOQ (CV tabanli), standart egri, PPV/NPV, orneklem buyuklugu
"""

from __future__ import annotations

import math
import random

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Basari orani icin Wilson skor guven araligi (0-1)."""
    if n <= 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def pct(x: float) -> float:
    return round(100.0 * float(x), 1)


def probit_nll(params, x, y):
    a, b = params
    eta = a + b * x
    p = np.clip(norm.cdf(eta), 1e-9, 1 - 1e-9)
    return -float(np.sum(y * np.log(p) + (1 - y) * np.log(1 - p)))


def probit_fit(log_copies, detected):
    """
    log10(kopya) ~ tespit olasiligi probit modeli.
    Donus: {'a','b','lod95_log10','lod95','n','ok','message'}
    """
    x = np.asarray(log_copies, dtype=float)
    y = np.asarray(detected, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 6 or len(set(y.tolist())) < 2:
        return {"ok": False, "message": "Probit icin en az 6 gozlem ve hem pozitif hem negatif sonuc gerekir.",
                "n": int(n), "a": None, "b": None, "lod95_log10": None, "lod95": None}

    best = None
    for start in ((0.0, 1.0), (-1.0, 1.5), (1.0, 0.8)):
        res = minimize(probit_nll, start, args=(x, y), method="Nelder-Mead",
                       options={"maxiter": 5000, "xatol": 1e-8, "fatol": 1e-8})
        if best is None or res.fun < best.fun:
            best = res
    a, b = float(best.x[0]), float(best.x[1])
    if abs(b) < 1e-9:
        return {"ok": False, "message": "Probit egimi sifir; LOD hesaplanamadi.",
                "n": int(n), "a": a, "b": b, "lod95_log10": None, "lod95": None}
    lod_log = (norm.ppf(0.95) - a) / b
    return {
        "ok": True,
        "a": round(a, 4),
        "b": round(b, 4),
        "lod95_log10": round(float(lod_log), 3),
        "lod95": round(float(10 ** lod_log), 2),
        "n": int(n),
        "nll": round(float(best.fun), 3),
        "message": "",
    }


def probit_bootstrap(log_copies, detected, n_boot=300, seed=42):
    """LOD95 icin bootstrap %95 GA."""
    x = np.asarray(log_copies, dtype=float)
    y = np.asarray(detected, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 6:
        return None
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(max(50, int(n_boot))):
        idx = rng.integers(0, len(x), len(x))
        if len(set(y[idx].tolist())) < 2:
            continue
        fit = probit_fit(x[idx], y[idx])
        if fit.get("ok") and fit.get("lod95") is not None:
            vals.append(fit["lod95"])
    if not vals:
        return None
    vals.sort()
    lo = vals[int(0.025 * (len(vals) - 1))]
    hi = vals[int(0.975 * (len(vals) - 1))]
    return {"low": round(lo, 2), "high": round(hi, 2), "n_boot": len(vals)}


def simulate_qpcr(copy_grid, replicates=24, efficiency=0.90, efficiency_sd=0.05,
                  max_cycles=40, threshold_copies=1e10, ct_sd=0.30,
                  inhibition_prob=0.0, seed=42):
    """
    Basit ama aciklanabilir qPCR Monte Carlo modeli.

    N0 ~ Poisson(kopya)
    Ct = ln(threshold_copies / N0) / ln(1 + verim) + Normal(0, ct_sd)
    Tespit: N0 > 0, inhibisyon yok ve Ct <= max_cycles
    """
    rng = np.random.default_rng(seed)
    rows = []
    for copies in copy_grid:
        copies = float(copies)
        for rep in range(int(replicates)):
            n0 = int(rng.poisson(copies)) if copies > 0 else 0
            eff = float(np.clip(rng.normal(efficiency, efficiency_sd), 0.50, 0.99))
            inhibited = bool(rng.random() < inhibition_prob)
            if n0 <= 0 or inhibited:
                ct = None
                detected = False
            else:
                ct = math.log(threshold_copies / n0) / math.log(1.0 + eff)
                ct += float(rng.normal(0.0, ct_sd))
                if ct < 0.5:
                    ct = 0.5
                detected = bool(ct <= max_cycles)
                if not detected:
                    ct = None
            rows.append({
                "copies": copies,
                "replicate": rep + 1,
                "n0": n0,
                "efficiency": round(eff, 3),
                "ct": round(float(ct), 2) if ct is not None else None,
                "detected": detected,
                "inhibited": inhibited,
            })
    return rows


def summarize_simulation(rows, copy_grid=None):
    """Kopya basina tespit orani, Ct ortalama/SD/CV."""
    if copy_grid is None:
        copy_grid = sorted({float(r["copies"]) for r in rows})
    summary = []
    for copies in copy_grid:
        subset = [r for r in rows if float(r["copies"]) == float(copies)]
        n = len(subset)
        detected = [r for r in subset if r["detected"] and r["ct"] is not None]
        cts = [r["ct"] for r in detected]
        mean_ct = float(np.mean(cts)) if cts else None
        sd_ct = float(np.std(cts, ddof=1)) if len(cts) > 1 else 0.0
        cv = (sd_ct / mean_ct) if (mean_ct and len(cts) > 1) else None
        lo, hi = wilson_ci(len(detected), n)
        summary.append({
            "copies": float(copies),
            "n": n,
            "detected": len(detected),
            "detection_rate": round(len(detected) / n, 4) if n else 0.0,
            "detection_ci_low": round(lo, 3),
            "detection_ci_high": round(hi, 3),
            "mean_ct": round(mean_ct, 2) if mean_ct is not None else None,
            "sd_ct": round(sd_ct, 2) if cts else None,
            "cv_ct": round(cv, 4) if cv is not None else None,
        })
    return summary


def lod_loq_analysis(copy_grid, replicates=24, efficiency=0.90, efficiency_sd=0.05,
                     max_cycles=40, threshold_copies=1e10, ct_sd=0.30,
                     inhibition_prob=0.0, seed=42, n_boot=300,
                     detection_target=0.95, loq_cv_target=0.05):
    """
    Monte Carlo + probit + LOQ.
    LOQ: tespit orani >= detection_target ve Ct CV <= loq_cv_target olan en dusuk kopya.
    """
    rows = simulate_qpcr(
        copy_grid, replicates=replicates, efficiency=efficiency,
        efficiency_sd=efficiency_sd, max_cycles=max_cycles,
        threshold_copies=threshold_copies, ct_sd=ct_sd,
        inhibition_prob=inhibition_prob, seed=seed,
    )
    summary = summarize_simulation(rows, copy_grid)
    x = [math.log10(float(r["copies"])) if float(r["copies"]) > 0 else -1.0 for r in rows]
    y = [1 if r["detected"] else 0 for r in rows]
    fit = probit_fit(x, y)
    ci = probit_bootstrap(x, y, n_boot=n_boot, seed=seed + 1) if fit.get("ok") else None

    loq = None
    for row in summary:
        if row["detection_rate"] >= detection_target and row["cv_ct"] is not None and row["cv_ct"] <= loq_cv_target:
            loq = row["copies"]
            break

    # Standart egri (tespit edilenler uzerinden)
    sc = standard_curve(rows)
    return {
        "rows": rows,
        "summary": summary,
        "probit": fit,
        "lod95": fit.get("lod95") if fit.get("ok") else None,
        "lod95_ci": ci,
        "loq": loq,
        "standard_curve": sc,
        "params": {
            "replicates": int(replicates),
            "efficiency": float(efficiency),
            "efficiency_sd": float(efficiency_sd),
            "max_cycles": int(max_cycles),
            "threshold_copies": float(threshold_copies),
            "ct_sd": float(ct_sd),
            "inhibition_prob": float(inhibition_prob),
            "detection_target": float(detection_target),
            "loq_cv_target": float(loq_cv_target),
        },
    }


def standard_curve(rows):
    """Ct = a + b * log10(kopya); PCR verimi = 10^(-1/b) - 1."""
    pts = [(math.log10(float(r["copies"])), r["ct"]) for r in rows
           if r["ct"] is not None and float(r["copies"]) > 0]
    if len(pts) < 3:
        return {"ok": False, "message": "Standart egri icin en az 3 nokta gerekir."}
    x = np.array([p[0] for p in pts], dtype=float)
    y = np.array([p[1] for p in pts], dtype=float)
    if len(set(x.tolist())) < 2:
        return {"ok": False, "message": "Standart egri icin en az 2 farkli kopya duzeyi gerekir."}
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    eff = 10 ** (-1.0 / slope) - 1.0 if slope != 0 else None
    return {
        "ok": True,
        "slope": round(float(slope), 3),
        "intercept": round(float(intercept), 3),
        "r2": round(float(r2), 4),
        "efficiency": round(float(eff), 4) if eff is not None else None,
        "n": len(pts),
    }


def ppv_npv(sensitivity, specificity, prevalence):
    """Bayes teoremi ile PPV/NPV."""
    sens = min(max(float(sensitivity), 0.0), 1.0)
    spec = min(max(float(specificity), 0.0), 1.0)
    prev = min(max(float(prevalence), 0.0), 1.0)
    ppv = (sens * prev) / (sens * prev + (1 - spec) * (1 - prev)) if (sens * prev + (1 - spec) * (1 - prev)) > 0 else 0.0
    npv = (spec * (1 - prev)) / (spec * (1 - prev) + (1 - sens) * prev) if (spec * (1 - prev) + (1 - sens) * prev) > 0 else 0.0
    return {
        "ppv": round(ppv, 4),
        "npv": round(npv, 4),
        "sensitivity": round(sens, 4),
        "specificity": round(spec, 4),
        "prevalence": round(prev, 4),
    }


def sample_size_sensitivity(expected=0.95, precision=0.05, z=1.96):
    """Duyarlilik icin yaklasik orneklem buyuklugu (Wald)."""
    expected = min(max(float(expected), 0.01), 0.99)
    precision = max(float(precision), 0.001)
    n = (z ** 2) * expected * (1 - expected) / (precision ** 2)
    return int(math.ceil(n))


def sample_size_specificity(expected=0.98, precision=0.02, z=1.96):
    return sample_size_sensitivity(expected, precision, z)


def likelihood_ratios(sensitivity, specificity):
    sens = min(max(float(sensitivity), 1e-9), 1 - 1e-9)
    spec = min(max(float(specificity), 1e-9), 1 - 1e-9)
    return {
        "lr_pos": round(sens / (1 - spec), 3),
        "lr_neg": round((1 - sens) / spec, 3),
    }


def cohen_kappa(a, b):
    """Iki kategori arasinda Cohen kappa (ornek: tasarim vs literatur)."""
    if len(a) != len(b) or not a:
        return None
    labels = sorted(set(a) | set(b))
    idx = {lab: i for i, lab in enumerate(labels)}
    n = len(a)
    mat = np.zeros((len(labels), len(labels)), dtype=float)
    for x, y in zip(a, b):
        mat[idx[x], idx[y]] += 1
    po = np.trace(mat) / n
    pe = float(np.sum(np.outer(mat.sum(axis=0), mat.sum(axis=1))) / (n * n))
    if pe >= 1:
        return 1.0
    return round(float((po - pe) / (1 - pe)), 4)
