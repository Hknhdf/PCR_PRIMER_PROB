"""
panel_opt.py
Multiplex panel uyumluluk ve kanal atama motoru:
  - assay'ler arasi heterodimer ΔG matrisi
  - 4 renkli kanal atamasi (amplikon boyu / prob Tm cakismasi)
  - panel uyarilari
"""

from __future__ import annotations

import itertools

import networkx as nx
import pandas as pd

from thermo import heterodimer

CHANNELS = ("FAM", "HEX", "ROX", "Cy5")
DIMER_WARN_DG = -6.0


def _assay_name(a, i):
    return str(a.get("name") or a.get("virus") or a.get("target") or f"Assay {i + 1}")


def _oligos(a):
    out = {}
    for key, label in (("forward", "F"), ("reverse", "R"), ("probe", "P")):
        seq = a.get(key)
        if seq:
            out[label] = str(seq).upper()
    return out


def panel_dimer_table(assays):
    """Tum assay'ler arasi oligo-oligo heterodimer etkilesimleri."""
    rows = []
    for i, a in enumerate(assays):
        for j, b in enumerate(assays):
            if j <= i:
                continue
            for (la, sa), (lb, sb) in itertools.product(_oligos(a).items(), _oligos(b).items()):
                r = heterodimer(sa, sb)
                if r["dg"] <= DIMER_WARN_DG:
                    rows.append({
                        "Assay 1": _assay_name(a, i),
                        "Oligo 1": la,
                        "Assay 2": _assay_name(b, j),
                        "Oligo 2": lb,
                        "ΔG (kcal/mol)": r["dg"],
                        "Tm (°C)": r["tm"],
                        "Risk": "yüksek" if r["dg"] <= -8 else "orta",
                    })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("ΔG (kcal/mol)")
    return df


def compatibility_matrix(assays):
    """Assay çiftleri arasindaki en kotu heterodimer ΔG matrisi."""
    names = [_assay_name(a, i) for i, a in enumerate(assays)]
    mat = pd.DataFrame(0.0, index=names, columns=names)
    for i, a in enumerate(assays):
        for j, b in enumerate(assays):
            if i == j:
                continue
            worst = 0.0
            for sa in _oligos(a).values():
                for sb in _oligos(b).values():
                    dg = heterodimer(sa, sb)["dg"]
                    if dg < worst:
                        worst = dg
            mat.iloc[i, j] = worst
    return mat


def _conflict(a, b):
    """Ayni kanalda karisir mi?"""
    amp_a = a.get("amplicon_size")
    amp_b = b.get("amplicon_size")
    if amp_a and amp_b and abs(float(amp_a) - float(amp_b)) <= 10:
        return True
    tm_a = a.get("probe_tm")
    tm_b = b.get("probe_tm")
    if tm_a and tm_b and abs(float(tm_a) - float(tm_b)) <= 2.0:
        return True
    return False


def channel_assignment(assays, channels=CHANNELS):
    """Graf renklendirme ile kanal atamasi; cakismalari minimuma indirir."""
    n = len(assays)
    if n == 0:
        return pd.DataFrame(), []
    names = [_assay_name(a, i) for i, a in enumerate(assays)]
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i, j in itertools.combinations(range(n), 2):
        if _conflict(assays[i], assays[j]):
            G.add_edge(i, j)
    coloring = nx.coloring.greedy_color(G, strategy="DSATUR")
    used = max(coloring.values()) + 1 if coloring else 0
    warnings = []
    if used > len(channels):
        warnings.append(f"{used} renk gerekiyor; yalnızca {len(channels)} kanal var. Bazı assay'ler aynı kanalı paylaşacak.")
    rows = []
    for i, a in enumerate(assays):
        color = coloring.get(i, 0)
        channel = channels[color] if color < len(channels) else f"Kanal {color + 1}"
        same_channel = [names[j] for j in G.neighbors(i) if coloring.get(j) == color]
        all_conflicts = [names[j] for j in G.neighbors(i)]
        rows.append({
            "Assay": names[i],
            "Kanal": channel,
            "Amplikon (bp)": a.get("amplicon_size"),
            "Prob Tm (°C)": a.get("probe_tm"),
            "Kanal çakışması": ", ".join(same_channel) if same_channel else "—",
            "Tüm çakışmalar": ", ".join(all_conflicts) if all_conflicts else "—",
        })
    return pd.DataFrame(rows), warnings


def panel_summary(assays):
    if not assays:
        return {}
    amps = [float(a.get("amplicon_size")) for a in assays if a.get("amplicon_size")]
    tms = [float(a.get("probe_tm")) for a in assays if a.get("probe_tm")]
    table, warnings = channel_assignment(assays)
    return {
        "n": len(assays),
        "amplicon_min": min(amps) if amps else None,
        "amplicon_max": max(amps) if amps else None,
        "probe_tm_min": min(tms) if tms else None,
        "probe_tm_max": max(tms) if tms else None,
        "channels": table["Kanal"].nunique() if not table.empty else 0,
        "warnings": warnings,
    }
