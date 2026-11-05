"""
v3_ui.py
v3.0 analiz / kütüphane / ileri sekmelerinin arayüzü.

Bu modül Streamlit'i doğrudan çalıştırmaz; realtime_pcr_designer.py içinden
bir bağlam (ctx) sözlüğü ile çağrılır. Ağır NCBI/BLAST işlemleri yalnızca
kullanıcı butona bastığında çalışır.
"""

from __future__ import annotations

import math
import statistics
from datetime import datetime

import pandas as pd
import streamlit as st

from advanced_design import (
    check_lamp_set,
    design_consensus_primers,
    design_lamp,
    design_rpa,
    design_snp_primers,
)
from browser_track import build_genome_track
from genotype import (
    analyze_genotypes,
    coverage_figure,
    heatmap_figure,
    logo_figure,
)
from literature import (
    compare_assays,
    comparison_table,
    library_dataframe,
    library_entry,
    library_from_json,
    library_to_json,
    locate_assay,
    validate_entry,
)
from panel_opt import (
    channel_assignment,
    compatibility_matrix,
    panel_dimer_table,
    panel_summary,
)
from specificity import blast_oligo, in_silico_pcr, risk_level, scan_sequences, specificity_score
from strain_check import (
    FAIL,
    build_order_gate,
    consensus_to_active,
    evaluate_turkey_strains,
    fetch_turkey_subjects,
    order_sheet_rows,
    parse_fasta_text,
    turkey_rows,
)
from stats_engine import (
    likelihood_ratios,
    lod_loq_analysis,
    ppv_npv,
    sample_size_sensitivity,
    sample_size_specificity,
)
from thermo import analyze_oligo, dimer_matrix, thermo_score, tm_uncertainty


# ---------------------------------------------------------------------------
# Ortak yardımcılar
# ---------------------------------------------------------------------------

def _active_header(active):
    st.markdown(f"### Aktif assay · {active.get('virus', '—')}")
    cols = st.columns(4)
    cols[0].metric("Forward Tm", f"{active.get('forward_tm') or '—'} °C")
    cols[1].metric("Reverse Tm", f"{active.get('reverse_tm') or '—'} °C")
    cols[2].metric("Prob Tm", f"{active.get('probe_tm') or '—'} °C")
    cols[3].metric("Amplikon", f"{active.get('amplicon_size') or '—'} bp")
    st.caption(
        f"{active.get('category', '—')} · {active.get('segment', '—')} · "
        f"{active.get('accession', '—')} · {active.get('region_text', '—')} · "
        f"F {active.get('forward_start')}–{active.get('forward_end')}, "
        f"R {active.get('reverse_start')}–{active.get('reverse_end')}"
    )
    with st.expander("Aktif primer / prob dizileri"):
        st.code(
            f"Forward 5'-{active.get('forward', '')}-3'\n"
            f"Reverse 5'-{active.get('reverse', '')}-3'\n"
            f"Prob    5'-{active.get('probe', '')}-3'",
            language="text",
        )


def _thermo_rows(active):
    rows = []
    for name, seq, role in (
        ("F", active.get("forward"), "primer"),
        ("R", active.get("reverse"), "primer"),
        ("P", active.get("probe"), "probe"),
    ):
        if not seq:
            continue
        a = analyze_oligo(seq, role=role)
        rows.append({
            "name": name,
            "tm": a["tm"],
            "hairpin_dg": a["hairpin_dg"],
            "homodimer_dg": a["homodimer_dg"],
            "three_prime_dg": a["three_prime_dg"],
            "warnings": a["warnings"],
            "analysis": a,
        })
    return rows


def _load_pathogen_segment(info, ctx):
    """Panel taraması için etkenin taranacak dizisini hazırlar."""
    has_gene = bool(info.get("gene_start") and info.get("gene_end"))
    if has_gene:
        fetched = ctx["fetch_genome"](info["genbank_id"])
        if fetched is None:
            return None
        seq, organism = fetched
        start = max(0, int(info["gene_start"]) - 1)
        end = min(len(seq), int(info["gene_end"]))
        return {
            "scan_seq": seq if len(seq) <= 200_000 else seq[start:end],
            "region_seq": seq[start:end],
            "organism": organism,
            "start": start,
            "end": end,
            "full_len": len(seq),
            "label": f"{start + 1}-{end}",
        }

    record = ctx["fetch_record"](info["genbank_id"])
    seq = ctx["record_sequence"](record)
    if seq is None:
        fetched = ctx["fetch_genome"](info["genbank_id"])
        if fetched is None:
            return None
        seq, organism = fetched
    else:
        organism = record.annotations.get("organism", "")
    start, end = ctx["parse_region"](info.get("target_region", ""), len(seq))
    return {
        "scan_seq": seq if len(seq) <= 200_000 else seq[start:end],
        "region_seq": seq[start:end],
        "organism": organism,
        "start": start,
        "end": end,
        "full_len": len(seq),
        "label": f"{start + 1}-{end}",
    }


def _coverage_rows(result):
    rows = []
    for key, c in (result.get("coverage") or {}).items():
        rows.append({
            "Oligo": key,
            "Değerlendirilebilir": f"{c['n_evaluable']}/{c['n_total']}",
            "Mükemmel (%)": c["perfect_pct"],
            "%95 GA": f"{c['perfect_ci'][0]}–{c['perfect_ci'][1]}",
            "Ortalama mismatch": c["mean_mm"],
            "3' uç riski": c["three_prime_risk"],
        })
    return pd.DataFrame(rows)


def _safe_float(value, default=0.0):
    try:
        return float(value)
    except Exception:
        return default


# ---------------------------------------------------------------------------
# ANALİZ
# ---------------------------------------------------------------------------

def render_analysis_tabs(ctx):
    st.markdown(
        "Aktif assay üzerinde **termodinamik**, **özgüllük**, **genotip kapsama**, "
        "**LOD/LOQ** ve **PDF rapor** analizleri. Ağır işlemler yalnızca butonla çalışır."
    )
    sub_thermo, sub_spec, sub_geno, sub_lod, sub_pdf = st.tabs(
        ["🌡️ Termodinamik", "🎯 Özgüllük", "🧬 Genotip Kapsama", "📈 LOD / LOQ", "📄 PDF Rapor"]
    )
    with sub_thermo:
        _render_thermo(ctx)
    with sub_spec:
        _render_specificity(ctx)
    with sub_geno:
        _render_genotype(ctx)
    with sub_lod:
        _render_lod(ctx)
    with sub_pdf:
        _render_pdf(ctx)


def _render_thermo(ctx):
    active = ctx["active_or_demo"]("thermo")
    if not active:
        return
    _active_header(active)
    st.markdown("#### Primer3 termodinamiği")
    st.caption(
        "Yapı analizi 60 °C'de, 3' uç stabilitesi 37 °C'de hesaplanır. "
        "ΔG değerleri kcal/mol'dür; daha negatif = daha güçlü istenmeyen yapı."
    )
    rows = _thermo_rows(active)
    if not rows:
        st.warning("Aktif assay'de oligonükleotid bulunamadı.")
        return

    a_cols = st.columns(len(rows))
    for col, row in zip(a_cols, rows):
        role = "probe" if row["name"] == "P" else "primer"
        score = thermo_score(row["analysis"], role=role)
        col.metric(f"{row['name']} termo skoru", f"{score}/100")
        col.caption(f"Tm {row['tm']} °C · GC {row['analysis']['gc']}%")

    table = []
    for row in rows:
        a = row["analysis"]
        table.append({
            "Oligo": row["name"],
            "Dizi": a["sequence"],
            "Uzunluk": a["length"],
            "Tm (°C)": a["tm"],
            "GC (%)": a["gc"],
            "Hairpin ΔG": a["hairpin_dg"],
            "Homodimer ΔG": a["homodimer_dg"],
            "3' uç ΔG": a["three_prime_dg"],
            "3' G/C": a["gc_clamp"],
            "Max run": a["max_run"],
            "Uyarı": "; ".join(a["warnings"])[:180] or "—",
        })
    st.dataframe(pd.DataFrame(table), width="stretch", hide_index=True)

    # Radar: açıklanabilir alt skorlar
    def _component_scores(a, role):
        if role == "primer":
            tm_ok = a["tm"] is not None and 58 <= a["tm"] <= 65
            gc_ok = 40 <= a["gc"] <= 60
        else:
            tm_ok = a["tm"] is not None and 66 <= a["tm"] <= 73
            gc_ok = 40 <= a["gc"] <= 65
        tm_score = 100 if tm_ok else max(0, 100 - 4 * abs((a["tm"] or 61) - (61.5 if role == "primer" else 69.5)))
        gc_score = 100 if gc_ok else max(0, 100 - 3 * abs(a["gc"] - (50 if role == "primer" else 52.5)))
        hp = 100 if a["hairpin_dg"] > -3 else max(0, 100 + 10 * (a["hairpin_dg"] + 3))
        hd = 100 if a["homodimer_dg"] > -6 else max(0, 100 + 10 * (a["homodimer_dg"] + 6))
        tp = 100 if a["three_prime_dg"] > -6 else max(0, 100 + 10 * (a["three_prime_dg"] + 6))
        run = 100 if a["max_run"] < 5 else max(0, 100 - 20 * (a["max_run"] - 4))
        return [tm_score, gc_score, hp, hd, tp, run]

    import plotly.graph_objects as go
    cats = ["Tm uyumu", "GC uyumu", "Hairpin", "Homodimer", "3' uç", "Homopolymer"]
    fig = go.Figure()
    for row in rows:
        role = "probe" if row["name"] == "P" else "primer"
        vals = _component_scores(row["analysis"], role)
        fig.add_trace(go.Scatterpolar(
            r=vals + [vals[0]], theta=cats + [cats[0]],
            fill="toself", name=row["name"], opacity=0.55,
        ))
    fig.update_layout(
        polar=dict(radialaxis=dict(visible=True, range=[0, 100])),
        height=360, margin=dict(l=20, r=20, t=30, b=20),
        legend=dict(orientation="h", y=1.08),
        title="Termodinamik alt skorları (100 = ideal)",
    )
    st.plotly_chart(fig, width="stretch")

    st.markdown("#### Oligo arası dimer matrisi")
    oligos = {r["name"]: r["analysis"]["sequence"] for r in rows}
    if len(oligos) >= 2:
        mat, flags = dimer_matrix(oligos)
        heat = go.Figure(go.Heatmap(
            z=mat.values.tolist(), x=list(mat.columns), y=list(mat.index),
            zmin=-12, zmax=0, colorscale=[[0, "#b71c1c"], [0.5, "#ffffff"], [1, "#e3f2fd"]],
            text=[[f"{v:.2f}" for v in row] for row in mat.values.tolist()],
            texttemplate="%{text}", hovertemplate="%{y}–%{x}<br>ΔG %{z:.2f} kcal/mol<extra></extra>",
        ))
        heat.update_layout(height=300, margin=dict(l=20, r=20, t=30, b=20),
                           title="Heterodimer / homodimer ΔG (60 °C)")
        st.plotly_chart(heat, width="stretch")
        if flags:
            st.warning(
                f"⚠️ {len(flags)} oligo çiftinde ΔG ≤ -6 kcal/mol. "
                "Sentez öncesi dimer analizini tekrarlayın."
            )
            st.dataframe(pd.DataFrame(flags), width="stretch", hide_index=True)
        else:
            st.success("✅ ΔG eşiğini aşan güçlü dimer yapısı yok.")
    else:
        st.info("Dimer matrisi için en az iki oligo gerekir.")

    with st.expander("🎲 Tm belirsizliği (tuz/Mg/dNTP Monte Carlo)"):
        n_mc = st.slider("Simülasyon sayısı", 200, 2000, 500, 100, key="tm_mc_n")
        if st.button("Tm belirsizliğini hesapla", key="tm_mc_run"):
            out = {r["name"]: tm_uncertainty(r["analysis"]["sequence"], n=n_mc) for r in rows}
            st.session_state["thermo_unc"] = out
        unc = st.session_state.get("thermo_unc")
        if unc:
            unc_rows = []
            for k, v in unc.items():
                unc_rows.append({
                    "Oligo": k, "n": v["n"], "Ortalama Tm": v["mean"],
                    "SD": v["sd"], "%2.5": v["p2.5"], "%97.5": v["p97.5"],
                })
            st.dataframe(pd.DataFrame(unc_rows), width="stretch", hide_index=True)
            st.caption("qPCR tuz/Mg/dNTP/oligo konsantrasyon aralıklarında Tm dağılımı.")

    st.markdown("#### 📄 PDF assay kartı")
    st.caption("Aşağıdaki PDF, mevcut termodinamik ve varsa önceki analiz özetlerini içerir.")
    last_g = st.session_state.get("last_genotype") or {}
    last_s = st.session_state.get("last_specificity") or {}
    last_l = st.session_state.get("last_lod") or {}
    spec_arg = None
    if last_s:
        spec_arg = {
            "score": last_s.get("score"),
            "off_targets": last_s.get("n_offtargets"),
            "notes": "; ".join(last_s.get("notes", []))[:300],
        }
    lod_arg = None
    if last_l:
        sc = last_l.get("standard_curve") or {}
        lod_arg = {
            "lod95": last_l.get("lod95"),
            "lod95_ci": last_l.get("lod95_ci"),
            "loq": last_l.get("loq"),
            "slope": sc.get("slope"),
            "efficiency": sc.get("efficiency"),
        }
    geno_arg = None
    if last_g:
        geno_arg = {
            "n_found": last_g.get("n_found"),
            "n_total": last_g.get("n_total"),
            "amplicon_found_pct": last_g.get("amplicon_found_pct"),
            "identity_mean": last_g.get("identity_mean"),
            "conservation_mean": last_g.get("conservation_mean"),
        }
    thermo_rows = [
        {k: v for k, v in r.items() if k != "analysis"} for r in rows
    ]
    pdf = ctx["assay_pdf_bytes"](
        active, thermo_rows=thermo_rows,
        coverage=last_g.get("coverage") if last_g else None,
        specificity=spec_arg, lod=lod_arg, genotype=geno_arg,
        notes="Bu kart, uygulamanın mevcut oturumundaki analizlerden üretilmiştir.",
    )
    st.download_button(
        "📥 PDF assay kartını indir",
        data=pdf,
        file_name=f"assay_karti_{active.get('accession', 'assay')}_{datetime.now():%Y%m%d_%H%M}.pdf",
        mime="application/pdf",
        key="thermo_pdf",
    )


def _render_specificity(ctx):
    active = ctx["active_or_demo"]("specificity")
    if not active:
        return
    _active_header(active)
    st.markdown("#### Deterministik in-silico PCR ve panel taraması")
    st.caption(
        "Seed-and-extend motoru ≤N mismatch ile primer/prob bağlanma bölgelerini arar. "
        "BLAST isteğe bağlıdır ve 30–150 s sürebilir."
    )
    c1, c2, c3 = st.columns(3)
    max_mm = c1.slider("İzin verilen mismatch", 0, 3, 3, 1, key="spec_mm")
    amp_default = int(active.get("amplicon_size") or 150)
    amp_min = c2.number_input("Min amplikon (bp)", min_value=20, max_value=5000,
                              value=max(20, amp_default - 50), step=10, key="spec_amp_min")
    amp_max = c3.number_input("Max amplikon (bp)", min_value=30, max_value=10000,
                              value=max(amp_default + 100, 200), step=10, key="spec_amp_max")
    if amp_max <= amp_min:
        st.error("Max amplikon, min amplikondan büyük olmalı.")
        return

    tcol1, tcol2 = st.columns([1, 2])
    if tcol1.button("🎯 Hedef dizide in-silico PCR", type="primary", key="spec_target_run",
                    width="stretch"):
        hits = in_silico_pcr(
            active["forward"], active["reverse"], active["ref_seq"],
            active.get("probe") or None, max_mismatches=max_mm,
            min_amp=int(amp_min), max_amp=int(amp_max),
        )
        st.session_state["spec_target"] = [h.as_dict() for h in hits]
    tcol2.caption(
        "Hedef referans: "
        + ("tam genom" if active.get("ref_is_full") else "hedef gen bölgesi")
        + f" ({len(active['ref_seq']):,} bp). Konumlar referans içindedir."
    )

    hits = st.session_state.get("spec_target")
    if hits is not None:
        if not hits:
            st.error("❌ Hedef dizide primer çifti amplikon oluşturmadı.")
        else:
            st.success(f"✅ Hedefte {len(hits)} amplikon bulundu.")
            st.dataframe(pd.DataFrame([{
                "Amplikon başlangıç": h["amplicon_start_1based"],
                "Amplikon bitiş": h["amplicon_end_1based"],
                "Boyut (bp)": h["size"],
                "Toplam mismatch": h["total_mismatches"],
                "3' uç mismatch": h["three_prime_mismatches"],
                "Prob var": "evet" if h.get("probe") else "hayır",
                "Prob mismatch": (h.get("probe") or {}).get("n_mismatches"),
            } for h in hits]), width="stretch", hide_index=True)

    st.markdown("#### 🧩 Panel taraması")
    other_names = [n for n in ctx["ordered_names"](ctx["VIRUS_DB"]) if n != active.get("virus")]
    selected = st.multiselect(
        "Taranacak diğer etkenler",
        other_names,
        default=[],
        format_func=lambda n: ctx["display_name"](ctx["VIRUS_DB"], n),
        key="spec_panel_names",
        help="Büyük bakteri/parazit genomlarında yalnızca katalogdaki hedef gen bölgesi taranır.",
    )
    if st.button("🔍 Seçili etkenlerde tara", key="spec_panel_run", width="stretch"):
        if not selected:
            st.warning("En az bir etken seçin.")
        else:
            target_hits = in_silico_pcr(
                active["forward"], active["reverse"], active["ref_seq"],
                active.get("probe") or None, max_mismatches=max_mm,
                min_amp=int(amp_min), max_amp=int(amp_max),
            )
            target_hit = target_hits[0] if target_hits else None
            sequences = {}
            warnings = []
            progress = st.progress(0.0, text="Diziler hazırlanıyor...")
            for i, name in enumerate(selected, start=1):
                progress.progress((i - 1) / len(selected), text=f"{name} hazırlanıyor...")
                info = ctx["VIRUS_DB"].get(name)
                if not info:
                    continue
                loaded = _load_pathogen_segment(info, ctx)
                if loaded is None:
                    warnings.append(f"{name}: dizi alınamadı")
                    continue
                sequences[name] = loaded["scan_seq"]
            progress.progress(1.0, text="Tarama tamamlandı")
            progress.empty()
            for w in warnings:
                st.warning("⚠️ " + w)
            if sequences:
                records = scan_sequences(
                    active["forward"], active["reverse"], active.get("probe") or None,
                    sequences, max_mismatches=max_mm,
                    min_amp=int(amp_min), max_amp=int(amp_max),
                )
                score, notes = specificity_score(target_hit, records)
                st.session_state["last_specificity"] = {
                    "records": records,
                    "score": score,
                    "notes": notes,
                    "n_offtargets": sum(1 for r in records if r["best"] is not None),
                }

    spec = st.session_state.get("last_specificity")
    if spec:
        st.metric("Özgüllük skoru", f"{spec['score']}/100")
        if spec["notes"]:
            for n in spec["notes"]:
                st.warning(n)
        else:
            st.success("✅ Hedef dışı amplikon saptanmadı.")
        rows = []
        for rec in spec["records"]:
            best = rec.get("best")
            rows.append({
                "Etken": rec["name"],
                "Dizi uzunluğu": rec["length"],
                "Hit sayısı": rec["n_hits"],
                "En iyi boyut": best.size if best else "—",
                "Toplam mismatch": best.total_mismatches if best else "—",
                "3' uç mismatch": best.three_prime_mismatches if best else "—",
                "Risk": risk_level(best, is_target=False),
            })
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
        with st.expander("Tüm off-target hit'ler"):
            detail = []
            for rec in spec["records"]:
                for h in rec.get("hits", []):
                    detail.append({
                        "Etken": rec["name"],
                        "Başlangıç": h.amplicon_start + 1,
                        "Bitiş": h.amplicon_end,
                        "Boyut": h.size,
                        "Mismatch": h.total_mismatches,
                        "3' mm": h.three_prime_mismatches,
                    })
            st.dataframe(pd.DataFrame(detail), width="stretch", hide_index=True)

    st.markdown("#### 🌐 NCBI BLAST (isteğe bağlı)")
    b1, b2, b3 = st.columns([1, 1, 2])
    oligo_name = b1.selectbox("Oligo", ["F", "R", "P"], key="blast_oligo_name")
    database = b2.selectbox("Veritabanı", ["nt", "refseq_rna"], key="blast_db")
    human_only = b3.checkbox("Yalnızca insan (txid9606)", key="blast_human")
    sequence = active.get({"F": "forward", "R": "reverse", "P": "probe"}[oligo_name], "")
    if st.button("🌐 BLAST ile doğrula", key="blast_run"):
        if not sequence:
            st.warning("Seçili oligo yok.")
        else:
            with st.spinner("NCBI BLAST çalışıyor; 30–150 s sürebilir..."):
                out = blast_oligo(sequence, database=database, human_only=human_only)
            st.session_state["blast_last"] = {"oligo": oligo_name, "result": out}
    blast = st.session_state.get("blast_last")
    if blast:
        out = blast["result"]
        if out.get("error"):
            st.error(out["error"])
        else:
            st.success(f"BLAST tamamlandı · RID {out.get('rid', '—')}")
            if out.get("hits"):
                st.dataframe(pd.DataFrame([{
                    "Accession": h["accession"], "Tanım": h["title"],
                    "E-value": h["evalue"], "Identite (%)": h["identity"],
                    "Hizalama": h["align_length"], "Mismatch": h["mismatches"],
                } for h in out["hits"]]), width="stretch", hide_index=True)
            else:
                st.info("BLAST hit'i döndürmedi.")


def _render_genotype(ctx):
    active = ctx["active_or_demo"]("genotype")
    if not active:
        return
    _active_header(active)
    st.markdown("#### Çoklu dizi kapsama analizi")
    st.caption(
        "NCBI'dan aynı organizmaya ait diziler çekilir, referans amplikona hizalanır ve "
        "her primer/prob için mismatch, 3' uç riski ve Wilson %95 GA hesaplanır."
    )
    c1, c2, c3 = st.columns(3)
    n_seq = c1.slider("Dizi sayısı", 3, 20, 10, 1, key="geno_n")
    max_mm = c2.slider("Maksimum mismatch", 2, 6, 4, 1, key="geno_mm")
    include_ref = c3.checkbox("Referans satırını matrise ekle", value=True, key="geno_ref")
    if st.button("🧬 Dizileri çek ve genotip kapsamayı hesapla", type="primary",
                 key="geno_run", width="stretch"):
        organism = active.get("organism") or ""
        if not organism:
            st.error("Organizma adı yok. NCBI modunda tasarım yapın veya manuel modda organizma girin.")
            return
        big = (active.get("seq_len", 0) > 200_000) or (not active.get("ref_is_full"))
        keywords = ()
        if active.get("info", {}).get("target_gene"):
            keywords = (active["info"]["target_gene"],)
        with st.spinner("İlgili diziler çekiliyor ve hizalanıyor..."):
            subjects = ctx["get_genotype_subjects"](
                organism, active.get("accession", ""), active["ref_seq"],
                n_seq, keywords, big,
            )
        if not subjects:
            st.warning("NCBI'dan hizalanabilir dizi bulunamadı; yalnızca referans gösterilecek.")
            subjects = []
        with st.spinner("Mismatch matrisi ve kapsama hesaplanıyor..."):
            result = analyze_genotypes(
                active["ref_seq"], active["amplicon_start"], active["amplicon_end"],
                active["forward"], active["reverse"], active.get("probe") or None,
                subjects, probe_offset=active.get("probe_offset"),
                max_mismatches=max_mm, include_reference=include_ref,
            )
        st.session_state["last_genotype"] = result

    result = st.session_state.get("last_genotype")
    if not result or result.get("error"):
        if result and result.get("error"):
            st.error(result["error"])
        return

    m = st.columns(5)
    m[0].metric("Amplikon bulunan", f"{result['n_found']}/{result['n_total']}")
    m[1].metric("Bulunma oranı", f"%{result['amplicon_found_pct']}")
    m[2].metric("Ortalama identite", result.get("identity_mean") or "—")
    m[3].metric("Ortalama korunmuşluk", result.get("conservation_mean") or "—")
    m[4].metric("Ortalama entropi", result.get("entropy_mean") or "—")

    st.markdown("#### Primer / prob kapsaması")
    st.dataframe(_coverage_rows(result), width="stretch", hide_index=True)
    st.caption(
        "Mükemmel = 0 mismatch. Değerlendirilebilir = dizide bağlanma bölgesi bulunanlar. "
        "%95 GA Wilson skor aralığıdır."
    )

    fig = heatmap_figure(result)
    if fig:
        st.plotly_chart(fig, width="stretch")
    cov_fig = coverage_figure(result)
    if cov_fig:
        st.plotly_chart(cov_fig, width="stretch")

    logo = logo_figure(result)
    if logo:
        st.markdown("#### Sequence logo")
        try:
            import matplotlib.pyplot as plt
            st.pyplot(logo)
            plt.close(logo)
        except Exception:
            st.info("Logo görseli çizilemedi.")

    with st.expander("Ham mismatch matrisini indir"):
        labels = result.get("matrix_labels") or []
        positions = list(range(result["amplicon_start"] + 1, result["amplicon_end"] + 1))
        mat_df = pd.DataFrame(result["matrix"], index=labels, columns=positions)
        st.dataframe(mat_df.iloc[:, :60], width="stretch")
        st.download_button(
            "📥 CSV indir",
            data=mat_df.to_csv().encode("utf-8"),
            file_name=f"genotip_matrisi_{active.get('accession', 'assay')}.csv",
            mime="text/csv",
            key="geno_csv",
        )

    found_rows = []
    for rec in result.get("records", []):
        found_rows.append({
            "Dizi": rec.get("id"),
            "Açıklama": rec.get("description"),
            "Amplikon": "var" if rec.get("found") else ("kısmi" if rec.get("partial") else "yok"),
            "Boyut": rec.get("amplicon_size"),
            "Identite": rec.get("identity"),
            "F mm": (rec.get("F") or {}).get("n_mismatches"),
            "R mm": (rec.get("R") or {}).get("n_mismatches"),
            "P mm": (rec.get("P") or {}).get("n_mismatches"),
            "3' mm": rec.get("three_prime_mm"),
            "Hata": rec.get("error") or "",
        })
    if found_rows:
        st.markdown("#### Dizi bazında sonuçlar")
        st.dataframe(pd.DataFrame(found_rows), width="stretch", hide_index=True)


def _render_lod(ctx):
    st.markdown("### LOD / LOQ ve biyoistatistik")
    active = st.session_state.get("active_assay")
    if active:
        st.caption(
            f"Aktif assay: {active.get('virus')} · {active.get('amplicon_size')} bp. "
            "Simülasyon, seçilen amplikon boyutundan bağımsızdır."
        )
    st.caption(
        "Poisson kopya dağılımı + verim değişkenliği + Ct gürültüsü ile Monte Carlo; "
        "LOD95 probit MLE ve bootstrap ile hesaplanır. Bu bir **simülasyondur**, ıslak "
        "laboratuvar LOD/LOQ validasyonunun yerine geçmez."
    )

    with st.expander("Simülasyon parametreleri", expanded=True):
        grid_text = st.text_input(
            "Kopya düzeyleri (virgülle)",
            value="1,2,3,5,10,25,50,100,1000,10000",
            key="lod_grid",
        )
        c1, c2, c3, c4 = st.columns(4)
        replicates = c1.slider("Tekrar sayısı", 12, 96, 24, 4, key="lod_rep")
        efficiency = c2.slider("PCR verimi", 0.70, 1.00, 0.90, 0.01, key="lod_eff")
        eff_sd = c3.slider("Verim SD", 0.00, 0.15, 0.05, 0.01, key="lod_eff_sd")
        ct_sd = c4.slider("Ct SD", 0.00, 0.60, 0.30, 0.05, key="lod_ct_sd")
        c5, c6, c7, c8 = st.columns(4)
        inhib = c5.slider("İnhibisyon olasılığı", 0.00, 0.20, 0.00, 0.01, key="lod_inhib")
        max_cycles = c6.slider("Maksimum döngü", 30, 50, 40, 1, key="lod_cycles")
        det_target = c7.slider("LOD hedef tespit oranı", 0.90, 0.99, 0.95, 0.01, key="lod_det")
        loq_cv = c8.slider("LOQ Ct CV eşiği", 0.02, 0.15, 0.05, 0.01, key="lod_loq_cv")

    if st.button("📈 Monte Carlo simülasyonunu çalıştır", type="primary", key="lod_run",
                 width="stretch"):
        copies = []
        for part in grid_text.replace(";", ",").split(","):
            part = part.strip()
            if not part:
                continue
            try:
                v = float(part)
                if v > 0:
                    copies.append(v)
            except ValueError:
                continue
        copies = sorted(set(copies))
        if len(copies) < 4:
            st.error("En az 4 geçerli kopya düzeyi girin.")
        else:
            with st.spinner("Simülasyon ve probit bootstrap çalışıyor..."):
                result = lod_loq_analysis(
                    copies, replicates=replicates, efficiency=efficiency,
                    efficiency_sd=eff_sd, max_cycles=max_cycles,
                    ct_sd=ct_sd, inhibition_prob=inhib, seed=42, n_boot=200,
                    detection_target=det_target, loq_cv_target=loq_cv,
                )
            st.session_state["last_lod"] = result

    result = st.session_state.get("last_lod")
    if not result:
        st.info("Parametreleri ayarlayıp simülasyonu çalıştırın.")
        return

    sc = result.get("standard_curve") or {}
    m = st.columns(5)
    m[0].metric("LOD95 (kopya)", result.get("lod95") or "—")
    ci = result.get("lod95_ci")
    m[1].metric("LOD95 %95 GA", f"{ci['low']}–{ci['high']}" if ci else "—")
    m[2].metric("LOQ (kopya)", result.get("loq") or "—")
    m[3].metric("Standart eğri eğimi", sc.get("slope") or "—")
    m[4].metric("PCR verimi", f"{100 * sc.get('efficiency', 0):.1f}%" if sc.get("efficiency") else "—")
    st.caption(
        f"R² = {sc.get('r2', '—')} · n = {sc.get('n', '—')} · "
        f"LOD95 probit MLE; LOQ tespit ≥ %{100 * det_target:.0f} ve Ct CV ≤ {loq_cv:.2f}."
    )

    summary = result.get("summary") or []
    if summary:
        st.markdown("#### Tespit olasılığı")
        import plotly.graph_objects as go
        xs = [math.log10(r["copies"]) for r in summary]
        ys = [r["detection_rate"] for r in summary]
        lo = [r["detection_rate"] - r["detection_ci_low"] for r in summary]
        hi = [r["detection_ci_high"] - r["detection_rate"] for r in summary]
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode="markers+lines", name="Gözlenen",
            error_y=dict(type="data", symmetric=False, array=hi, arrayminus=lo),
            marker=dict(size=8, color="#1565c0"),
        ))
        fit = result.get("probit") or {}
        if fit.get("ok"):
            grid = [min(xs) - 0.2 + i * (max(xs) - min(xs) + 0.4) / 80 for i in range(81)]
            probs = [statistics.NormalDist().cdf(fit["a"] + fit["b"] * x) for x in grid]
            fig.add_trace(go.Scatter(
                x=grid, y=probs, mode="lines", name="Probit model",
                line=dict(color="#c62828", width=2),
            ))
            if result.get("lod95"):
                fig.add_vline(x=math.log10(result["lod95"]), line_dash="dash",
                              line_color="#2e7d32", annotation_text="LOD95")
        fig.update_layout(
            height=380, template="plotly_white",
            xaxis_title="log10(kopya)", yaxis_title="Tespit olasılığı",
            yaxis=dict(range=[0, 1.05]), legend=dict(orientation="h", y=1.1),
            margin=dict(l=20, r=20, t=30, b=40),
        )
        st.plotly_chart(fig, width="stretch")

        st.markdown("#### Ct – kopya ilişkisi ve standart eğri")
        rows = result.get("rows") or []
        det = [r for r in rows if r.get("ct") is not None]
        fig2 = go.Figure()
        fig2.add_trace(go.Scatter(
            x=[math.log10(r["copies"]) for r in det],
            y=[r["ct"] for r in det],
            mode="markers", name="Ct",
            marker=dict(size=5, color="#6a1b9a", opacity=0.65),
        ))
        if sc.get("ok"):
            xmin = min(math.log10(r["copies"]) for r in det)
            xmax = max(math.log10(r["copies"]) for r in det)
            fig2.add_trace(go.Scatter(
                x=[xmin, xmax],
                y=[sc["intercept"] + sc["slope"] * xmin, sc["intercept"] + sc["slope"] * xmax],
                mode="lines", name="Regresyon", line=dict(color="#c62828", width=2),
            ))
        fig2.update_layout(
            height=380, template="plotly_white",
            xaxis_title="log10(kopya)", yaxis_title="Ct",
            legend=dict(orientation="h", y=1.1),
            margin=dict(l=20, r=20, t=30, b=40),
        )
        st.plotly_chart(fig2, width="stretch")

        st.markdown("#### Özet tablo")
        st.dataframe(pd.DataFrame(summary), width="stretch", hide_index=True)

    st.markdown("#### PPV / NPV ve örneklem büyüklüğü")
    p1, p2, p3 = st.columns(3)
    sens = p1.slider("Duyarlılık", 0.50, 1.00, 0.95, 0.01, key="ppv_sens")
    spec_v = p2.slider("Özgüllük", 0.50, 1.00, 0.98, 0.01, key="ppv_spec")
    prev = p3.slider("Prevalans", 0.001, 0.50, 0.01, 0.001, key="ppv_prev", format="%.3f")
    metrics = ppv_npv(sens, spec_v, prev)
    lr = likelihood_ratios(sens, spec_v)
    q = st.columns(4)
    q[0].metric("PPV", f"%{100 * metrics['ppv']:.1f}")
    q[1].metric("NPV", f"%{100 * metrics['npv']:.1f}")
    q[2].metric("LR+", lr["lr_pos"])
    q[3].metric("LR−", lr["lr_neg"])

    s1, s2 = st.columns(2)
    exp_sens = s1.slider("Beklenen duyarlılık", 0.80, 0.99, 0.95, 0.01, key="ss_sens")
    prec = s2.slider("Hedef %95 GA yarı genişliği", 0.01, 0.10, 0.05, 0.01, key="ss_prec")
    st.caption(
        f"Duyarlılık için yaklaşık n = **{sample_size_sensitivity(exp_sens, prec)}**; "
        f"özgüllük için yaklaşık n = **{sample_size_specificity(0.98, prec)}**."
    )

    st.download_button(
        "📥 LOD/LOQ özetini indir (CSV)",
        data=pd.DataFrame(summary).to_csv(index=False).encode("utf-8") if summary else b"",
        file_name=f"lod_loq_{datetime.now():%Y%m%d_%H%M}.csv",
        mime="text/csv",
        key="lod_csv",
    )


def _render_pdf(ctx):
    active = ctx["active_or_demo"]("pdf")
    if not active:
        return
    _active_header(active)
    st.markdown("### PDF assay kartı")
    st.caption(
        "Termodinamik, varsa genotip kapsama, özgüllük ve LOD/LOQ özetleri tek PDF'te birleştirilir."
    )
    thermo_rows = [
        {k: v for k, v in r.items() if k != "analysis"} for r in _thermo_rows(active)
    ]
    last_g = st.session_state.get("last_genotype") or {}
    last_s = st.session_state.get("last_specificity") or {}
    last_l = st.session_state.get("last_lod") or {}
    spec_arg = None
    if last_s:
        spec_arg = {
            "score": last_s.get("score"),
            "off_targets": last_s.get("n_offtargets"),
            "notes": "; ".join(last_s.get("notes", []))[:300],
        }
    lod_arg = None
    if last_l:
        sc = last_l.get("standard_curve") or {}
        lod_arg = {
            "lod95": last_l.get("lod95"), "lod95_ci": last_l.get("lod95_ci"),
            "loq": last_l.get("loq"), "slope": sc.get("slope"),
            "efficiency": sc.get("efficiency"),
        }
    geno_arg = None
    if last_g:
        geno_arg = {
            "n_found": last_g.get("n_found"), "n_total": last_g.get("n_total"),
            "amplicon_found_pct": last_g.get("amplicon_found_pct"),
            "identity_mean": last_g.get("identity_mean"),
            "conservation_mean": last_g.get("conservation_mean"),
        }
    notes = st.text_area("Rapor notu", key="pdf_notes",
                         placeholder="Laboratuvar notu, hedef bölge, kullanım amacı...")
    blocks = [
        "Termodinamik: " + ("var" if thermo_rows else "yok"),
        "Genotip kapsama: " + ("var" if last_g else "yok"),
        "Özgüllük: " + ("var" if last_s else "yok"),
        "LOD/LOQ: " + ("var" if last_l else "yok"),
    ]
    st.caption(" · ".join(blocks))
    pdf = ctx["assay_pdf_bytes"](
        active, thermo_rows=thermo_rows,
        coverage=last_g.get("coverage") if last_g else None,
        specificity=spec_arg, lod=lod_arg, genotype=geno_arg, notes=notes,
    )
    st.download_button(
        "📥 PDF raporunu indir",
        data=pdf,
        file_name=f"assay_raporu_{active.get('accession', 'assay')}_{datetime.now():%Y%m%d_%H%M}.pdf",
        mime="application/pdf",
        key="pdf_report",
        width="stretch",
    )
    st.info(
        "PDF araştırma amaçlıdır. Klinik kullanım öncesi Primer-BLAST, dimer analizi, "
        "LOD/LOQ ve ıslak laboratuvar doğrulaması gereklidir."
    )


# ---------------------------------------------------------------------------
# SİPARİŞ KONTROLÜ
# ---------------------------------------------------------------------------

def _status_text(status):
    return {
        "PASS": "✅ PASS",
        "WARN": "⚠️ WARN",
        "FAIL": "❌ FAIL",
        "BLOCKED": "⛔ BLOCKED",
    }.get(status, str(status))


def render_order_tabs(ctx):
    st.markdown(
        "Türkiye'de dolaşan suşlara uyumu kontrol edin ve sipariş öncesi son QC "
        "kapısını çalıştırın. **Bu bölüm in-silico ön kontroldür; ıslak laboratuvar "
        "validasyonunun yerine geçmez.**"
    )
    sub_tr, sub_qc = st.tabs(["🇹🇷 Türkiye Suş Kontrolü", "✅ Sipariş Öncesi Kontrol"])
    with sub_tr:
        _render_turkey_strains(ctx)
    with sub_qc:
        _render_order_gate(ctx)


def _dedupe_subjects(subjects):
    seen = set()
    unique = []
    for item in subjects:
        seq = str(item.get("seq", "")).upper()
        if not seq or seq in seen:
            continue
        seen.add(seq)
        unique.append(item)
    return unique


def _render_turkey_strains(ctx):
    active = ctx["active_or_demo"]("turkey")
    if not active:
        return
    _active_header(active)
    st.markdown("#### Türkiye suşu uyum kontrolü")
    st.caption(
        "NCBI'da **Turkey[Country]** filtresiyle bulunan diziler ve/veya yüklediğiniz "
        "yerel FASTA dizileri referans amplikona hizalanır. Primer/prob bağlanma "
        "bölgelerindeki mismatch ve 3' uç riski raporlanır. Ülke metadata'sı olmayan "
        "kayıtlar NCBI aramasına gelmez; bu durumda FASTA yükleyin."
    )

    c1, c2, c3 = st.columns(3)
    n_seq = c1.slider("NCBI Türkiye dizisi", 5, 40, 20, 1, key="tr_n")
    max_mm = c2.slider("İzin verilen mismatch", 2, 6, 4, 1, key="tr_mm")
    country = c3.text_input("NCBI ülke filtresi", value="Turkey", key="tr_country",
                            help="NCBI'da ülke alanı genellikle 'Turkey' olarak geçer.")
    include_ncbi = st.checkbox(
        "NCBI'da Turkey[Country] aramasını dahil et", value=True, key="tr_include_ncbi"
    )

    uploaded = st.file_uploader(
        "Türkiye suşu FASTA (isteğe bağlı, çoklu kayıt olabilir)",
        type=["fasta", "fa", "fna", "txt"], key="tr_fasta",
        help="Kendi dizinizi veya laboratuvarınızın sekansını yükleyin. FASTA başlıkları accession/izolat adı olarak kullanılır.",
    )
    manual = st.text_area(
        "veya FASTA / ham dizi yapıştırın",
        key="tr_manual", height=130,
        placeholder=">TR_isolate_1\nACGT...\n>TR_isolate_2\nACGT...",
    )

    if st.button("🇹🇷 Türkiye dizilerini değerlendir", type="primary",
                 key="tr_run", width="stretch"):
        subjects = []
        if manual.strip():
            subjects.extend(parse_fasta_text(manual))
        if uploaded is not None:
            try:
                subjects.extend(parse_fasta_text(uploaded.getvalue().decode("utf-8", errors="replace")))
            except Exception as exc:
                st.warning(f"Yüklenen FASTA okunamadı: {exc}")
        if include_ncbi:
            organism = active.get("organism") or ""
            if not organism:
                st.warning("Organizma adı yok; NCBI araması atlandı. Manuel FASTA kullanabilirsiniz.")
            else:
                ref_amp = active["ref_seq"][active["amplicon_start"]:active["amplicon_end"]]
                keywords = ()
                if active.get("info", {}).get("target_gene"):
                    keywords = (active["info"]["target_gene"],)
                with st.spinner("NCBI'da Türkiye dizileri aranıyor..."):
                    try:
                        ncbi_subjects = fetch_turkey_subjects(
                            organism, ref_amp, n=n_seq, keywords=keywords,
                            accession=active.get("accession", ""), country=country,
                        )
                    except Exception as exc:
                        ncbi_subjects = []
                        st.warning(f"NCBI araması başarısız: {exc}")
                if not ncbi_subjects:
                    st.warning(
                        "NCBI'da bu organizma + ülke filtresiyle dizi bulunamadı. "
                        "Ülke metadata'sı eksik olabilir; yerel FASTA yükleyin."
                    )
                subjects.extend(ncbi_subjects)
        subjects = _dedupe_subjects(subjects)
        if not subjects:
            st.error("Değerlendirilecek dizi yok. NCBI aramasını açın veya FASTA yükleyin.")
        else:
            with st.spinner("Türkiye dizileri referans amplikona hizalanıyor..."):
                result = evaluate_turkey_strains(active, subjects, max_mismatches=max_mm)
            st.session_state["last_turkey"] = result
            st.session_state["last_turkey_subjects"] = subjects
            st.session_state["last_turkey_consensus"] = None

    result = st.session_state.get("last_turkey")
    if not result:
        st.info("Henüz Türkiye suşu kontrolü çalıştırılmadı.")
        return
    if result.get("error"):
        st.error(result["error"])
        return

    m = st.columns(5)
    m[0].metric("Değerlendirilen", result.get("n_total", 0))
    m[1].metric("PASS", result.get("n_pass", 0))
    m[2].metric("WARN", result.get("n_warn", 0))
    m[3].metric("FAIL", result.get("n_fail", 0))
    m[4].metric("3' uç riski", result.get("n_three_prime_risk", 0))

    if result.get("overall") == "PASS":
        st.success("✅ Tüm değerlendirilen Türkiye dizileri primer/prob setiyle uyumlu.")
    elif result.get("overall") == "WARN":
        st.warning("⚠️ Bazı dizilerde mismatch/kısmi eşleşme var; sipariş öncesi inceleyin.")
    else:
        st.error("❌ En az bir Türkiye dizisinde kritik uyumsuzluk var; sipariş vermeyin.")

    st.dataframe(pd.DataFrame(turkey_rows(result)), width="stretch", hide_index=True)
    st.caption(
        "PASS = tam eşleşme; WARN = 1–2 mismatch veya kısmi eşleşme; "
        "FAIL = 3' uç mismatch, ≥3 primer mismatch veya ≥2 prob mismatch."
    )

    if result.get("coverage"):
        st.markdown("#### Primer / prob kapsaması")
        st.dataframe(_coverage_rows(result), width="stretch", hide_index=True)

    subjects = st.session_state.get("last_turkey_subjects") or []
    if subjects:
        st.markdown("#### Türkiye dizileriyle yeni aday tasarımı")
        st.caption(
            "Mevcut set Türkiye dizilerinde FAIL/WARN veriyorsa, bu dizilerin majority/consensus "
            "dizisinden yeni primer-prob adayları üretilebilir. Adaylar **in-silico**dir; "
            "sentez öncesi Primer-BLAST ve ıslak laboratuvar doğrulaması gerekir."
        )
        cc1, cc2, cc3 = st.columns(3)
        tr_cov = cc1.slider("Minimum Türkiye kapsaması (%)", 50, 100, 80, 5, key="tr_cons_cov") / 100.0
        tr_deg = cc2.selectbox("Maksimum dejenere kat sayısı", [64, 256, 1024], index=1, key="tr_cons_deg")
        tr_region = cc3.selectbox(
            "Tasarım bölgesi",
            ["Aktif amplikon", "Tüm referans dizisi"],
            key="tr_cons_region",
            help="Aktif amplikon dışında daha korunmuş bir bölge aramak için tüm referansı deneyin.",
        )
        if st.button("🇹🇷 Türkiye uyumlu adayları üret", key="tr_cons_run", width="stretch"):
            if tr_region == "Aktif amplikon":
                region_start = int(active["amplicon_start"])
                region_end = int(active["amplicon_end"])
            else:
                region_start, region_end = 0, len(active["ref_seq"])
            ref_region = active["ref_seq"][region_start:region_end]
            with st.spinner("Türkiye dizilerinden consensus adaylar hesaplanıyor..."):
                cons = design_consensus_primers(
                    ref_region, subjects, region_start=0,
                    target_amp=int(active.get("amplicon_size") or 150),
                    top_n=20, min_coverage=tr_cov, max_degeneracy=int(tr_deg),
                )
            if cons.get("ok") and cons.get("candidates"):
                base = dict(active)
                base["amplicon_start"] = region_start
                base["amplicon_end"] = region_end
                evaluated = []
                for c in cons["candidates"]:
                    try:
                        cand_active = consensus_to_active(base, c)
                        tr = evaluate_turkey_strains(cand_active, subjects, max_mismatches=max_mm)
                        c = dict(c)
                        c["turkey_status"] = tr.get("overall", FAIL)
                        c["turkey_pass"] = tr.get("n_pass", 0)
                        c["turkey_warn"] = tr.get("n_warn", 0)
                        c["turkey_fail"] = tr.get("n_fail", 0)
                        c["turkey_three"] = tr.get("n_three_prime_risk", 0)
                        evaluated.append(c)
                    except Exception:
                        continue
                cons["evaluated_candidates"] = evaluated
                cons["region_start"] = region_start
                cons["region_end"] = region_end
            st.session_state["last_turkey_consensus"] = cons

        cons = st.session_state.get("last_turkey_consensus")
        if cons:
            if not cons.get("ok"):
                st.error(cons.get("message") or "Consensus adayı üretilemedi.")
            elif not cons.get("candidates"):
                st.warning("Bu kapsama/dejenerasyon eşiğinde aday bulunamadı.")
            else:
                evaluated = cons.get("evaluated_candidates") or []
                candidate_rows = [{
                    "F primer": c["forward"],
                    "R primer": c["reverse"],
                    "Prob": c["probe"],
                    "F kapsama (%)": c["coverage_f"],
                    "R kapsama (%)": c["coverage_r"],
                    "P kapsama (%)": c["coverage_p"],
                    "Boyut (bp)": c["amplicon_size"],
                    "F Tm": c["forward_tm"],
                    "R Tm": c["reverse_tm"],
                    "P Tm": c["probe_tm"],
                    "Türkiye": _status_text(c.get("turkey_status", FAIL)),
                    "TR PASS/WARN/FAIL": f"{c.get('turkey_pass', 0)}/{c.get('turkey_warn', 0)}/{c.get('turkey_fail', 0)}",
                    "Dejenerasyon F/R/P": f"{c['forward_deg']}/{c['reverse_deg']}/{c['probe_deg']}",
                } for c in evaluated]
                if candidate_rows:
                    st.dataframe(pd.DataFrame(candidate_rows), width="stretch", hide_index=True)
                with st.expander("Dejenere (IUPAC) karşılıkları"):
                    deg_rows = [{
                        "F primer (dejenere)": c["forward_degenerate"],
                        "R primer (dejenere)": c["reverse_degenerate"],
                        "Prob (dejenere)": c["probe_degenerate"],
                    } for c in evaluated]
                    st.dataframe(pd.DataFrame(deg_rows), width="stretch", hide_index=True)

                usable = [
                    c for c in evaluated
                    if c.get("turkey_fail", 1) == 0 and c.get("turkey_three", 1) == 0
                ]
                if not usable:
                    st.error(
                        "Üretilen non-dejenere adayların hiçbiri bu Türkiye dizilerinin tamamında "
                        "3' uç/prob kriterlerini geçmedi. Adayı sipariş etmeyin; dejenere IUPAC "
                        "tasarımını veya farklı hedef bölgeyi deneyin."
                    )
                else:
                    labels = [
                        f"{i + 1}. {c['forward']} / {c['reverse']} / {c['probe']}"
                        for i, c in enumerate(usable)
                    ]
                    idx = st.selectbox("Aktif assay yapılacak aday", range(len(labels)),
                                       format_func=lambda i: labels[i], key="tr_cons_choice")
                    if st.button("Seçili Türkiye adayını aktif assay yap", key="tr_cons_use",
                                 width="stretch"):
                        try:
                            new_active = consensus_to_active(
                                {**active,
                                 "amplicon_start": cons.get("region_start", active["amplicon_start"]),
                                 "amplicon_end": cons.get("region_end", active["amplicon_end"])},
                                usable[idx],
                            )
                            st.session_state["active_assay"] = new_active
                            for _k in ("last_turkey", "last_turkey_subjects", "last_turkey_consensus",
                                       "last_genotype", "last_specificity", "last_lod",
                                       "last_order_gate", "blast_last",
                                       "order_blast", "order_synthesis", "order_wetlab", "order_ack"):
                                st.session_state[_k] = False if _k.startswith("order_") else None
                            st.success("Türkiye consensus adayı aktif assay yapıldı. "
                                       "Türkiye kontrolünü ve sipariş QC'yi yeniden çalıştırın.")
                            st.rerun()
                        except Exception as exc:
                            st.error(f"Aday aktif assay'a çevrilemedi: {exc}")

    with st.expander("Ham Türkiye sonuçlarını indir"):
        st.download_button(
            "📥 Türkiye kontrol CSV'sini indir",
            data=pd.DataFrame(turkey_rows(result)).to_csv(index=False).encode("utf-8-sig"),
            file_name=f"turkiye_sus_kontrolu_{active.get('accession', 'assay')}.csv",
            mime="text/csv",
            key="tr_download",
        )


def _render_order_gate(ctx):
    active = ctx["active_or_demo"]("order")
    if not active:
        return
    _active_header(active)
    st.markdown("#### Sipariş öncesi QC kapısı")
    st.caption(
        "Termodinamik, özgüllük, Türkiye suşları ve manuel teyitler birlikte değerlendirilir. "
        "FAIL veya BLOCKED varsa sipariş listesi kilitlenir."
    )

    turkey_result = st.session_state.get("last_turkey")
    specificity = st.session_state.get("last_specificity")
    c1, c2, c3 = st.columns(3)
    primer_blast_done = c1.checkbox("Primer-BLAST / BLASTn yapıldı", key="order_blast")
    synthesis_confirmed = c2.checkbox("Sentez/modifikasyon teyit edildi", key="order_synthesis")
    wetlab_planned = c3.checkbox("LOD/LOQ ve ıslak lab planı hazır", key="order_wetlab")

    o1, o2, o3, o4 = st.columns(4)
    reporter = o1.selectbox("Prob 5' raportör", ["FAM", "HEX", "ROX", "Cy5", "ATT O550", "Diğer"], key="order_reporter")
    quencher = o2.selectbox("Prob 3' quencher", ["BHQ1", "BHQ2", "BHQ3", "TAMRA", "Diğer"], key="order_quencher")
    scale = o3.selectbox("Sentez ölçeği", ["0.05 µmol", "0.1 µmol", "0.25 µmol", "1 µmol"], index=1, key="order_scale")
    purification = o4.selectbox(
        "Saflaştırma",
        ["Standart desalt (primer) / HPLC (prob)", "HPLC", "PAGE"],
        key="order_purification",
    )

    gate = build_order_gate(
        active,
        turkey_result=turkey_result,
        specificity=specificity,
        primer_blast_done=primer_blast_done,
        synthesis_confirmed=synthesis_confirmed,
        wetlab_planned=wetlab_planned,
    )

    m = st.columns(3)
    m[0].metric("Toplam PASS", gate.get("n_pass", 0))
    m[1].metric("WARN", gate.get("n_warn", 0))
    m[2].metric("FAIL", gate.get("n_fail", 0))

    if gate["overall"] == "PASS":
        st.success("✅ Sipariş öncesi in-silico kontroller geçti. Yine de ıslak laboratuvar doğrulaması gerekir.")
    elif gate["overall"] == "WARN":
        st.warning("⚠️ Sipariş öncesi uyarılar var. Aşağıdaki maddeleri inceleyin.")
    else:
        st.error("❌ Sipariş listesi kilitli: kritik kontrol başarısız.")

    gate_rows = [{
        "Kontrol": c["label"],
        "Durum": _status_text(c["status"]),
        "Ayrıntı": c["detail"],
    } for c in gate["checks"]]
    st.dataframe(pd.DataFrame(gate_rows), width="stretch", hide_index=True)

    order_rows = order_sheet_rows(
        active, reporter=reporter, quencher=quencher, scale=scale, purification=purification
    )
    st.markdown("#### Sipariş listesi önizlemesi")
    st.dataframe(pd.DataFrame(order_rows), width="stretch", hide_index=True)

    allow_warn = True
    if gate["overall"] == "WARN":
        allow_warn = st.checkbox(
            "WARN uyarılarını inceledim; bu risklerle devam ediyorum.", key="order_ack"
        )
    can_download = gate["overall"] == "PASS" or (gate["overall"] == "WARN" and allow_warn)

    if can_download:
        d1, d2 = st.columns(2)
        d1.download_button(
            "📥 Sipariş listesini indir",
            data=pd.DataFrame(order_rows).to_csv(index=False).encode("utf-8-sig"),
            file_name=f"siparis_listesi_{active.get('accession', 'assay')}.csv",
            mime="text/csv",
            key="order_sheet_download",
        )
        d2.download_button(
            "📥 QC karar raporunu indir",
            data=pd.DataFrame(gate_rows).to_csv(index=False).encode("utf-8-sig"),
            file_name=f"siparis_qc_{active.get('accession', 'assay')}.csv",
            mime="text/csv",
            key="order_qc_download",
        )
    else:
        st.info("Sipariş listesi, FAIL/BLOCKED kontrolleri giderilmeden indirilemez.")

    if gate.get("blocking"):
        with st.expander("Siparişi bloke eden maddeler"):
            for item in gate["blocking"]:
                st.write("• " + item)


# ---------------------------------------------------------------------------
# KÜTÜPHANE
# ---------------------------------------------------------------------------

def render_library_tabs(ctx):
    st.markdown(
        "Yayınlanmış assay'leri yapıştırıp referans üzerinde konumlandırın; "
        "aktif assay ile karşılaştırın ve JSON kütüphanesinde saklayın. "
        "**Kütüphaneye doğrulanmamış hazır sekans gömülmez.**"
    )
    sub_lit, sub_track = st.tabs(["📖 Literatür / Assay Karşılaştırma", "🗺️ Genome Track"])
    with sub_lit:
        _render_literature(ctx)
    with sub_track:
        _render_track(ctx)


def _render_literature(ctx):
    active = ctx["active_or_demo"]("literature")
    if not active:
        return
    _active_header(active)
    st.markdown("#### Yayınlanmış assay'ı konumlandır ve karşılaştır")
    c1, c2 = st.columns(2)
    name = c1.text_input("Assay adı", key="lit_name", placeholder="ör. WHO CCHFV RT-qPCR")
    source = c2.text_input("Kaynak / yayın", key="lit_source", placeholder="ör. WHO 2022")
    f = st.text_input("Forward primer (5'→3')", key="lit_f")
    r = st.text_input("Reverse primer (5'→3')", key="lit_r")
    p = st.text_input("Prob (5'→3', opsiyonel)", key="lit_p")
    notes = st.text_area("Not", key="lit_notes", height=70)
    ref_choice = st.radio(
        "Referans",
        ["Aktif assay referansı", "GenBank ID"],
        horizontal=True, key="lit_ref_choice",
    )
    external_accession = ""
    if ref_choice == "GenBank ID":
        external_accession = st.text_input("GenBank ID", key="lit_accession",
                                           placeholder="ör. NC_005302.1")

    if st.button("📍 Konumlandır ve aktif assay ile karşılaştır", type="primary",
                 key="lit_run", width="stretch"):
        if not (f and r):
            st.error("Forward ve reverse primer zorunlu.")
        else:
            ref_seq, ref_offset, ref_label = active["ref_seq"], active.get("ref_offset", 0), "aktif assay referansı"
            if ref_choice == "GenBank ID" and external_accession.strip():
                fetched = ctx["fetch_genome"](external_accession.strip())
                if fetched is None:
                    st.error("GenBank dizisi alınamadı.")
                    return
                ref_seq, _org = fetched
                ref_offset = 0
                ref_label = external_accession.strip()
                if len(ref_seq) > 200_000:
                    st.warning("Büyük genom: yalnızca aktif assay'ın hedef bölgesi taranacak.")
                    ref_seq = active["ref_seq"]
                    ref_offset = active.get("ref_offset", 0)
            with st.spinner("Assay konumlandırılıyor..."):
                loc = locate_assay(ref_seq, f, r, p or None, max_mismatches=4,
                                   min_amp=30, max_amp=2000)
            if not loc.get("found"):
                st.error(loc.get("message", "Assay konumlandırılamadı."))
            else:
                best = loc["best"]
                comparable = ref_offset == active.get("ref_offset", 0) and ref_seq == active["ref_seq"]
                cmp = None
                if comparable:
                    active_cmp = {
                        "forward": active["forward"], "reverse": active["reverse"],
                        "probe": active.get("probe", ""),
                        "amplicon_start": active["amplicon_start"],
                        "amplicon_end": active["amplicon_end"],
                        "amplicon": active["amplicon"],
                    }
                    cmp = compare_assays(active_cmp, best)
                st.session_state["lit_last"] = {
                    "loc": loc, "best": best, "compare": cmp,
                    "ref_label": ref_label, "ref_offset": ref_offset,
                    "name": name, "source": source, "f": f, "r": r, "p": p, "notes": notes,
                }

    lit = st.session_state.get("lit_last")
    if lit:
        best = lit["best"]
        st.markdown("#### Konum")
        st.dataframe(pd.DataFrame([{
            "Referans": lit["ref_label"],
            "Amplikon başlangıç": best["amplicon_start_1based"] + lit["ref_offset"],
            "Amplikon bitiş": best["amplicon_end_1based"] + lit["ref_offset"],
            "Boyut (bp)": best["size"],
            "F mismatch": best["forward"]["n_mismatches"],
            "R mismatch": best["reverse"]["n_mismatches"],
            "Prob mismatch": (best.get("probe") or {}).get("n_mismatches"),
            "Toplam mismatch": best["total_mismatches"],
        }]), width="stretch", hide_index=True)
        with st.expander("Yayın assay'ının dizileri ve Tm"):
            st.code(
                f"Forward 5'-{best.get('forward_seq', '')}-3'  Tm {best.get('forward_tm')} °C\n"
                f"Reverse 5'-{best.get('reverse_seq', '')}-3'  Tm {best.get('reverse_tm')} °C\n"
                f"Prob    5'-{best.get('probe_seq') or ''}-3'  Tm {best.get('probe_tm')} °C",
                language="text",
            )
        if lit["compare"] and lit["compare"].get("ok"):
            cmp = lit["compare"]
            st.markdown("#### Aktif assay ile uyum")
            st.metric("Uyum skoru", f"{cmp['concordance_score']}/100")
            st.dataframe(comparison_table(cmp), width="stretch", hide_index=True)
            st.caption(
                "Uyum skoru: %30 amplikon örtüşmesi + %30 amplikon identitesi + "
                "%40 oligo identitesi. Farklı hedef bölge ise ceza uygulanır."
            )
        elif lit["compare"]:
            st.warning(lit["compare"].get("message", "Karşılaştırma yapılamadı."))
        else:
            st.info("Farklı referans seçildiği için aktif assay ile doğrudan karşılaştırma yapılmadı.")

        if st.button("➕ Bu assay'ı kütüphaneye ekle", key="lit_add"):
            entry = library_entry(
                lit.get("name") or f"Assay {len(st.session_state['assay_library']) + 1}",
                lit["ref_label"], lit["f"], lit["r"], lit["p"],
                notes=lit.get("notes", ""), source=lit.get("source", ""),
            )
            ok, msg = validate_entry(entry)
            existing = {(e.get("name"), e.get("forward"), e.get("reverse"))
                        for e in st.session_state["assay_library"]}
            if not ok:
                st.error(msg)
            elif (entry["name"], entry["forward"], entry["reverse"]) in existing:
                st.info("Bu assay zaten kütüphanede.")
            else:
                st.session_state["assay_library"].append(entry)
                st.success("Kütüphaneye eklendi.")
                st.rerun()

    st.markdown("#### 📚 Assay kütüphanesi")
    entries = st.session_state.get("assay_library", [])
    if entries:
        st.dataframe(library_dataframe(entries), width="stretch", hide_index=True)
    else:
        st.info("Kütüphane boş. Konumlandırılmış bir assay ekleyin veya aktif assay'ı ekleyin.")
    b1, b2, b3 = st.columns(3)
    if b1.button("➕ Aktif assay'ı kütüphaneye ekle", key="lib_add_active"):
        entry = library_entry(
            active.get("virus", "Aktif assay"), active.get("accession", ""),
            active["forward"], active["reverse"], active.get("probe", ""),
            notes="Uygulama tasarımı", source="PCR Primer & Probe Designer v3.0",
        )
        existing = {(e.get("name"), e.get("forward"), e.get("reverse"))
                    for e in st.session_state["assay_library"]}
        if (entry["name"], entry["forward"], entry["reverse"]) in existing:
            st.info("Bu assay zaten kütüphanede.")
        else:
            st.session_state["assay_library"].append(entry)
            st.success("Aktif assay kütüphaneye eklendi.")
            st.rerun()
    if entries:
        b2.download_button(
            "⬇️ Kütüphaneyi indir (JSON)",
            data=library_to_json(entries),
            file_name="assay_kutuphanesi.json",
            mime="application/json",
            key="lib_json",
        )
        if b3.button("🗑️ Kütüphaneyi temizle", key="lib_clear"):
            st.session_state["assay_library"] = []
            st.rerun()
    uploaded = st.file_uploader("JSON kütüphane yükle", type=["json"], key="lib_upload")
    if uploaded is not None:
        try:
            text = uploaded.read().decode("utf-8")
            loaded, err = library_from_json(text)
            if err:
                st.error(err)
            elif loaded:
                existing = {(e.get("name"), e.get("forward"), e.get("reverse"))
                            for e in st.session_state["assay_library"]}
                added = 0
                for e in loaded:
                    key = (e.get("name"), e.get("forward"), e.get("reverse"))
                    if key not in existing:
                        st.session_state["assay_library"].append(e)
                        existing.add(key)
                        added += 1
                if added:
                    st.success(f"{added} yeni kayıt yüklendi.")
                    st.rerun()
                else:
                    st.info("Dosyadaki kayıtlar kütüphanede zaten var.")
        except Exception as exc:
            st.error(f"Dosya okunamadı: {exc}")


def _render_track(ctx):
    active = ctx["active_or_demo"]("track")
    if not active:
        return
    _active_header(active)
    st.markdown("#### Genom bağlamı")
    flank = st.slider("Gösterilecek kanat (bp)", 50, 1000, 200, 50, key="track_flank")
    reg = ctx["active_design_region"](active, flank)
    local_start, local_end = reg["local_start"], reg["local_end"]
    genomic_start, genomic_end = reg["genomic_start"], reg["genomic_end"]

    track_len = local_end - local_start
    features = []
    info = active.get("info", {})
    if info.get("gene_start") and info.get("gene_end"):
        features.append({
            "name": info.get("target_gene") or "Hedef gen",
            "start": int(info["gene_start"]) - 1,
            "end": int(info["gene_end"]),
            "color": "#90caf9",
        })
    else:
        features.append({
            "name": "Hedef bölge",
            "start": active.get("start", local_start) - 1,
            "end": active.get("end", local_end),
            "color": "#b0bec5",
        })

    primers = []
    for name, s, e, color in (
        ("F", active.get("forward_start"), active.get("forward_end"), "#1565c0"),
        ("P", active.get("probe_start"), active.get("probe_end"), "#f9a825"),
        ("R", active.get("reverse_start"), active.get("reverse_end"), "#6a1b9a"),
    ):
        if s and e:
            primers.append({"name": name, "start": int(s) - 1, "end": int(e), "color": color})

    cons = None
    ci = active.get("cons_info")
    if ci and ci.get("conservation"):
        arr = ci["conservation"]
        # Korunmuşluk dizisi aktif hedef bölgenin (region_start_local) başından başlar.
        base_genomic = active.get("ref_offset", 0) + active.get("region_start_local", 0) + 1
        cons = []
        for i in range(track_len):
            genomic = active.get("ref_offset", 0) + local_start + i + 1
            idx = genomic - base_genomic
            cons.append(arr[idx] if 0 <= idx < len(arr) else None)
        if not any(v is not None for v in cons):
            cons = None

    density = None
    geno = st.session_state.get("last_genotype")
    if geno and geno.get("matrix"):
        d = ctx["genotype_mismatch_density"](geno)
        amp_lo = geno["amplicon_start"]
        amp_hi = geno["amplicon_end"]
        if d and amp_hi > amp_lo:
            # d[0] amplikonun ilk bazına (amp_lo) karşılık gelir.
            skip = max(0, local_start - amp_lo)
            take = max(0, min(len(d) - skip, track_len))
            density = [0] * max(0, amp_lo - local_start) + d[skip:skip + take]
            density = (density + [0] * track_len)[:track_len]

    fig = build_genome_track(
        genomic_start - 1, genomic_end,
        features=features, primers=primers, conservation=cons,
        mismatch_density=density,
        title=f"{active.get('virus')} · {genomic_start}–{genomic_end}",
        ref_len=active.get("seq_len"),
    )
    st.plotly_chart(fig, width="stretch")
    st.caption(
        "Korunmuşluk ve mismatch izleri yalnızca mevcut oturumdaki analizlerden gelir. "
        "Büyük genomlarda gösterilen bölge aktif assay çevresidir."
    )


# ---------------------------------------------------------------------------
# İLERİ
# ---------------------------------------------------------------------------

def render_advanced_tabs(ctx):
    st.markdown(
        "Consensus/dejenere primerler, SNP ayrımı, RPA/LAMP önerileri ve panel optimizasyonu. "
        "RPA/LAMP **beta** modüllerdir; ıslak laboratuvar optimizasyonu gerekir."
    )
    sub_cons, sub_snp, sub_iso, sub_panel = st.tabs(
        ["🧬 Consensus / Dejenere", "🔬 SNP Ayrımı", "🧪 LAMP / RPA", "🧩 Panel Optimizer"]
    )
    with sub_cons:
        _render_consensus(ctx)
    with sub_snp:
        _render_snp(ctx)
    with sub_iso:
        _render_isothermal(ctx)
    with sub_panel:
        _render_panel_opt(ctx)


def _render_consensus(ctx):
    active = ctx["active_or_demo"]("consensus")
    if not active:
        return
    _active_header(active)
    reg = ctx["active_design_region"](active, 100)
    st.caption(
        f"Tasarım bölgesi: {reg['genomic_start']}–{reg['genomic_end']} "
        f"({len(reg['region'])} bp). NCBI'dan homolog diziler çekilir."
    )
    c1, c2, c3 = st.columns(3)
    n_seq = c1.slider("Dizi sayısı", 3, 20, 8, 1, key="cons_n")
    min_cov = c2.slider("Minimum kapsama (%)", 50, 100, 80, 5, key="cons_cov") / 100.0
    max_deg = c3.selectbox("Maksimum dejenere kat sayısı", [64, 256, 1024], index=1,
                           key="cons_deg")
    if st.button("🧬 Consensus ve dejenere primerleri tasarla", type="primary",
                 key="cons_run", width="stretch"):
        organism = active.get("organism") or ""
        if not organism:
            st.error("Organizma adı yok; NCBI modunda tasarım yapın.")
            return
        big = (active.get("seq_len", 0) > 200_000) or (not active.get("ref_is_full"))
        keywords = ()
        if active.get("info", {}).get("target_gene"):
            keywords = (active["info"]["target_gene"],)
        with st.spinner("Homolog diziler çekiliyor..."):
            subjects = ctx["get_genotype_subjects"](
                organism, active.get("accession", ""), reg["region"],
                n_seq, keywords, big,
            )
        if not subjects:
            st.warning("Hizalanabilir dizi bulunamadı.")
            return
        with st.spinner("Consensus ve dejenere primerler hesaplanıyor..."):
            result = design_consensus_primers(
                reg["region"], subjects,
                region_start=active.get("ref_offset", 0) + reg["local_start"],
                target_amp=ctx.get("target_amplicon", 150), top_n=10,
                min_coverage=min_cov, max_degeneracy=int(max_deg),
            )
            # Kullanıcıyı boş ekranla bırakmamak için bir kez gevşetilmiş deneme.
            if result.get("ok") and not result.get("candidates"):
                relaxed = design_consensus_primers(
                    reg["region"], subjects,
                    region_start=active.get("ref_offset", 0) + reg["local_start"],
                    target_amp=ctx.get("target_amplicon", 150), top_n=10,
                    min_coverage=max(0.50, min_cov - 0.25),
                    max_degeneracy=min(1024, max(4, int(max_deg) * 4)),
                )
                if relaxed.get("candidates"):
                    relaxed["relaxed"] = True
                    result = relaxed
        st.session_state["consensus_last"] = result

    result = st.session_state.get("consensus_last")
    if not result:
        st.info("Dizileri çekip consensus tasarımını çalıştırın.")
        return
    if not result.get("ok"):
        st.error(result.get("message", "Consensus üretilemedi."))
        return
    if result.get("relaxed"):
        st.warning(
            "⚠️ Seçilen kapsama/dejenerasyon eşiğinde aday yoktu; eşikler otomatik "
            "gevşetilerek gösterilen adaylar üretildi. Laboratuvarda doğrulayın."
        )
    cons = result["consensus"]
    m = st.columns(4)
    m[0].metric("Hizalanan dizi", cons.get("n"))
    m[1].metric("Ortalama identite", cons.get("mean_identity") or "—")
    m[2].metric("Consensus uzunluğu", len(result["majority"]))
    m[3].metric("Dejenere pozisyon", sum(1 for b in result["degenerate"] if b not in "ACGT"))
    with st.expander("Consensus dizileri"):
        st.markdown("**Majority (ACGT)**")
        st.code(result["majority"], language="text")
        st.markdown("**Dejenere (IUPAC)**")
        st.code(result["degenerate"], language="text")
        fasta = f">consensus_majority\n{result['majority']}\n>consensus_degenerate\n{result['degenerate']}\n"
        st.download_button(
            "📥 Consensus FASTA indir",
            data=fasta.encode("utf-8"),
            file_name=f"consensus_{active.get('accession', 'assay')}.fasta",
            mime="text/plain",
            key="cons_fasta",
        )
    candidates = result.get("candidates") or []
    if not candidates:
        st.warning("Kapsama/dejenerasyon kriterlerini sağlayan aday yok. Eşikleri gevşetin.")
        return
    rows = []
    for c in candidates:
        rows.append({
            "Forward": c["forward"],
            "F dejenere": c["forward_degenerate"],
            "F kat": c["forward_deg"],
            "Reverse": c["reverse"],
            "R dejenere": c["reverse_degenerate"],
            "R kat": c["reverse_deg"],
            "Prob": c["probe"],
            "P dejenere": c["probe_degenerate"],
            "P kat": c["probe_deg"],
            "F kapsama %": c["coverage_f"],
            "R kapsama %": c["coverage_r"],
            "P kapsama %": c["coverage_p"],
            "F Tm": c["forward_tm"],
            "R Tm": c["reverse_tm"],
            "P Tm": c["probe_tm"],
            "Amplikon": c["amplicon_size"],
            "Başlangıç": c["forward_start"],
            "Bitiş": c["reverse_end"],
        })
    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
    st.caption(
        "Reverse dejenere dizi primer yönündedir (5'→3'). Kapsama, referans bölgeye hizalanan "
        "dizilerde ≤1 mismatch oranıdır."
    )


def _render_snp(ctx):
    active = ctx["active_or_demo"]("snp")
    if not active:
        return
    _active_header(active)
    reg = ctx["active_design_region"](active, 100)
    st.caption(f"Tasarım bölgesi: {reg['genomic_start']}–{reg['genomic_end']} (genomik).")
    default_snp = int(active.get("probe_start") or ((reg["genomic_start"] + reg["genomic_end"]) // 2))
    c1, c2, c3 = st.columns(3)
    snp_pos = c1.number_input("SNP pozisyonu (1-based genomik)", min_value=1,
                              max_value=max(1, active.get("seq_len", 1)),
                              value=max(1, default_snp), step=1, key="snp_pos")
    alt = c2.selectbox("Alternatif baz", ["A", "C", "G", "T"], key="snp_alt")
    target_amp = c3.slider("Amplikon hedefi", 70, 250, 120, 10, key="snp_amp")
    reverse = st.text_input("Reverse primer (5'→3', opsiyonel)", key="snp_reverse")
    local_snp = int(snp_pos) - 1 - active.get("ref_offset", 0)
    ref_base = ""
    if 0 <= local_snp < len(active["ref_seq"]):
        ref_base = active["ref_seq"][local_snp]
    st.caption(f"Referans baz: **{ref_base or '—'}** · Alternatif: **{alt}**")
    if st.button("🔬 SNP ayrım primerlerini tasarla", type="primary",
                 key="snp_run", width="stretch"):
        if not ref_base:
            st.error("SNP pozisyonu aktif referans aralığında değil.")
            return
        if alt == ref_base:
            st.error("Alternatif baz referans baz ile aynı.")
            return
        result = design_snp_primers(
            active["ref_seq"], reg["local_start"], reg["local_end"],
            local_snp, alt, reverse=reverse.strip() or None,
            target_amp=target_amp, top_n=12,
        )
        st.session_state["snp_last"] = result

    result = st.session_state.get("snp_last")
    if not result:
        st.info("SNP pozisyonunu ve alternatif bazı seçip tasarımı çalıştırın.")
        return
    if not result.get("ok"):
        st.error(result.get("message", "SNP tasarımı üretilemedi."))
        return
    candidates = result.get("candidates") or []
    if not candidates:
        st.warning("Uygun SNP primer çifti bulunamadı; bölgeyi/amplikon boyunu değiştirin.")
        return
    offset = active.get("ref_offset", 0)
    rows = []
    for c in candidates:
        rows.append({
            "Forward (mutant)": c["forward"],
            "Forward (WT kontrol)": c["forward_wt"],
            "Reverse": c["reverse"],
            "F Tm": c["forward_tm"],
            "R Tm": c["reverse_tm"],
            "F GC": c["forward_gc"],
            "R GC": c["reverse_gc"],
            "3' SNP": c["alt_base"],
            "WT son5 mismatch": c["wt_last5_mismatch"],
            "Mut son5 mismatch": c["alt_last5_mismatch"],
            "3' mismatch (WT)": c["terminal_mismatch_wt"],
            "Ayrım": c["discrimination"],
            "Amplikon": c["amplicon_size"],
            "SNP (genomik)": c["snp_pos"] + offset,
            "F başlangıç (genomik)": c["forward_start"] + offset,
            "R bitiş (genomik)": c["reverse_end"] + offset,
        })
    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
    st.caption(
        "Mutant forward primerin 3' ucu SNP bazıdır. WT kontrol primerı aynı bölgeyi "
        "referans bazla temsil eder. Ayrım gücü, 3' uç ve son 5 baz mismatch sayısına göre "
        "sınıflandırılır; gerçek ayrımı qPCR standardı ve polimeraz hassasiyeti belirler."
    )
    with st.expander("Primer3 heterodimer bilgilendirici ΔG"):
        st.dataframe(pd.DataFrame([{
            "Forward": c["forward"], "ΔG farkı": c["delta_dg"],
            "Mut vs WT Tm farkı": c["tm_delta_wt"],
        } for c in candidates]), width="stretch", hide_index=True)
        st.caption("Bu ΔG, primer3'ün en iyi lokal hizalamasıdır; 3' ucunu zorlamaz. Ana metrik değildir.")
    st.download_button(
        "📥 SNP adaylarını indir (CSV)",
        data=pd.DataFrame(rows).to_csv(index=False).encode("utf-8"),
        file_name=f"snp_adaylari_{active.get('accession', 'assay')}.csv",
        mime="text/csv",
        key="snp_csv",
    )


def _render_isothermal(ctx):
    st.warning(
        "⚠️ RPA ve LAMP tasarımları **beta**dır. Primer/prob dizileri laboratuvarda "
        "optimize edilmeden kullanılmamalıdır; LAMP için loop primerleri ve reaksiyon "
        "sıcaklığı ayrıca doğrulanmalıdır."
    )
    active = ctx["active_or_demo"]("isothermal")
    reg = ctx["active_design_region"](active, 100) if active else None
    sub_rpa, sub_lamp, sub_check = st.tabs(["RPA", "LAMP", "LAMP Kontrol"])
    with sub_rpa:
        st.markdown("#### RPA primer + prob önerisi")
        if not active:
            st.info("RPA için aktif assay gerekir.")
        else:
            st.caption(f"Bölge: {reg['genomic_start']}–{reg['genomic_end']} ({len(reg['region'])} bp)")
            amp = st.slider("RPA amplikon hedefi", 120, 300, 200, 10, key="rpa_amp")
            if st.button("🧪 RPA adaylarını üret", type="primary", key="rpa_run",
                         width="stretch"):
                result = design_rpa(
                    reg["region"],
                    region_start=active.get("ref_offset", 0) + reg["local_start"],
                    target_amp=amp, top_n=6,
                )
                st.session_state["rpa_last"] = result
            result = st.session_state.get("rpa_last")
            if result:
                if not result.get("ok"):
                    st.error(result.get("message", "RPA tasarımı üretilemedi."))
                elif not result.get("candidates"):
                    st.warning("RPA adayı bulunamadı; bölgeyi genişletin.")
                else:
                    st.dataframe(pd.DataFrame([{
                        "Forward": c["forward"], "Reverse": c["reverse"],
                        "F Tm": c["forward_tm"], "R Tm": c["reverse_tm"],
                        "F GC": c["forward_gc"], "R GC": c["reverse_gc"],
                        "Amplikon": c["amplicon_size"],
                        "Başlangıç": c["forward_start"], "Bitiş": c["reverse_end"],
                    } for c in result["candidates"]]), width="stretch", hide_index=True)
                    for i, c in enumerate(result["candidates"], start=1):
                        with st.expander(f"Aday {i} · {c['amplicon_size']} bp"):
                            st.code(
                                f"Forward 5'-{c['forward']}-3'\n"
                                f"Reverse 5'-{c['reverse']}-3'",
                                language="text",
                            )
                            pr = c.get("probe")
                            if pr:
                                st.markdown(
                                    f"**Prob**: 5'-[{pr['label_5']}]-{pr['sequence']}-"
                                    f"[{pr['quencher_3']}]-3'  \n"
                                    f"THF pozisyonu: {pr['thf_position']} · "
                                    f"Tm {pr['tm']} °C · GC {pr['gc']}%"
                                )
                            st.caption("Standart NN Tm RPA için bilgilendiricidir; 37–42 °C reaksiyon optimize edilmelidir.")

    with sub_lamp:
        st.markdown("#### LAMP F3/B3/FIP/BIP önerisi")
        if not active:
            st.info("LAMP için aktif assay gerekir.")
        else:
            st.caption(f"Bölge: {reg['genomic_start']}–{reg['genomic_end']} ({len(reg['region'])} bp)")
            win = st.slider("LAMP hedef pencere", 160, 300, 220, 10, key="lamp_win")
            if st.button("🧪 LAMP setlerini üret", type="primary", key="lamp_run",
                         width="stretch"):
                with st.spinner("LAMP adayları taranıyor..."):
                    result = design_lamp(
                        reg["region"],
                        region_start=active.get("ref_offset", 0) + reg["local_start"],
                        target_amp=win, top_n=5,
                    )
                st.session_state["lamp_last"] = result
            result = st.session_state.get("lamp_last")
            if result:
                if not result.get("ok"):
                    st.error(result.get("message", "LAMP tasarımı üretilemedi."))
                elif not result.get("candidates"):
                    st.warning("LAMP seti bulunamadı; bölgeyi/pencereyi değiştirin.")
                else:
                    st.dataframe(pd.DataFrame([{
                        "F3": c["f3"], "B3": c["b3"], "FIP": c["fip"], "BIP": c["bip"],
                        "F3 Tm": c["f3_tm"], "B3 Tm": c["b3_tm"],
                        "FIP Tm": c["fip_tm"], "BIP Tm": c["bip_tm"],
                        "Pencere": f"{c['window_start']}–{c['window_end']}",
                        "Skor": c["score"],
                    } for c in result["candidates"]]), width="stretch", hide_index=True)
                    for i, c in enumerate(result["candidates"], start=1):
                        with st.expander(f"LAMP seti {i} · {c['window_start']}–{c['window_end']}"):
                            st.code(
                                f"F3  5'-{c['f3']}-3'\n"
                                f"B3  5'-{c['b3']}-3'\n"
                                f"FIP 5'-{c['fip']}-3'\n"
                                f"BIP 5'-{c['bip']}-3'",
                                language="text",
                            )
                            st.caption(
                                f"FIP = F1c + F2; BIP = B1c + B2. "
                                f"F3 Tm {c['f3_tm']} °C, B3 Tm {c['b3_tm']} °C, "
                                f"FIP Tm {c['fip_tm']} °C, BIP Tm {c['bip_tm']} °C."
                            )
                    st.caption("Loop primerleri (LF/LB) bu beta tasarıma dahil değildir.")

    with sub_check:
        st.markdown("#### Yayınlanmış LAMP setini kontrol et")
        st.caption("Dizileri 5'→3' yönünde girin. Dimer ve 3' uç özellikleri raporlanır.")
        f3 = st.text_input("F3", key="chk_f3")
        b3 = st.text_input("B3", key="chk_b3")
        fip = st.text_input("FIP", key="chk_fip")
        bip = st.text_input("BIP", key="chk_bip")
        lf = st.text_input("LF (opsiyonel)", key="chk_lf")
        lb = st.text_input("LB (opsiyonel)", key="chk_lb")
        if st.button("Kontrol et", key="chk_run"):
            if not all([f3, b3, fip, bip]):
                st.error("F3, B3, FIP ve BIP zorunlu.")
            else:
                out = check_lamp_set(f3, b3, fip, bip, lf, lb)
                st.dataframe(out["table"], width="stretch", hide_index=True)
                if not out["conflicts"].empty:
                    st.warning("Dimer/komplementerlik uyarıları:")
                    st.dataframe(out["conflicts"], width="stretch", hide_index=True)
                else:
                    st.success("Belirgin dimer uyarısı yok.")


def _render_panel_opt(ctx):
    st.markdown("#### Panel uyumluluğu ve kanal ataması")
    active = st.session_state.get("active_assay")
    panel = list(st.session_state.get("panel_assays") or [])
    if not panel:
        st.info(
            "Panel assay'leri **Multiplex Panel** sekmesinde tasarlandığında otomatik olarak "
            "buraya aktarılır. Hemen denemek için aktif assay'ı ekleyebilirsiniz."
        )
        if active:
            if st.button("➕ Aktif assay'ı panele ekle", key="panel_add_active"):
                panel.append(ctx["active_to_panel_assay"](active))
                st.session_state["panel_assays"] = panel
                st.rerun()
        else:
            st.caption("Aktif assay yok. Önce Tek Etken Tasarımı veya örnek assay yükleyin.")
        return

    summary = panel_summary(panel)
    m = st.columns(4)
    m[0].metric("Assay sayısı", summary.get("n"))
    m[1].metric("Amplikon aralığı", f"{summary.get('amplicon_min')}–{summary.get('amplicon_max')} bp")
    m[2].metric("Prob Tm aralığı", f"{summary.get('probe_tm_min')}–{summary.get('probe_tm_max')} °C")
    m[3].metric("Kullanılan kanal", summary.get("channels"))
    for w in summary.get("warnings", []):
        st.warning(w)

    st.markdown("#### Kanal ataması")
    table, warnings = channel_assignment(panel)
    st.dataframe(table, width="stretch", hide_index=True)
    st.download_button(
        "📥 Kanal tablosunu indir (CSV)",
        data=table.to_csv(index=False).encode("utf-8"),
        file_name="panel_kanal_atamasi.csv",
        mime="text/csv",
        key="panel_csv",
    )

    st.markdown("#### Assay'ler arası dimer riski")
    dimers = panel_dimer_table(panel)
    if dimers.empty:
        st.success("✅ ΔG ≤ -6 kcal/mol olan assay'ler arası dimer yok.")
    else:
        st.dataframe(dimers, width="stretch", hide_index=True)
        st.warning("Kırmızı eşik: ΔG ≤ -6 kcal/mol. Sentez öncesi dimer analizi yapın.")

    st.markdown("#### Uyumluluk matrisi")
    mat = compatibility_matrix(panel)
    try:
        import plotly.graph_objects as go
        heat = go.Figure(go.Heatmap(
            z=mat.values.tolist(), x=list(mat.columns), y=list(mat.index),
            zmin=-12, zmax=0, colorscale=[[0, "#b71c1c"], [0.5, "#ffffff"], [1, "#e3f2fd"]],
            text=[[f"{v:.2f}" for v in row] for row in mat.values.tolist()],
            texttemplate="%{text}", hovertemplate="%{y}–%{x}<br>ΔG %{z:.2f}<extra></extra>",
        ))
        heat.update_layout(height=max(320, 32 * len(panel)), margin=dict(l=20, r=20, t=30, b=20))
        st.plotly_chart(heat, width="stretch")
    except Exception:
        st.dataframe(mat, width="stretch")

    c1, c2 = st.columns(2)
    if active:
        names = {a.get("name") for a in panel}
        if active.get("virus") not in names and c1.button("➕ Aktif assay'ı panele ekle", key="panel_add_active2"):
            panel.append(ctx["active_to_panel_assay"](active))
            st.session_state["panel_assays"] = panel
            st.rerun()
    if c2.button("🗑️ Paneli temizle", key="panel_clear"):
        st.session_state["panel_assays"] = []
        st.rerun()
