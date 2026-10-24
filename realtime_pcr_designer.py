import os
import re
from datetime import datetime

import pandas as pd
import streamlit as st
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pcr_engine import (
    calculate_gc_content,
    calculate_tm,
    design_primers_best,
    design_probe,
    revcomp,
    validate_primer_primer,
)
from conservation import (
    best_conserved_windows,
    compute_conservation,
    fetch_related_sequences,
    region_keywords,
    window_means,
)
from report import build_panel_workbook

# --- 1. AYARLAR ---
st.set_page_config(
    page_title="Real-Time PCR Primer & Probe Designer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# NCBI nezaket kuralı: kendi e-postanızı yazın, NCBI_EMAIL ortam değişkeni
# veya Streamlit Cloud'da st.secrets["NCBI_EMAIL"] tanımlayın.
def _ncbi_email() -> str:
    try:
        if "NCBI_EMAIL" in st.secrets:
            return str(st.secrets["NCBI_EMAIL"])
    except Exception:
        pass
    return os.environ.get("NCBI_EMAIL", "hakan.hedef@example.com")


Entrez.email = _ncbi_email()

DEFAULT_VIRUS_DATABASE = {
    "Kırım Kongo Kanamalı Ateşi (CCHFV)": {
        "segment": "S Segmenti (Nükleokapsid Proteini)",
        "genbank_id": "NC_005302.1",
        "target_region": "1100 - 1250 bp",
        "reference": "WHO & Corman Arbovirus Panel",
        "database": "NCBI",
    },
    "Batı Nil Virüsü (WNV)": {
        "segment": "Envelope (E) Geni",
        "genbank_id": "NC_001563.2",
        "target_region": "1150 - 1280 bp",
        "reference": "CDC Diagnostic Panel",
        "database": "NCBI",
    },
    "Hantavirüs (Dobrava/Puumala)": {
        "segment": "S Segmenti (Nükleokapsid)",
        "genbank_id": "NC_005228.1",
        "target_region": "450 - 580 bp",
        "reference": "European RefSeq Validation",
        "database": "NCBI",
    },
    "Dengue Virüsü (DENV 1-4)": {
        "segment": "3' UTR / Capsid Geni",
        "genbank_id": "NC_001474.2",
        "target_region": "150 - 270 bp",
        "reference": "Lancet Infectious Diseases",
        "database": "NCBI",
    },
    "Chikungunya Virüsü (CHIKV)": {
        "segment": "nsP1 Geni",
        "genbank_id": "NC_004162.2",
        "target_region": "850 - 970 bp",
        "reference": "PAHO/WHO Arbovirus",
        "database": "NCBI",
    },
    "Zika Virüsü (ZIKV)": {
        "segment": "Membrane/Envelope Bölgesi",
        "genbank_id": "NC_012532.1",
        "target_region": "1020 - 1140 bp",
        "reference": "CDC Real-Time RT-PCR",
        "database": "NCBI",
    },
    "Ebola Virüsü (EBOV)": {
        "segment": "NP (Nucleoprotein) Geni",
        "genbank_id": "NC_002549.1",
        "target_region": "800 - 920 bp",
        "reference": "USAMRIID / WHO",
        "database": "NCBI",
    },
    "Marburg Virüsü (MARV)": {
        "segment": "GP (Glycoprotein) Geni",
        "genbank_id": "NC_001608.3",
        "target_region": "1200 - 1320 bp",
        "reference": "WHO Filovirus",
        "database": "NCBI",
    },
}

if "virus_db" not in st.session_state:
    st.session_state.virus_db = dict(DEFAULT_VIRUS_DATABASE)
VIRUS_DB = st.session_state.virus_db


# --- 2. NCBI FONKSİYONLARI (önbellekli) ---
@st.cache_data(ttl=3600, show_spinner=False)
def fetch_record(accession: str):
    """GenBank kaydını çeker. Hata durumunda None döner."""
    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            return SeqIO.read(handle, "genbank")
    except Exception:
        return None


@st.cache_data(ttl=3600, show_spinner=False)
def get_conservation(accession, seq_str, start, end, organism, n_seqs, keywords=()):
    """İlgili dizileri çekip pozisyon bazlı korunmuşluk hesaplar."""
    subjects = fetch_related_sequences(
        organism,
        exclude_accessions=[accession],
        max_n=int(n_seqs),
        reference_seq=seq_str,
        keywords=tuple(keywords),
    )
    info = compute_conservation(seq_str, start, end, subjects)
    info["requested"] = int(n_seqs)
    return info


@st.cache_data(ttl=600, show_spinner=False)
def check_ncbi_connection():
    try:
        handle = Entrez.esearch(db="nucleotide", term="Arbovirus", retmax=1)
        found = Entrez.read(handle)
        handle.close()
        return True, f"{len(found.get('IdList', []))} kayıt önizlemesi"
    except Exception as exc:
        return False, str(exc)[:150]


# --- 3. YARDIMCILAR ---
def parse_region(region_text, seq_len):
    """'1100 - 1250 bp' -> (1099, 1250) 0-based yarı açık."""
    numbers = re.findall(r"\d+", str(region_text))
    if len(numbers) >= 2:
        start = max(0, int(numbers[0]) - 1)
        end = min(int(numbers[1]), seq_len)
    else:
        start, end = 0, min(1500, seq_len)
    if end - start < 60:
        start, end = 0, min(1500, seq_len)
    return start, end


def clean_fasta(text):
    """FASTA veya ham diziden sadece ACGTN bırakır."""
    if not text:
        return ""
    cleaned = "".join(
        line.strip() for line in text.splitlines() if not line.strip().startswith(">")
    ).upper()
    return re.sub(r"[^ACGTN]", "", cleaned)


def design_one(seq_str, start, end, amp, cons_info=None, top_n=8):
    """Tek bir hedef bölge için primer+prob tasarımı."""
    conservation = cons_info["conservation"] if cons_info else None
    sets, strict = design_primers_best(
        seq_str, start, end, amp, top_n=top_n, conservation=conservation
    )
    if not sets:
        return None

    best = dict(sets[0])
    offset = max(0, best["forward_start"] - start - 1)
    amp_cons = None
    if conservation:
        amp_cons = conservation[offset:offset + best["amplicon_size"]]
    probe = design_probe(
        best["amplicon_seq"],
        best["forward_len"],
        best["reverse_len"],
        conservation=amp_cons,
    )
    return {"best": best, "probe": probe, "candidates": sets, "strict": strict}


def cross_reactivity(results, sequences):
    """Panel primerlerini diğer etken genomlarında arar."""
    rows = []
    for res in results:
        best = res.get("best") or {}
        for direction, primer in (("Forward", best.get("forward")), ("Reverse", best.get("reverse"))):
            if not primer:
                continue
            p = primer.upper()
            rc = revcomp(p)
            for target_virus, target_seq in sequences.items():
                if target_virus == res["virus"]:
                    continue
                seq = str(target_seq).upper()
                if p in seq or rc in seq:
                    match, risk = "Tam eşleşme", "YÜKSEK"
                elif p[-15:] in seq or rc[-15:] in seq:
                    match, risk = "3' uç 15-mer", "ORTA"
                elif p[:15] in seq or rc[:15] in seq:
                    match, risk = "5' uç 15-mer", "DÜŞÜK"
                else:
                    continue
                rows.append({
                    "primer_of": res["virus"],
                    "direction": direction,
                    "primer": p,
                    "target": target_virus,
                    "match": match,
                    "risk": risk,
                })
    order = {"YÜKSEK": 0, "ORTA": 1, "DÜŞÜK": 2}
    rows.sort(key=lambda r: order.get(r["risk"], 3))
    return rows


def result_to_export(result, meta):
    """UI sonucunu Excel raporunun beklediği yapıya çevirir."""
    return {
        "virus": meta["virus"],
        "segment": meta["segment"],
        "accession": meta["accession"],
        "reference": meta.get("reference", ""),
        "region_text": meta["region_text"],
        "seq_len": meta.get("seq_len"),
        "strict": result["strict"],
        "best": result["best"],
        "probe": result["probe"],
        "cons_info": result.get("cons_info"),
        "candidates": result["candidates"],
    }


def excel_download(results, cross_rows, settings, key):
    data = build_panel_workbook(results, cross_rows=cross_rows, settings=settings)
    st.download_button(
        "📥 Excel raporunu indir (.xlsx)",
        data=data,
        file_name=f"pcr_tasarim_raporu_{datetime.now():%Y%m%d_%H%M}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key=key,
    )


def candidate_table(candidates):
    rows = []
    for i, c in enumerate(candidates, start=1):
        rows.append({
            "#": i,
            "Forward": c["forward"],
            "F Tm": c["forward_tm"],
            "F GC%": c["forward_gc"],
            "Reverse": c["reverse"],
            "R Tm": c["reverse_tm"],
            "R GC%": c["reverse_gc"],
            "Amplikon (bp)": c["amplicon_size"],
            "Skor": c["score"],
        })
    return pd.DataFrame(rows)


# --- 4. KENAR ÇUBUĞU ---
with st.sidebar:
    st.header("⚙️ Ayarlar")

    st.subheader("Virüs Veritabanı")
    with st.expander("➕ Yeni Virüs Ekle"):
        new_name = st.text_input("Virüs Adı")
        new_segment = st.text_input("Hedef Segment")
        new_id = st.text_input("NCBI GenBank ID")
        new_region = st.text_input("Hedef Bölge (örn: 1000-1200)")
        if st.button("Virüsü Ekle", width="stretch"):
            if new_name and new_id:
                VIRUS_DB[new_name] = {
                    "segment": new_segment or "Belirtilmedi",
                    "genbank_id": new_id,
                    "target_region": new_region or "",
                    "reference": "Kullanıcı tanımlı",
                    "database": "Custom",
                }
                st.success(f"{new_name} eklendi.")
                st.rerun()
            else:
                st.error("Virüs adı ve GenBank ID zorunlu.")

    st.markdown("---")
    st.subheader("Tasarım Parametreleri")
    target_amplicon = st.slider("Hedef amplikon boyutu (bp)", 100, 300, 150, 10)
    n_candidates = st.slider("Gösterilecek aday sayısı", 3, 15, 8, 1)

    st.markdown("---")
    st.subheader("Korunmuşluk Analizi")
    use_cons = st.checkbox("Çoklu hizalama ile korunmuşluk hesapla", value=False)
    n_seqs = st.slider(
        "Karşılaştırılacak dizi sayısı", 3, 25, 8, 1,
        help="NCBI'dan çekilen ilgili diziler; yüksek değer daha yavaş ama daha güvenilir.",
    )
    min_cons = st.slider(
        "Minimum korunmuşluk (%)", 70, 100, 85, 1,
        help="Primer bölgesi bu oranın altındaysa aday elenir.",
    ) / 100.0

    st.markdown("---")
    st.subheader("NCBI Bağlantısı")
    if st.button("Bağlantıyı test et", width="stretch"):
        ok, info = check_ncbi_connection()
        st.success(f"✅ Bağlantı aktif ({info})") if ok else st.error(f"❌ {info}")
    st.caption("Erişim yoksa 'Manuel dizi' modunu kullanın.")

# --- 5. BAŞLIK ---
st.title("🧬 Real-Time PCR Primer & Prob Tasarım Motoru")
st.markdown(
    "NCBI GenBank'tan canlı dizi çekme, **çoklu hizalama ile korunmuş bölge tespiti**, "
    "TaqMan primer/prob tasarımı ve **multiplex panel** raporlaması."
)

tab_single, tab_panel, tab_docs = st.tabs(
    ["🧬 Tek Etken Tasarımı", "🧩 Multiplex Panel", "ℹ️ Yöntem ve Notlar"]
)

# --- 6. TEK ETKEN ---
with tab_single:
    if not VIRUS_DB:
        st.warning("Veritabanı boş. Kenar çubuğundan virüs ekleyin.")
    else:
        virus = st.selectbox("Etken seçin", list(VIRUS_DB.keys()), key="single_virus")
        info = VIRUS_DB[virus]

        c1, c2, c3 = st.columns(3)
        c1.metric("Hedef Segment", info["segment"])
        c2.metric("GenBank ID", info["genbank_id"])
        c3.metric("Önerilen Bölge", info["target_region"] or "—")

        st.info(
            f"**Referans:** {info.get('reference', '—')}  \n"
            f"**Bölge:** {info['target_region'] or 'belirtilmedi (tüm dizi taranır)'}  \n"
            f"**Ampikon hedefi:** {target_amplicon} bp"
        )

        source = st.radio(
            "Dizi kaynağı",
            ["NCBI GenBank (canlı)", "Manuel dizi / FASTA"],
            horizontal=True,
            key="single_source",
        )

        manual_text = ""
        manual_region = ""
        manual_organism = ""
        if source.startswith("Manuel"):
            manual_text = st.text_area(
                "Dizi (FASTA başlığı veya ham dizi)",
                height=150,
                placeholder=">sekans_1\nATGCGT...",
                key="manual_seq",
            )
            mc1, mc2 = st.columns(2)
            manual_region = mc1.text_input("Taranacak bölge (örn: 1-300, boş = tümü)", key="manual_region")
            manual_organism = mc2.text_input(
                "Organizma adı (korunmuşluk için, opsiyonel)",
                placeholder="ör. Crimean-Congo hemorrhagic fever virus",
                key="manual_org",
            )

        go = st.button("🚀 Tasarımı çalıştır", type="primary", width="stretch", key="single_run")

        if go:
            with st.spinner("Dizi hazırlanıyor, primerler tasarlanıyor..."):
                if source.startswith("Manuel"):
                    seq_str = clean_fasta(manual_text)
                    if len(seq_str) < 60:
                        st.error("❌ Geçerli bir dizi girin (en az 60 baz, sadece ACGTN).")
                        st.stop()
                    accession = "manuel"
                    organism = manual_organism.strip()
                    region_text = manual_region.strip() or "tüm dizi"
                    reference = "Kullanıcı girdisi"
                else:
                    record = fetch_record(info["genbank_id"])
                    if record is None:
                        st.error(
                            "❌ NCBI'dan dizi alınamadı. Bağlantıyı kontrol edin veya "
                            "'Manuel dizi / FASTA' moduna geçin."
                        )
                        st.stop()
                    seq_str = str(record.seq).upper()
                    accession = info["genbank_id"]
                    organism = record.annotations.get("organism", "")
                    region_text = info["target_region"] or "tüm dizi"
                    reference = info.get("reference", "")

                start, end = parse_region(region_text if region_text != "tüm dizi" else "", len(seq_str))
                keywords = ()
                if source.startswith("Manuel"):
                    keywords = ()
                else:
                    try:
                        keywords = tuple(region_keywords(record, start, end))
                    except Exception:
                        keywords = ()
                st.caption(
                    f"Dizi: **{len(seq_str):,} baz** · Taranan bölge: **{start + 1}-{end}** "
                    f"({end - start} baz) · Organizma: {organism or 'bilinmiyor'}"
                    + (f" · Gen: {', '.join(keywords)}" if keywords else "")
                )

                cons_info = None
                if use_cons:
                    if not organism:
                        st.warning("⚠️ Organizma adı olmadan korunmuşluk analizi yapılamaz; atlandı.")
                    else:
                        with st.spinner(f"İlgili diziler çekiliyor ve hizalanıyor ({n_seqs} dizi)..."):
                            cons_info = get_conservation(
                                accession, seq_str, start, end, organism, n_seqs, keywords
                            )
                        if cons_info["n_used"] == 0:
                            st.warning(
                                "⚠️ Hizalanabilir ilgili dizi bulunamadı; korunmuşluk filtresi uygulanmadı."
                            )
                            cons_info = None
                        else:
                            st.success(
                                f"✅ {cons_info['n_used']} dizi hizalandı · ortalama korunmuşluk "
                                f"%{cons_info['mean'] * 100:.1f}"
                            )

                with st.spinner("Primer ve prob adayları hesaplanıyor..."):
                    result = design_one(
                        seq_str, start, end, target_amplicon,
                        cons_info=cons_info, top_n=n_candidates,
                    )

                if result is None:
                    st.error(
                        "❌ Bu bölgede kriterlere uyan primer çifti bulunamadı. "
                        "Bölgeyi genişletin, amplikon boyutunu değiştirin veya korunmuşluk eşiğini düşürün."
                    )
                    st.stop()

                result["cons_info"] = cons_info
                best = result["best"]
                probe = result["probe"]

                if not result["strict"]:
                    st.warning(
                        "⚠️ Katı kriterleri (Tm 58-65 °C, GC %40-60, dimer yok, 3' G/C) sağlayan aday yok. "
                        "Aşağıdakiler puanlanmış en iyi adaylardır; laboratuvarda optimizasyon gerekebilir."
                    )
                else:
                    st.success("✅ Katı kriterleri sağlayan adaylar bulundu.")

                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Forward Tm", f"{best['forward_tm']} °C")
                m2.metric("Reverse Tm", f"{best['reverse_tm']} °C")
                m3.metric("Amplikon", f"{best['amplicon_size']} bp")
                m4.metric("Prob Tm", f"{probe['tm']} °C" if probe else "—")

                col1, col2, col3 = st.columns(3)
                with col1:
                    st.markdown("**Forward primer (5'→3')**")
                    st.code(best["forward"], language="text")
                    st.caption(
                        f"{best['forward_len']} nt · Tm {best['forward_tm']} °C · "
                        f"GC %{best['forward_gc']} · pozisyon {best['forward_start']}"
                    )
                with col2:
                    st.markdown("**Reverse primer (5'→3')**")
                    st.code(best["reverse"], language="text")
                    st.caption(
                        f"{best['reverse_len']} nt · Tm {best['reverse_tm']} °C · "
                        f"GC %{best['reverse_gc']} · bitiş {best['reverse_end']}"
                    )
                with col3:
                    if probe:
                        st.markdown("**TaqMan prob (5'→3')**")
                        st.code(probe["sequence"], language="text")
                        st.caption(
                            f"{probe['length']} nt · Tm {probe['tm']} °C · GC %{probe['gc']} · "
                            f"pozisyon {probe['position']}"
                        )
                        st.caption("Önerilen etiket: [FAM] … [BHQ1]")
                    else:
                        st.warning("Bu amplikon için uygun prob bulunamadı (SYBR Green kullanılabilir).")

                with st.expander("🧫 Amplikon dizisi ve doğrulama detayları"):
                    st.code(best["amplicon_seq"], language="text")
                    st.caption(
                        f"Amplikon {best['forward_start']}–{best['reverse_end']} nt "
                        f"({best['amplicon_size']} bp)"
                    )
                    ok, msg = validate_primer_primer(best["forward"], best["reverse"])
                    st.write(("✅ " if ok else "❌ ") + msg)
                    if probe:
                        gap = probe["tm"] - max(best["forward_tm"], best["reverse_tm"])
                        st.write(
                            f"Prob–primer Tm farkı: {gap:+.1f} °C "
                            f"({'uygun' if gap >= 3 else 'düşük; probu uzatmayı deneyin'})"
                        )
                    st.write(
                        "Bu çıktı in silico ön değerlendirmedir; Primer-BLAST ve dimer analizi "
                        "ile doğrulanmalıdır."
                    )

                if cons_info:
                    st.markdown("#### 🛡️ Korunmuşluk profili")
                    st.caption(
                        f"{cons_info['n_used']} dizi · ortalama %{cons_info['mean'] * 100:.1f} · "
                        f"pencere ortalaması (30 nt)"
                    )
                    chart_df = pd.DataFrame(
                        {"Korunmuşluk (%)": [v * 100 for v in window_means(cons_info["conservation"], 30)]}
                    )
                    st.line_chart(chart_df, height=180)
                    windows = best_conserved_windows(cons_info["conservation"], 60, 3)
                    if windows:
                        st.write(
                            "En korunmuş 60 nt'lik pencereler: "
                            + ", ".join(f"{start + s + 1}-{start + e}" for s, e, _ in windows)
                        )
                    with st.expander("Hizalanan diziler"):
                        st.dataframe(
                            pd.DataFrame(cons_info["details"]),
                            width="stretch", hide_index=True,
                        )

                st.markdown("#### 🏆 Aday setler")
                st.dataframe(candidate_table(result["candidates"]), width="stretch", hide_index=True)

                st.markdown("#### 📥 Rapor")
                export = [result_to_export(result, {
                    "virus": virus,
                    "segment": info["segment"],
                    "accession": accession,
                    "reference": reference,
                    "region_text": f"{start + 1}-{end}",
                    "seq_len": len(seq_str),
                })]
                excel_download(
                    export, None,
                    {"amplicon": target_amplicon, "use_cons": bool(cons_info), "n_seqs": n_seqs},
                    key="single_xlsx",
                )

# --- 7. MULTIPLEX PANEL ---
with tab_panel:
    st.markdown(
        "Aynı çalışmada kullanılacak birden fazla etken için primer/prob setlerini tasarlar, "
        "ardından **primerlerin diğer etken genomlarına bağlanma riskini** tarar."
    )

    selected = st.multiselect(
        "Panele eklenecek etkenler",
        list(VIRUS_DB.keys()),
        default=list(VIRUS_DB.keys())[:3],
    )
    st.caption(
        f"Amplikon hedefi: {target_amplicon} bp · Korunmuşluk analizi: "
        f"{'açık (' + str(n_seqs) + ' dizi)' if use_cons else 'kapalı'}"
    )

    run_panel = st.button(
        "🧩 Paneli tasarla ve çapraz reaksiyonu tara",
        type="primary", width="stretch", key="panel_run",
    )

    if run_panel:
        if len(selected) < 2:
            st.error("Panel için en az 2 etken seçin.")
            st.stop()

        results, sequences = [], {}
        failed = []
        progress = st.progress(0.0, text="Başlıyor...")

        for idx, virus in enumerate(selected, start=1):
            progress.progress((idx - 1) / len(selected), text=f"{virus} işleniyor...")
            info = VIRUS_DB[virus]
            record = fetch_record(info["genbank_id"])
            if record is None:
                failed.append((virus, "NCBI'dan dizi alınamadı"))
                continue
            seq_str = str(record.seq).upper()
            sequences[virus] = seq_str
            start, end = parse_region(info["target_region"], len(seq_str))

            cons_info = None
            if use_cons:
                organism = record.annotations.get("organism", "")
                if organism:
                    try:
                        keywords = tuple(region_keywords(record, start, end))
                    except Exception:
                        keywords = ()
                    cons_info = get_conservation(
                        info["genbank_id"], seq_str, start, end, organism, n_seqs, keywords
                    )
                    if cons_info["n_used"] == 0:
                        cons_info = None

            result = design_one(
                seq_str, start, end, target_amplicon, cons_info=cons_info, top_n=n_candidates
            )
            if result is None:
                failed.append((virus, "kriterlere uyan primer bulunamadı"))
                continue
            result["cons_info"] = cons_info
            results.append(result_to_export(result, {
                "virus": virus,
                "segment": info["segment"],
                "accession": info["genbank_id"],
                "reference": info.get("reference", ""),
                "region_text": f"{start + 1}-{end}",
                "seq_len": len(seq_str),
            }))

        progress.progress(1.0, text="Tamamlandı")
        progress.empty()

        for virus, reason in failed:
            st.warning(f"⚠️ {virus}: {reason}")

        if not results:
            st.error("Hiçbir etken için tasarım üretilemedi.")
            st.stop()

        # Özet tablo
        summary_rows = []
        for res in results:
            best = res["best"]
            probe = res["probe"] or {}
            summary_rows.append({
                "Etken": res["virus"],
                "Accession": res["accession"],
                "Bölge": res["region_text"],
                "F Tm": best["forward_tm"],
                "R Tm": best["reverse_tm"],
                "Prob Tm": probe.get("tm", "—"),
                "Amplikon (bp)": best["amplicon_size"],
                "Korunmuşluk %": round(res["cons_info"]["mean"] * 100, 1) if res["cons_info"] else "—",
                "Mod": "Katı" if res["strict"] else "Gevşek",
            })
        st.markdown("#### Panel özeti")
        st.dataframe(pd.DataFrame(summary_rows), width="stretch", hide_index=True)

        # Dizi tabloları
        seq_rows = []
        for res in results:
            best = res["best"]
            probe = res["probe"]
            seq_rows.append({
                "Etken": res["virus"],
                "Forward (5'→3')": best["forward"],
                "Reverse (5'→3')": best["reverse"],
                "Prob (5'→3')": probe["sequence"] if probe else "—",
            })
        st.markdown("#### Oligonükleotid setleri")
        st.dataframe(pd.DataFrame(seq_rows), width="stretch", hide_index=True)

        # Çapraz reaksiyon
        st.markdown("#### 🔍 Çapraz reaksiyon taraması")
        cross_rows = cross_reactivity(results, sequences)
        if cross_rows:
            st.dataframe(
                pd.DataFrame([{
                    "Primer (etken)": f"{r['primer_of']} ({r['direction']})",
                    "Dizi": r["primer"],
                    "Aranan genom": r["target"],
                    "Eşleşme": r["match"],
                    "Risk": r["risk"],
                } for r in cross_rows]),
                width="stretch", hide_index=True,
            )
            high = [r for r in cross_rows if r["risk"] == "YÜKSEK"]
            if high:
                st.error(
                    f"❌ {len(high)} yüksek riskli eşleşme var. Bu primerler diğer etkenlerde "
                    "amplifikasyona yol açabilir; alternatif aday seçin."
                )
            else:
                st.info("Yüksek riskli tam eşleşme yok; orta/düşük riskler 3' uç optimizasyonu ile azaltılabilir.")
        else:
            st.success("✅ Seçilen etkenler arasında primer eşleşmesi bulunmadı.")

        st.markdown("#### 📥 Rapor")
        excel_download(
            results, cross_rows,
            {
                "amplicon": target_amplicon,
                "use_cons": use_cons,
                "n_seqs": n_seqs,
                "panel": ", ".join(selected),
            },
            key="panel_xlsx",
        )

# --- 8. YÖNTEM ---
with tab_docs:
    st.markdown(
        """
### Tasarım kuralları
| Parametre | Aralık |
|---|---|
| Primer uzunluğu | 18–25 nt |
| Primer Tm | 58–65 °C (hedef 61 °C), F/R farkı ≤ 2 °C |
| Primer GC | %40–60 |
| 3' uç | G veya C tercih edilir |
| Homopolymer | ≥ 6 aynı baz istenmez |
| Dimer | 5+ baz tamamlayıcılık yüksek risk, 4 baz ceza puanı |
| Prob | 18–30 nt, Tm 66–73 °C (~69 °C), 5' ucu G değil |
| Amplikon | 100–300 bp (seçilen hedef ± 35 bp) |

### Tm hesabı
SantaLucia nearest-neighbour modeli (Biopython `MeltingTemp.Tm_NN`), qPCR koşulları:
50 mM Na⁺, 1.5 mM Mg²⁺, 0.2 mM dNTP, 250 nM oligonükleotid.

### Korunmuşluk analizi
1. Referans bölge + 120 nt kanat alınır.
2. NCBI'dan aynı organizmaya ait diziler çekilir.
3. Her dizide 24-mer çapa aranır, bulunan pencere global olarak hizalanır (bulunamazsa lokal hizalama).
4. Hizalama blokları referans koordinatlarına iz düşürülür; her pozisyon için
   "aynı bazda olan dizi oranı" hesaplanır.
5. Primer bölgeleri eşiğin (varsayılan %85) altındaysa aday elenir; korunmuşluk skora bonus ekler.

### Sınırlamalar
- Bu araç **in silico ön değerlendirme** yapar; klinik validasyon yerine geçmez.
- Özgüllük için mutlaka **Primer-BLAST / BLASTn** ile doğrulama yapın.
- Çapraz reaksiyon taraması k-mer tabanlıdır; tam hizalama değildir.
- NCBI erişimi olmadan "Manuel dizi" modunu kullanın (korunmuşluk analizi organizma adı gerektirir).
        """
    )

st.caption(f"v2.0 · {datetime.now():%Y-%m-%d %H:%M} · in silico ön değerlendirme aracı")
