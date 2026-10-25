import json
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
from pathogens import PATHOGENS, CATEGORY_ORDER, CATEGORY_ICON

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

if "virus_db" not in st.session_state:
    st.session_state.virus_db = {name: dict(entry) for name, entry in PATHOGENS.items()}
VIRUS_DB = st.session_state.virus_db


def ordered_names(db):
    """Etken adlarını kategoriye göre sıralar (Virüs → Bakteri → Parazit)."""
    names = []
    for category in CATEGORY_ORDER:
        names.extend([n for n, v in db.items() if v.get("category", "Virüs") == category])
    names.extend([n for n in db if n not in names])
    return names


def display_name(db, name):
    category = db[name].get("category", "Virüs")
    return f"{CATEGORY_ICON.get(category, '•')} {name}"


def conservation_organism(organism, category):
    """Korunmuşluk araması için organizma adını uygun düzeye indirir.
    Bakterilerde suş eki (ör. 'NCTC 8325') atılır; virüslerde tam ad korunur."""
    if not organism:
        return organism
    if category == "Bakteri":
        return " ".join(organism.split()[:2])
    return organism


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
def fetch_genome(accession: str):
    """Büyük genomlar için tam dizi (FASTA) + organizma adı (esummary).
    Dönüş: (sequence_str, organism) veya None."""
    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
            record = SeqIO.read(handle, "fasta")
        seq = str(record.seq).upper()
    except Exception:
        return None

    organism = ""
    try:
        with Entrez.esummary(db="nucleotide", id=accession, retmode="json") as handle:
            data = json.loads(handle.read())
        result = data.get("result", {})
        uid = (result.get("uids") or [None])[0]
        if uid and uid in result:
            organism = result[uid].get("organism", "") or ""
    except Exception:
        pass
    return seq, organism


@st.cache_data(ttl=3600, show_spinner=False)
def get_conservation(accession, seq_str, start, end, organism, n_seqs, keywords=(), filter_seq=None):
    """İlgili dizileri çekip pozisyon bazlı korunmuşluk hesaplar."""
    big_genome = len(seq_str) > 200_000
    # Büyük genomlarda k-mer filtresi hedef bölge üzerinden yapılır (hafıza/hız).
    reference_for_filter = filter_seq if filter_seq else seq_str
    subjects = fetch_related_sequences(
        organism,
        exclude_accessions=[accession],
        max_n=int(n_seqs),
        reference_seq=reference_for_filter,
        keywords=tuple(keywords),
        max_candidates=(int(n_seqs) * 3 if big_genome else None),
        min_subject_len=(int((end - start) * 0.5) if big_genome else None),
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
        start, end = 0, min(seq_len, 30000)
    if end - start < 60:
        start, end = 0, min(seq_len, 30000)
    return start, end


def clean_fasta(text):
    """FASTA veya ham diziden sadece ACGTN bırakır."""
    if not text:
        return ""
    cleaned = "".join(
        line.strip() for line in text.splitlines() if not line.strip().startswith(">")
    ).upper()
    return re.sub(r"[^ACGTN]", "", cleaned)


def record_sequence(record):
    """SeqRecord'dan diziyi güvenle çıkarır; tanımsız (CON) dizilerde None döner."""
    if record is None:
        return None
    try:
        seq = str(record.seq).upper()
    except Exception:
        return None
    return seq if seq else None


def design_one(seq_str, start, end, amp, cons_info=None, top_n=8):
    """Tek bir hedef bölge için primer+prob tasarımı."""
    conservation = cons_info["conservation"] if cons_info else None
    sets, strict = design_primers_best(
        seq_str, start, end, amp, top_n=top_n, conservation=conservation
    )
    if not sets:
        return None

    def _probe_for(candidate, **kwargs):
        offset = max(0, candidate["forward_start"] - start - 1)
        amp_cons = None
        if conservation:
            amp_cons = conservation[offset:offset + candidate["amplicon_size"]]
        return design_probe(
            candidate["amplicon_seq"],
            candidate["forward_len"],
            candidate["reverse_len"],
            conservation=amp_cons,
            **kwargs,
        )

    best, probe, relaxed = None, None, False
    for candidate in sets:
        found = _probe_for(candidate)
        if found:
            best, probe = dict(candidate), found
            break

    if best is None:
        # AT-zengin bölgelerde (ör. parazit mitokondrisi) prob kriterlerini gevşet
        best = dict(sets[0])
        probe = _probe_for(best, tm_min=60.0, tm_max=75.0, gc_min=25.0, gc_max=70.0)
        relaxed = probe is not None

    return {
        "best": best,
        "probe": probe,
        "probe_relaxed": relaxed,
        "candidates": sets,
        "strict": strict,
    }


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

    st.subheader("Etken Kataloğu")
    st.caption(f"{len(VIRUS_DB)} etken yüklü · Virüs / Bakteri / Parazit")
    with st.expander("➕ Yeni Etken Ekle"):
        new_name = st.text_input("Etken adı", key="new_name")
        new_category = st.selectbox("Kategori", CATEGORY_ORDER, key="new_category")
        new_segment = st.text_input("Hedef gen / bölge açıklaması", key="new_segment")
        new_id = st.text_input("NCBI GenBank ID", key="new_id")
        new_region = st.text_input("Hedef bölge (örn: 1000-1200)", key="new_region")
        new_gene = st.text_input(
            "Hedef gen adı (opsiyonel)", key="new_gene",
            help="Büyük genomlarda (bakteri/parazit) gen adı + koordinat verin.",
        )
        gc1, gc2 = st.columns(2)
        new_gene_start = gc1.number_input("Gen başlangıcı", min_value=0, value=0, step=1, key="new_gstart")
        new_gene_end = gc2.number_input("Gen bitişi", min_value=0, value=0, step=1, key="new_gend")
        if st.button("Etkeni Ekle", width="stretch", key="add_pathogen"):
            if new_name and new_id:
                VIRUS_DB[new_name] = {
                    "category": new_category,
                    "segment": new_segment or "Belirtilmedi",
                    "genbank_id": new_id,
                    "target_region": new_region or "",
                    "target_gene": new_gene.strip(),
                    "gene_start": int(new_gene_start) or None,
                    "gene_end": int(new_gene_end) or None,
                    "reference": "Kullanıcı tanımlı",
                }
                st.success(f"{new_name} eklendi.")
                st.rerun()
            else:
                st.error("Etken adı ve GenBank ID zorunlu.")
    st.download_button(
        "⬇️ Kataloğu indir (JSON)",
        data=json.dumps(VIRUS_DB, ensure_ascii=False, indent=1),
        file_name="etken_katalogu.json",
        mime="application/json",
        key="catalog_json",
    )
    st.caption("Kenar çubuğundan eklenenler oturum boyunca kalır; kalıcı olması için `pathogens.py` dosyasına ekleyin.")

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
    "TaqMan primer/prob tasarımı ve **multiplex panel** raporlaması.  \n"
    "Katalogda **9 virüs, 21 bakteri ve 10 parazit** hazır gelir; kenar çubuğundan yeni etken ekleyebilirsiniz."
)

tab_single, tab_panel, tab_docs = st.tabs(
    ["🧬 Tek Etken Tasarımı", "🧩 Multiplex Panel", "ℹ️ Yöntem ve Notlar"]
)

# --- 6. TEK ETKEN ---
with tab_single:
    if not VIRUS_DB:
        st.warning("Veritabanı boş. Kenar çubuğundan virüs ekleyin.")
    else:
        virus = st.selectbox(
            "Etken seçin",
            ordered_names(VIRUS_DB),
            format_func=lambda n: display_name(VIRUS_DB, n),
            key="single_virus",
        )
        info = VIRUS_DB[virus]
        has_gene = bool(info.get("gene_start") and info.get("gene_end"))
        target_label = (
            f"{info['gene_start']}-{info['gene_end']} (gen: {info.get('target_gene', '')})"
            if has_gene
            else (info["target_region"] or "tüm dizi")
        )

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Kategori", info.get("category", "—"))
        c2.metric("Hedef", info["segment"])
        c3.metric("GenBank ID", info["genbank_id"])
        c4.metric("Taranacak bölge", target_label)

        st.info(
            f"**Referans:** {info.get('reference', '—')}  \n"
            f"**Hedef:** {target_label}  \n"
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
                    accession = info["genbank_id"]
                    reference = info.get("reference", "")
                    if has_gene:
                        fetched = fetch_genome(accession)
                        if fetched is None:
                            st.error(
                                "❌ NCBI'dan genom alınamadı. Bağlantıyı kontrol edin veya "
                                "'Manuel dizi / FASTA' moduna geçin."
                            )
                            st.stop()
                        seq_str, organism = fetched
                        region_text = f"{info['gene_start']}-{info['gene_end']}"
                        record = None
                    else:
                        record = fetch_record(accession)
                        seq_str = record_sequence(record)
                        if seq_str is None:
                            fetched = fetch_genome(accession)
                            if fetched is None:
                                st.error(
                                    "❌ NCBI'dan dizi alınamadı. Bağlantıyı kontrol edin veya "
                                    "'Manuel dizi / FASTA' moduna geçin."
                                )
                                st.stop()
                            seq_str, organism = fetched
                        else:
                            organism = record.annotations.get("organism", "")
                        region_text = info["target_region"] or "tüm dizi"

                start, end = parse_region(region_text if region_text != "tüm dizi" else "", len(seq_str))
                if not source.startswith("Manuel") and has_gene:
                    start = max(0, int(info["gene_start"]) - 1)
                    end = min(len(seq_str), int(info["gene_end"]))
                big_genome = len(seq_str) > 200_000
                keywords = ()
                if info.get("target_gene"):
                    keywords = (info["target_gene"],)
                elif not source.startswith("Manuel") and record is not None:
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
                        n_eff = min(n_seqs, 4) if big_genome else n_seqs
                        if big_genome and n_eff != n_seqs:
                            st.info(
                                f"ℹ️ Büyük genom: korunmuşluk için dizi sayısı {n_eff} ile sınırlandı."
                            )
                        with st.spinner(f"İlgili diziler çekiliyor ve hizalanıyor ({n_eff} dizi)..."):
                            cons_info = get_conservation(
                                accession, seq_str, start, end,
                                conservation_organism(organism, info.get("category", "Virüs")),
                                n_eff, keywords,
                                filter_seq=seq_str[start:end] if big_genome else None,
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
                        if result.get("probe_relaxed"):
                            st.warning(
                                "ℹ️ Bu bölge AT-zengin olduğu için prob Tm/GC hedefleri gevşetildi "
                                f"(Tm {probe['tm']} °C, GC %{probe['gc']}). Laboratuvarda prob "
                                "konsantrasyonunu buna göre ayarlayın."
                            )
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
        ordered_names(VIRUS_DB),
        default=ordered_names(VIRUS_DB)[:3],
        format_func=lambda n: display_name(VIRUS_DB, n),
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
            has_gene = bool(info.get("gene_start") and info.get("gene_end"))
            if has_gene:
                fetched = fetch_genome(info["genbank_id"])
                if fetched is None:
                    failed.append((virus, "NCBI'dan genom alınamadı"))
                    continue
                seq_str, organism = fetched
                record = None
                start = max(0, int(info["gene_start"]) - 1)
                end = min(len(seq_str), int(info["gene_end"]))
                region_label = f"{start + 1}-{end}"
            else:
                record = fetch_record(info["genbank_id"])
                seq_str = record_sequence(record)
                if seq_str is None:
                    fetched = fetch_genome(info["genbank_id"])
                    if fetched is None:
                        failed.append((virus, "NCBI'dan dizi alınamadı"))
                        continue
                    seq_str, organism = fetched
                else:
                    organism = record.annotations.get("organism", "")
                start, end = parse_region(info["target_region"], len(seq_str))
                region_label = f"{start + 1}-{end}"

            sequences[virus] = seq_str

            cons_info = None
            if use_cons and organism:
                big = len(seq_str) > 200_000
                n_eff = min(n_seqs, 4) if big else n_seqs
                keywords = ()
                if info.get("target_gene"):
                    keywords = (info["target_gene"],)
                elif record is not None:
                    try:
                        keywords = tuple(region_keywords(record, start, end))
                    except Exception:
                        keywords = ()
                cons_info = get_conservation(
                    info["genbank_id"], seq_str, start, end,
                    conservation_organism(organism, info.get("category", "Virüs")),
                    n_eff, keywords,
                    filter_seq=seq_str[start:end] if big else None,
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

### Etken kataloğu ve gen hedefleme
- Katalogda virüsler (CCHFV, WNV, Hantavirüs, Dengue, Chikungunya, Zika, Ebola, Marburg, Aigai),
  bakteriler ve parazitler hazır gelir.
- Bakteri/parazit genomları çok büyük olduğu için tüm genom taranmaz: katalogda tanımlı
  **hedef gen koordinatları** (ör. *S. aureus* `nuc`, *M. tuberculosis* `IS6110`) taranır.
- Yeni etken eklerken büyük genomlar için **gen adı + başlangıç/bitiş koordinatı** girin;
  küçük genomlarda (virüs) yalnızca bölge aralığı yeterlidir.
- AT-zengin bölgelerde (ör. *Leishmania* maxicircle) uygun prob bulunamazsa Tm/GC hedefleri
  otomatik gevşetilir ve arayüzde uyarı gösterilir.

### Sınırlamalar
- Bu araç **in silico ön değerlendirme** yapar; klinik validasyon yerine geçmez.
- Özgüllük için mutlaka **Primer-BLAST / BLASTn** ile doğrulama yapın.
- Çapraz reaksiyon taraması k-mer tabanlıdır; tam hizalama değildir.
- NCBI erişimi olmadan "Manuel dizi" modunu kullanın (korunmuşluk analizi organizma adı gerektirir).
        """
    )

st.caption(f"v2.0 · {datetime.now():%Y-%m-%d %H:%M} · in silico ön değerlendirme aracı")
