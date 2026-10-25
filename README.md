# Real-Time PCR Primer & Prob Tasarım Motoru

NCBI GenBank'tan canlı dizi çeken, çoklu hizalama ile korunmuş bölge tespit eden ve
TaqMan primer/prob setleri tasarlayan Streamlit uygulaması.

## Özellikler

- **Canlı NCBI erişimi** (Biopython Entrez) — her çalıştırmada güncel diziler
- **Korunmuşluk analizi** — k-mer ön filtre + global/lokal hizalama, pozisyon bazlı korunmuşluk
- **TaqMan tasarımı** — SantaLucia nearest-neighbour Tm, dimer/homopolymer kontrolü
- **Özgüllük filtresi** — tek etken tasarımında diğer katalog virüslerine karşı çapraz reaksiyon elemesi
- **Prob–primer dimer kontrolü** — ≥6 baz tamamlayıcılık yapan prob adayları elenir
- **Multiplex panel** — çoklu etken + primerler arası çapraz reaksiyon taraması
- **Excel raporu** — 7 sayfa: özet, primerler, problar, korunmuşluk, çapraz reaksiyon, protokol, adaylar
- **Manuel dizi modu** — internet yoksa FASTA/ham dizi ile çalışır
- **Etken kataloğu** — virüs, bakteri ve parazit; büyük genomlarda hedef gen taraması
- **Hazır etken şablonu** — yeni etken/varyant formunda listeden seçince alanlar otomatik dolar

## Hedef etkenler

**Virüsler (9):** CCHFV · Batı Nil · Hantavirüs · Dengue · Chikungunya · Zika · Ebola · Marburg · Aigai

**Bakteriler (21):** *M. tuberculosis* · *E. coli* · *S. aureus* · *S. pneumoniae* ·
*K. pneumoniae* · *Salmonella* · *L. monocytogenes* · *V. cholerae* · *N. meningitidis* ·
*B. pertussis* · *L. pneumophila* · *H. pylori* · *C. difficile* · *Brucella* ·
*F. tularensis* · *Y. pestis* · *B. anthracis* · *C. trachomatis* · *M. pneumoniae* ·
*C. burnetii* · *P. aeruginosa*

**Parazitler (10):** *P. falciparum* · *T. gondii* · *L. donovani* · *T. brucei* · *T. cruzi* ·
*Giardia* · *E. histolytica* · *C. parvum* · *S. mansoni* · *E. multilocularis*

Kenar çubuğundan yeni etken eklenebilir; kalıcı ekleme için `pathogens.py` içindeki
`PATHOGENS` sözlüğüne kayıt ekleyin.

## Yerel çalıştırma

```bash
pip install -r requirements.txt
streamlit run realtime_pcr_designer.py
```

Windows'ta `BASLAT.bat` dosyasına çift tıklamak yeterli.

## Streamlit Cloud'da yayınlama

1. Bu depoyu GitHub'a yükle
2. https://share.streamlit.io → **New app** → depo + `realtime_pcr_designer.py`
3. **Advanced settings → Secrets** alanına ekle:

```toml
NCBI_EMAIL = "ad.soyad@kurum.gov.tr"
```

4. Deploy — uygulama `https://xxx.streamlit.app` adresinde yayınlanır.

## Uyarı

Bu araç **in silico ön değerlendirme** içindir. Klinik/laboratuvar kullanımı öncesi
Primer-BLAST özgüllük kontrolü, dimer analizi ve LOD/LOQ validasyonu zorunludur.
