# Real-Time PCR Primer & Prob Tasarım Motoru

NCBI GenBank'tan canlı dizi çeken, çoklu hizalama ile korunmuş bölge tespit eden ve
TaqMan primer/prob setleri tasarlayan Streamlit uygulaması.

## Özellikler

- **Koyu lacivert tema** — göz yormayacak koyu lacivert zemin ve açık metin
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

### v3.0 analiz ve ileri tasarım sekmeleri

- **Termodinamik** — Primer3 hairpin / homo-heterodimer ΔG, 3' uç stabilitesi, dimer matrisi, Tm Monte Carlo belirsizliği
- **Özgüllük** — deterministik in-silico PCR, panel taraması, risk skoru, isteğe bağlı NCBI BLAST
- **Genotip kapsama** — çoklu dizi hizalaması, mismatch matrisi, Wilson %95 GA, sequence logo, genome track
- **LOD / LOQ** — Poisson + verim + Ct gürültüsü Monte Carlo, probit MLE/bootstrap, standart eğri, PPV/NPV, örneklem büyüklüğü
- **PDF assay kartı** — termodinamik, kapsama, özgüllük ve LOD/LOQ özeti
- **Literatür / assay karşılaştırma** — yayınlanmış primerleri referansta konumlandırma, aktif assay ile uyum skoru, JSON kütüphanesi
- **Yayınlanmış Assay'ler** — Semantic Scholar/Crossref üzerinden makale bibliyografisi arama; sonuçlar atıf sayısına göre sıralanır. Makaledeki primer/prob dizileri aktarılıp referans üzerinde konumlandırılır ve kütüphaneye eklenir; sistem yayınlanmış dizi uydurmaz.
- **Consensus / dejenere primer** — IUPAC dejenere primerler, kapsama ve dejenere kat sayısı
- **SNP ayrımı** — 3' uç SNP primerleri + WT kontrol primerı
- **RPA / LAMP (beta)** — RPA primer/prob ve LAMP F3/B3/FIP/BIP önerileri, LAMP set kontrolü
- **Panel Optimizer** — assay'ler arası dimer matrisi, uyumluluk matrisi, kanal ataması

### 🛡️ Sipariş Kontrolü (Türkiye suşları + sipariş öncesi QC)

- **Türkiye Suş Kontrolü** — NCBI `Turkey[Country]` araması, yerel FASTA/ham dizi yükleme;
  primer/prob bağlanma bölgelerinde mismatch, 3' uç mismatch'i, prob mismatch'i ve amplikon
  varlığını PASS / WARN / FAIL olarak raporlar. Hedef bölgeyi içermeyen kısa/partial
  kayıtlar **N/A** sayılır ve mismatch olarak değerlendirilmez.
- **Sipariş Öncesi Kontrol** — termodinamik, özgüllük taraması, Türkiye suşları, Primer-BLAST,
  sentez/modifikasyon ve ıslak laboratuvar planı tek karar tablosunda birleştirilir.
  FAIL/BLOCKED varsa sipariş listesi indirilemez; WARN varsa açık onay gerekir.
- **Türkiye consensus adayı** — Türkiye dizilerinden yeni aday üretir ve adayları tekrar
  Türkiye dizilerine karşı test eder. Non-dejenere aday tüm dizileri geçmiyorsa "sipariş
  etmeyin" uyarısı verir.

## Hedef etkenler

**Virüsler (11):** CCHFV · Batı Nil · Hantavirüs · Dengue · Chikungunya · Zika · Ebola · Marburg · Aigai ·
Tatarcık humması (Toscana/SFNV) · Kene kaynaklı ensefalit (TBEV)

**Bakteriler (28):** *M. tuberculosis* · *E. coli* · *S. aureus* · *S. pneumoniae* ·
*K. pneumoniae* · *Salmonella* · *L. monocytogenes* · *V. cholerae* · *N. meningitidis* ·
*B. pertussis* · *L. pneumophila* · *H. pylori* · *C. difficile* · *Brucella* ·
*F. tularensis* · *Y. pestis* · *B. anthracis* · *C. trachomatis* · *M. pneumoniae* ·
*C. burnetii* · *P. aeruginosa* · *B. burgdorferi* · *B. hermsii* · *R. rickettsii* ·
*A. phagocytophilum* · *E. chaffeensis* · *N. mikurensis* · *B. henselae*

**Parazitler (11):** *P. falciparum* · *T. gondii* · *L. donovani* · *T. brucei* · *T. cruzi* ·
*Giardia* · *E. histolytica* · *C. parvum* · *S. mansoni* · *E. multilocularis* · *B. microti*

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
Sipariş öncesi **Türkiye Suş Kontrolü** ve **Sipariş Öncesi Kontrol** kapısı çalıştırılmalıdır;
NCBI'da ülke metadata'sı eksik olabileceği için yerel/Türkiye izolatı FASTA dizileriyle
kontrolü tamamlayın. RPA/LAMP çıktıları **beta**dır; ıslak laboratuvar optimizasyonu olmadan
kullanılmamalıdır. PDF ve kütüphane çıktıları araştırma amaçlıdır, tıbbi karar desteği değildir.

**Oligo doğrulama notu:** uygulama hazır literatür dizileri sunmaz; referans üzerinden yeni aday
üretir. Kısa (19–24 nt) oligo dizilerini Google'da aramak doğrulama değildir. Forward ve prob
referans ileri iplikte aynen bulunur; reverse primerin sipariş dizisi ileri iplikte ters
komplement olarak bulunur. Doğrulama için NCBI BLASTn/Primer-BLAST ve accession + koordinat
kullanılmalıdır.

**Öncelikli arbovirüsler:** CCHFV, Batı Nil, Dengue, Zika ve Chikungunya için katalogdaki
hazır setler doğrudan sipariş edilmemelidir. CCHFV ve Batı Nil'de NCBI Türkiye kayıtları
mevcut olup mevcut setler uyumsuzdur; Dengue, Zika ve Chikungunya için NCBI'da `Turkey[Country]`
kaydı bulunmadığından yerel/importe izolat FASTA dizileriyle kontrol zorunludur.
