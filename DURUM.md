# PROJE DURUMU — kaldığım yerden devam

Son güncelleme: 07.11.2026

## Proje nedir

Virüs, bakteri ve parazitler için NCBI'dan canlı dizi çekip TaqMan primer/prob tasarlayan
Streamlit uygulaması. Hazır katalog: **11 virüs, 28 bakteri, 11 parazit (50 etken)**.

## Şu an nerede

| Konu | Durum |
|---|---|
| Kod | `C:\Users\Hakan Hedef\Desktop\PCR_PRIMER_PROB` |
| GitHub | `Hknhdf/PCR_PRIMER_PROB` — `main` dalı |
| Streamlit Cloud | Yayında — `https://pcr-primer-prob.streamlit.app` |
| Yerel uygulama | `BASLAT.bat` → `http://localhost:8501` |
| Sürüm | **v3.0** |

## Tamamlanan özellikler

### v2 temel motoru
- **Koyu lacivert tema:** Streamlit teması, grafik şablonları ve sequence logo zeminleri koyu lacivert
- Canlı NCBI erişimi (Biopython Entrez), önbellekli
- SantaLucia nearest-neighbour Tm (qPCR tuz koşulları)
- Gerçek amplikon mantığı (forward = sense başlangıç, reverse = amplikon sonunun ters komplementi)
- Katı / gevşek iki kademeli filtre, dimer ve homopolymer kontrolü
- **Korunmuşluk analizi:** k-mer ön filtre → 24-mer çapa → global/lokal hizalama → pozisyon bazlı oran
- **Multiplex panel:** çoklu etken, primerler arası çapraz reaksiyon taraması
- **Excel raporu:** Özet, Primerler, Problar, Korunmuşluk, Çapraz Reaksiyon, Protokol, Adaylar
- Manuel FASTA/ham dizi modu (internet yoksa)
- **Etken kataloğu:** virüs / bakteri / parazit, kategoriye göre gruplu seçim
  - Bakteri ve parazitlerde tüm genom yerine **hedef gen** taranır
  - Büyük genomlarda korunmuşluk için dizi sayısı otomatik sınırlandırılır
  - AT-zengin bölgelerde prob Tm/GC hedefleri otomatik gevşetilir
  - *T. gondii* hedefi 529 bp tekrar elementi (AF146527.1)
- **Özgüllük filtresi:** tek etken tasarımında diğer katalog virüslerine bağlanabilecek adaylar elenir
- **Prob–primer dimer filtresi:** prob, primerlerle ≥6 baz tamamlayıcılık yapmayacak şekilde seçilir
- Kenar çubuğundan yeni etken ekleme (oturumluk) + kataloğu JSON olarak indirme
- **Hazır etken şablonu:** “Yeni Etken / Varyant Ekle” formunda listeden etken seçince alanlar otomatik dolar

### v3.0 yeni sekmeler
- **📊 Analiz**
  - **Termodinamik:** Primer3 hairpin, homo/heterodimer ΔG, 3' uç stabilitesi, dimer matrisi,
    Tm Monte Carlo belirsizliği, PDF assay kartı
  - **Özgüllük:** deterministik in-silico PCR, panel taraması, risk skoru, isteğe bağlı NCBI BLAST
  - **Genotip Kapsama:** çoklu dizi hizalaması, mismatch matrisi, Wilson %95 GA, entropi,
    sequence logo, CSV indirme
  - **LOD / LOQ:** Poisson + verim + Ct gürültüsü Monte Carlo, probit MLE/bootstrap,
    standart eğri, PPV/NPV, örneklem büyüklüğü
  - **PDF Rapor:** termodinamik + kapsama + özgüllük + LOD/LOQ özeti
- **📚 Kütüphane**
  - **Literatür / Assay Karşılaştırma:** yayınlanmış primerleri referansta konumlandırma,
    aktif assay ile uyum skoru, JSON kütüphanesi; doğrulanmamış hazır dizi gömülmez
  - **Genome Track:** Plotly tabanlı gen/primer/korunmuşluk/mismatch izleyicisi
- **📄 Yayınlanmış Assay'ler**
  - Semantic Scholar/Crossref ile makale bibliyografisini atıf sayısına göre sıralar; erişilebilir
    tam metin varsa olası primer/prob dizileri otomatik aday olarak çıkarılır. Onaylanan diziler
    referans üzerinde konumlandırılır ve assay kütüphanesine eklenir. Sistem yayınlanmış dizi
    uydurmaz.
- **🚀 İleri**
  - **Consensus / Dejenere:** IUPAC dejenere primerler, kapsama ve dejenere kat sayısı
  - **SNP Ayrımı:** 3' uç SNP primerleri + WT kontrol primerı
  - **LAMP / RPA (beta):** RPA primer/prob, LAMP F3/B3/FIP/BIP, LAMP set kontrolü
  - **Panel Optimizer:** assay'ler arası dimer matrisi, uyumluluk matrisi, kanal ataması
- **🛡️ Sipariş Kontrolü**
  - **Türkiye Suş Kontrolü:** NCBI `Turkey[Country]` araması, yerel FASTA/ham dizi yükleme,
    primer/prob mismatch, 3' uç ve prob mismatch analizi; PASS/WARN/FAIL sınıflandırması
  - **Sipariş Öncesi Kontrol:** termodinamik + özgüllük + Türkiye sonucu + Primer-BLAST +
    sentez teyidi + ıslak laboratuvar planı; FAIL/BLOCKED indirmeyi kilitler, WARN onay ister
  - **Türkiye consensus adayı:** Türkiye dizilerinden yeni aday üretir, adayları tekrar
    Türkiye dizilerine karşı test eder; hiçbiri geçmiyorsa sipariş uyarısı verir

## Kataloğa eklenen yeni etkenler (05.11.2026)

WHO vektör/olası patojen listesini kapsayacak şekilde eklendi. Her kayıt NCBI RefSeq
anotasyonundan alınan hedef gen koordinatlarıyla doğrulandı.

| Etken | Kategori | Erişim | Hedef |
|---|---|---|---|
| Tatarcık humması virüsü (Toscana/SFNV) | Virüs | NC_006318.1 | S segmenti N, 1100–1600 |
| Kene kaynaklı ensefalit virüsü (TBEV) | Virüs | NC_001672.1 | NS5, 9000–9400 |
| *Borrelia burgdorferi* (Lyme) | Bakteri | NC_001318.1 | flaB, 147648–148659 |
| *Borrelia hermsii* (tekrarlayan ateş) | Bakteri | NZ_CP073136.1 | glpQ, 250003–251026 |
| *Rickettsia rickettsii* (benekli ateş) | Bakteri | NC_010263.1 | ompA, 1176529–1181728 |
| *Anaplasma phagocytophilum* | Bakteri | NC_007797.1 | ankA, 784799–788498 |
| *Ehrlichia chaffeensis* | Bakteri | NC_007799.1 | groEL, 353837–355490 |
| *Neoehrlichia mikurensis* | Bakteri | NZ_CP089285.1 | groEL, 158219–159869 |
| *Bartonella henselae* | Bakteri | NC_005956.1 | groEL, 1547626–1549270 |
| *Babesia microti* | Parazit | XR_002459986.1 | 18S rRNA |

> Zaten katalogda olan CCHFV, Batı Nil, Dengue, Chikungunya, Zika, sıtma (*P. falciparum*),
> *Leishmania*, *Coxiella burnetii* ve tularemi (*F. tularensis*) yeniden eklenmedi.

## Son doğrulama (05.11.2026)

| Test | Sonuç |
|---|---|
| `py_compile` (tüm modüller) | Geçti |
| AppTest — boş uygulama | 0 istisna, 6 üst sekme + alt sekmeler |
| AppTest — CCHFV örnek assay | 0 istisna, aktif assay 145 bp |
| AppTest — hedef in-silico PCR | 0 istisna, 145 bp tam eşleşme |
| AppTest — LOD/LOQ | 0 istisna, LOD95 ≈ 3,09 kopya, LOQ 25 kopya |
| AppTest — genotip kapsama | 10/10 dizi, tam amplikon |
| AppTest — literatür karşılaştırma | 0 istisna, uyum skoru 100/100 |
| AppTest — SNP | 8 aday, 3' uç ayrımı |
| AppTest — RPA | 6 aday |
| AppTest — LAMP | 5 set |
| AppTest — consensus (gevşetilmiş) | 10 aday, kapsama uyarısı gösterildi |
| AppTest — TBEV tek tasarım | 0 istisna, 154 bp, prob üretildi |
| AppTest — *B. burgdorferi* büyük genom | 0 istisna, 163 bp, ref 7.012 bp'ye kırpıldı |
| PDF assay kartı | 45.762 bayt, geçerli `%PDF` |
| 10 yeni etken bağımsız tasarım testi | 10/10 katı kriter, prob dimer yok |

**Doğrulanmış CCHFV seti** (NC_005302.1, S segmenti, 1101–1245):

| Oligo | Dizi (5'→3') | Uzunluk | Tm | GC |
|---|---|---|---|---|
| Forward | TGAGCACCCCAATGAAGTG | 19 nt | 61,8 °C | %52,6 |
| Reverse | ACAGGGATTGTCCCAAAGC | 19 nt | 61,8 °C | %52,6 |
| Prob | ACATGCATCCTGCCGTGCTTACAG | 24 nt | 69,2 °C | %54,2 |

Amplikon 145 bp · prob amplikonun 76. bazında · prob–primer dimer yok · diğer 8 virüse karşı çapraz reaksiyon yok.

> **Türkiye suşu uyarısı (05.11.2026):** Bu örnek set, NCBI'daki 20 CCHFV Türkiye S segmenti
> dizisiyle karşılaştırıldığında 20 dizinin 18'inde kritik mismatch gösterdi (özellikle prob
> bölgesinde 2–4 mismatch; bazı dizilerde primer 3' ucunda mismatch). **Bu seti Türkiye
> suşları için olduğu gibi sipariş etmeyin.** Sipariş Kontrolü sekmesindeki Türkiye kontrolünü
> çalıştırın; yeni/dejenere tasarım veya farklı hedef bölge gerekir. Bu sonuç örneklem
> dizilere dayalı in-silico ön kontroldür, tüm Türkiye izolatlarını temsil etmez.

## Öncelikli 5 etken — Türkiye taraması (06.11.2026)

NCBI `Turkey[Country]` araması ve hedef bölgeye lokal homoloji filtresiyle yapılan kontrol.
Hedef bölgeyi içermeyen kısa/partial kayıtlar **N/A** sayılır; PASS/WARN/FAIL ve sipariş
kapısı yalnızca değerlendirilebilir dizileri kullanır.

| Etken | Türkiye kaydı | Mevcut katalog seti | Durum |
|---|---|---|---|
| CCHFV | 22 S segmenti | 0 PASS · 2 WARN · 20 FAIL | **Sipariş engelli.** Türkiye referanslı en iyi aday: 1 PASS · 20 WARN · 1 FAIL; hâlâ kapıyı geçmiyor. |
| Batı Nil | 2 tam genom + 1 kısmi hedef bölge | Lineage-2 set: 0 PASS · 0 WARN · 3 FAIL | **Sipariş engelli.** OP617270.1 NS5 adayı 2/2 tam genomda PASS; örneklem yalnızca 2, sipariş için yeterli değil. |
| Dengue | 0 | Değerlendirilemedi | NCBI'da Türkiye kaydı yok; yerel/importe izolat FASTA gerekli. |
| Zika | 0 | Değerlendirilemedi | NCBI'da Türkiye kaydı yok; yerel/importe izolat FASTA gerekli. |
| Chikungunya | 0 | Değerlendirilemedi | NCBI'da Türkiye kaydı yok; yerel/importe izolat FASTA gerekli. |

> Bu tablo sipariş listesi değildir. CCHFV ve Batı Nil için bulunan adaylar in-silico
> ön değerlendirme aşamasındadır; Primer-BLAST, dimer analizi, LOD/LOQ ve ıslak laboratuvar
> doğrulaması olmadan kullanılmamalıdır.

## Baştan sona referans bağlanma kontrolü (07.11.2026)

Kısa oligo dizilerinin Google'da bulunmaması, yanlış tasarım sonucu değildir. Uygulama hazır
literatür dizileri sunmaz; referans diziden yeni aday üretir. 19–24 nt'lik oligo dizileri arama
motorlarında düzenli indekslenmez. Ayrıca reverse primerin sipariş dizisi (`5'→3'`) GenBank/NCBI
ileri iplikte ters komplement olarak bulunur; ileri iplikte aranan dizi ters komplementtir.

| Kontrol | Sonuç |
|---|---|
| Kod derleme | `PY_COMPILE_OK` |
| AppTest — boş uygulama | 0 istisna |
| CCHFV referans eşleşmesi | `NC_005302.1`: F 1101–1119, R bağlanma bölgesi 1227–1245, prob 1176–1199; amplikon 1101–1245, 145 bp, 0 mismatch |
| 5 öncelikli etken referans eşleşmesi | F/Prob ileri iplikte tam eşleşme; reverse primer ters komplementi tam eşleşme |
| Arayüz | Referans bağlanma haritası + Google/Primer-BLAST doğrulama notu eklendi |
| Yayınlanmış Assay'ler sekmesi | Semantic Scholar/Crossref arama 25 sonuç döndürdü; atıf sayısına göre sıralama çalıştı; AppTest 0 istisna |

## Kalan / isteğe bağlı işler

1. **Islak laboratuvar doğrulaması** — tüm setler için Primer-BLAST, dimer analizi, LOD/LOQ
2. (İsteğe bağlı) **Şifre koruması** — `st.secrets` içinde şifre, giriş ekranı
3. (İsteğe bağlı) Yeni etkenlerin kalıcı olması — `pathogens.py` içindeki `PATHOGENS` sözlüğüne eklemek
   (kenar çubuğundan eklenenler yalnızca oturum boyunca kalır)
4. (İsteğe bağlı) BLAST URL API için kurumsal anahtar/limit ayarı

## Nasıl devam edilir

- **Yerel çalıştırma:** `BASLAT.bat` → `http://localhost:8501`
- **Ağdan erişim:** `BASLAT_AG.bat`
- **Kod değişikliğinden sonra GitHub'a gönder:**
  ```powershell
  git add -A
  git commit -m "aciklama"
  git push
  ```
  Streamlit Cloud otomatik yeniden yayınlar.

## Önemli notlar

- `Entrez.email` sırası: `st.secrets["NCBI_EMAIL"]` → `NCBI_EMAIL` ortam değişkeni → varsayılan
- Virüs korunmuşluğu: 8 dizi ≈ 2–5 sn/etken; bakteri genomları büyük olduğu için ≈ 10 sn
- Tek etken özgüllük taraması 8 virüs genomunu kapsar; bakteri/parazitler için panel sekmesi kullanılmalı
- Prob–primer dimer eşiği ≥6 baz; primer–primer dimer eşiği ≥5 baz
- Genotip/consensus için NCBI'da homolog dizi bulunması gerekir; bulunamazsa yalnızca referans gösterilir
- LOD/LOQ **simülasyondur**; ıslak laboratuvar validasyonunun yerine geçmez
- RPA/LAMP **beta**dır; loop primerleri ve reaksiyon optimizasyonu ayrıca yapılmalıdır
- PDF ve kütüphane çıktıları araştırma amaçlıdır, tıbbi karar desteği değildir
- Streamlit Community Cloud ücretsiz; uygulama linki herkese açıktır (kod private olsa bile)
