# PROJE DURUMU — kaldığım yerden devam

Son güncelleme: 25.10.2026

## Proje nedir

Virüs, bakteri ve parazitler için NCBI'dan canlı dizi çekip TaqMan primer/prob tasarlayan
Streamlit uygulaması. Hazır katalog: **9 virüs, 21 bakteri, 10 parazit**.

## Şu an nerede

| Konu | Durum |
|---|---|
| Kod | `C:\Users\Hakan Hedef\Desktop\PCR_PRIMER_PROB` |
| GitHub | `Hknhdf/PCR_PRIMER_PROB` — `main` dalı |
| Streamlit Cloud | Yayında — `https://pcr-primer-prob.streamlit.app` |
| Yerel uygulama | `BASLAT.bat` → `http://localhost:8501` |

## Tamamlanan özellikler

- Canlı NCBI erişimi (Biopython Entrez), önbellekli
- SantaLucia nearest-neighbour Tm (qPCR tuz koşulları)
- Gerçek amplikon mantığı (forward = sense başlangıç, reverse = amplikon sonunun ters komplementi)
- Katı / gevşek iki kademeli filtre, dimer ve homopolymer kontrolü
- **Korunmuşluk analizi:** k-mer ön filtre → 24-mer çapa → global/lokal hizalama → pozisyon bazlı oran
- **Multiplex panel:** çoklu etken, primerler arası çapraz reaksiyon taraması (tam / 3' uç / 5' uç)
- **Excel raporu:** Özet, Primerler, Problar, Korunmuşluk, Çapraz Reaksiyon, Protokol, Adaylar
- Manuel FASTA/ham dizi modu (internet yoksa)
- **Etken kataloğu:** virüs / bakteri / parazit, kategoriye göre gruplu seçim
  - Bakteri ve parazitlerde tüm genom yerine **hedef gen** taranır (ör. `nuc`, `IS6110`, 18S rRNA)
  - Büyük genomlarda korunmuşluk için dizi sayısı otomatik sınırlandırılır
  - AT-zengin bölgelerde prob Tm/GC hedefleri otomatik gevşetilir (uyarı gösterilir)
  - T. gondii hedefi 529 bp tekrar elementi (AF146527.1) — 6 diziyle %96,6 korunmuşluk
- **Özgüllük filtresi:** tek etken tasarımında diğer katalog virüslerine bağlanabilecek adaylar elenir (varsayılan açık)
- **Prob–primer dimer filtresi:** prob, primerlerle ≥6 baz tamamlayıcılık yapmayacak şekilde seçilir
- Varsayılan aday sayısı 10 (3–15 arası ayarlanabilir)
- Kenar çubuğundan yeni etken ekleme (oturumluk) + kataloğu JSON olarak indirme
- **Hazır etken şablonu:** “Yeni Etken / Varyant Ekle” formunda listeden etken seçince kategori, GenBank ID, bölge ve gen alanları otomatik dolar
- Taşınabilir paket: `PCR_PRIMER_PROB_tasinabilir.zip`

## Son doğrulama (25.10.2026)

Bağımsız doğrulama betiği (dizileri referansa geri eşleme, Tm/GC yeniden hesabı, dimer,
korunmuşluk, çapraz reaksiyon, Excel) ve arayüz testleri:

| Test | Sonuç |
|---|---|
| Bağımsız kontrol | **96 PASS · 0 FAIL · 0 WARN** |
| 40 etken (primer + prob) | **40/40** · dimer uyarısı 0 · 1 gevşek prob (L. donovani, AT-zengin) |
| Arayüz: CCHFV korunmuşluklu tasarım | Hatasız |
| Arayüz: 3'lü karışık panel (CCHFV + WNV + S. aureus) | Hatasız, çapraz eşleşme yok |
| Arayüz: manuel FASTA modu (800 bp) | Hatasız, amplikon 151 bp |
| Excel raporu | 7 sayfa, içerik doğrulandı |

**Doğrulanmış CCHFV seti** (NC_005302.1, S segmenti, 1101–1245):

| Oligo | Dizi (5'→3') | Uzunluk | Tm | GC | Korunmuşluk |
|---|---|---|---|---|---|
| Forward | TGAGCACCCCAATGAAGTG | 19 nt | 61,8 °C | %52,6 | %98,7 |
| Reverse | ACAGGGATTGTCCCAAAGC | 19 nt | 61,8 °C | %52,6 | %94,7 |
| Prob | ACATGCATCCTGCCGTGCTTACAG | 24 nt | 69,2 °C | %54,2 | %92,7 |

Amplikon 145 bp · prob amplikonun 76. bazında · prob–primer dimer yok · diğer 8 virüse karşı çapraz reaksiyon yok.

## Kalan / isteğe bağlı işler

1. (İsteğe bağlı) **Şifre koruması** — `st.secrets` içinde şifre, giriş ekranı
2. (İsteğe bağlı) Yeni etkenlerin kalıcı olması — `pathogens.py` içindeki `PATHOGENS` sözlüğüne eklemek
   (kenar çubuğundan eklenenler yalnızca oturum boyunca kalır)

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
- Virüs korunmuşluğu: 8 dizi ≈ 2-5 sn/etken; bakteri genomları büyük olduğu için ≈ 10 sn
- Tek etken özgüllük taraması 8 virüs genomunu kapsar; bakteri/parazitler için panel sekmesi kullanılmalı
- Prob–primer dimer eşiği ≥6 baz; primer–primer dimer eşiği ≥5 baz
- Excel çapraz reaksiyon taraması k-mer tabanlıdır, tam hizalama değildir
- Üretilen setler klinik kullanım öncesi **Primer-BLAST + dimer analizi + LOD/LOQ** ile doğrulanmalı
- Streamlit Community Cloud ücretsiz; uygulama linki herkese açıktır (kod private olsa bile)
