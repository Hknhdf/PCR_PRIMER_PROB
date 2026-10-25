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
- Kenar çubuğundan yeni etken ekleme (oturumluk) + kataloğu JSON olarak indirme
- Taşınabilir paket: `PCR_PRIMER_PROB_tasinabilir.zip`

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
- Excel çapraz reaksiyon taraması k-mer tabanlıdır, tam hizalama değildir
- Üretilen setler klinik kullanım öncesi **Primer-BLAST + dimer analizi + LOD/LOQ** ile doğrulanmalı
- Streamlit Community Cloud ücretsiz; uygulama linki herkese açıktır (kod private olsa bile)
