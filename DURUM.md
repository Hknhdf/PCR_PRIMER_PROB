# PROJE DURUMU — kaldığım yerden devam

Son güncelleme: 24.10.2026

## Proje nedir

CCHFV, Batı Nil, Hantavirüs, Dengue, Chikungunya, Zika, Ebola ve Marburg için
NCBI'dan canlı dizi çekip TaqMan primer/prob tasarlayan Streamlit uygulaması.

## Şu an nerede

| Konu | Durum |
|---|---|
| Kod | `C:\Users\Hakan Hedef\Desktop\PCR_PRIMER_PROB` (tüm dosyalar hazır) |
| Git | `main` dalı, commit `47a2241`, çalışma ağacı temiz |
| GitHub | **Yüklendi** — `Hknhdf/PCR_PRIMER_PROB` (private), uzak `main` = yerel `main` |
| Streamlit Cloud | **Yarım** — deploy formu açıldı, Deploy'a basılıp basılmadığı belirsiz |
| Yerel uygulama | Çalışıyordu, `http://localhost:8501` (PC kapanınca durur, sorun değil) |

## Tamamlanan özellikler

- Canlı NCBI erişimi (Biopython Entrez), önbellekli
- SantaLucia nearest-neighbour Tm (qPCR tuz koşulları)
- Gerçek amplikon mantığı (forward = sense başlangıç, reverse = amplikon sonunun ters komplementi)
- Katı / gevşek iki kademeli filtre, dimer ve homopolymer kontrolü
- **Korunmuşluk analizi:** k-mer ön filtre → 24-mer çapa → global/lokal hizalama → pozisyon bazlı oran
  (CCHFV testinde 8 dizi, ortalama %95,5; primer/prob seçimi buna göre ağırlıklandırılıyor)
- **Multiplex panel:** çoklu etken, primerler arası çapraz reaksiyon taraması (tam / 3' uç / 5' uç)
- **Excel raporu:** Özet, Primerler, Problar, Korunmuşluk, Çapraz Reaksiyon, Protokol, Adaylar
- Manuel FASTA/ham dizi modu (internet yoksa)
- Taşınabilir paket: `PCR_PRIMER_PROB_tasinabilir.zip`

## Kalan işler

1. **Streamlit Cloud deploy'unu bitir**
   - https://share.streamlit.io → New app
   - Repository: `Hknhdf/PCR_PRIMER_PROB` · Branch: `main` · Main file: `realtime_pcr_designer.py`
   - Advanced settings → Secrets: `NCBI_EMAIL = "..."` (kendi e-postası)
   - Private depo listede çıkmazsa: **Configure GitHub access** → depoyu seç → authorize
2. (İsteğe bağlı) **Şifre koruması** — `st.secrets` içinde şifre, giriş ekranı
3. (İsteğe bağlı) Yeni virüs eklemelerinin kalıcı olması — `DEFAULT_VIRUS_DATABASE` sözlüğüne eklemek
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
- Korunmuşluk analizi NCBI'dan ek istek atar; 8 dizi ≈ 2-5 sn/etken (sonrası önbellek)
- Excel çapraz reaksiyon taraması k-mer tabanlıdır, tam hizalama değildir
- Üretilen setler klinik kullanım öncesi **Primer-BLAST + dimer analizi + LOD/LOQ** ile doğrulanmalı
- Streamlit Community Cloud ücretsiz; uygulama linki herkese açıktır (kod private olsa bile)
