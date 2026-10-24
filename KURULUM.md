# PCR Primer & Prob Tasarım Motoru — Kurulum ve Taşıma Rehberi

Bu klasörün tamamı taşınabilir. Python kurulu her Windows bilgisayarda çalışır.

---

## 1) Aynı bilgisayar, başka bir gün

En kolay senaryo. Dosyalar silinmediyse:

1. `BASLAT.bat` dosyasına **çift tıkla**
2. Tarayıcıda `http://localhost:8501` açılır
3. Açılan siyah pencereyi **kapatma** (kapatırsan uygulama durur)

Kapatmak için: pencereyi kapat veya `Ctrl + C`.

---

## 2) Başka bir bilgisayara taşımak

### Adım 1 — Dosyaları kopyala
`PCR_PRIMER_PROB_tasinabilir.zip` dosyasını USB / e-posta / OneDrive ile diğer bilgisayara taşı ve bir klasöre çıkart.
(Örn: `C:\PCR_PRIMER_PROB`)

> Zip içinde: 4 Python dosyası, 4 `.bat` dosyası, `requirements.txt`, `.streamlit` klasörü, bu rehber.

### Adım 2 — Python kur (yoksa)
https://www.python.org/downloads/ → Python 3.10 veya üzeri.
Kurulum ekranında **"Add Python to PATH"** kutusunu mutlaka işaretle.

### Adım 3 — Paketleri kur (tek seferlik)
`KURULUM.bat` dosyasına çift tıkla. Şunları kurar:

```
streamlit, biopython, pandas, openpyxl
```

### Adım 4 — Başlat
`BASLAT.bat` → `http://localhost:8501`

### Adım 5 — NCBI e-postası
`realtime_pcr_designer.py` dosyasını Not Defteri ile aç, şu satırı kendi e-postanla değiştir:

```python
Entrez.email = os.environ.get("NCBI_EMAIL", "hakan.hedef@example.com")
```

veya Windows'ta ortam değişkeni tanımla:

```powershell
setx NCBI_EMAIL "ad.soyad@kurum.gov.tr"
```

NCBI, ücretsiz API kullanımında geçerli e-posta ister (nezaket kuralı).

---

## 3) Aynı ağdaki başka cihazdan girmek (telefon / tablet / ikinci PC)

Sunucu bu bilgisayarda çalışır, diğer cihazlar tarayıcıdan bağlanır.

1. `BASLAT_AG.bat` dosyasına çift tıkla
2. Pencerede yazan adresi diğer cihazın tarayıcısına yaz:
   `http://192.168.x.x:8501`
3. **Erişemiyorsan** Windows Güvenlik Duvarı'nda izin ver:

```powershell
# Yönetici olarak PowerShell aç
netsh advfirewall firewall add rule name="PCR PCR 8501" dir=in action=allow protocol=TCP localport=8501
```

Not: İki cihaz aynı Wi-Fi / ağda olmalı. Kurumsal ağlarda port kapalı olabilir.

---

## 4) İnternetten erişim (her yerden)

İki seçenek var:

### A. Streamlit Community Cloud (ücretsiz, kalıcı)
1. Projeyi GitHub'a yükle (public repo)
2. https://share.streamlit.io → "New app" → repo + `realtime_pcr_designer.py`
3. `requirements.txt` otomatik okunur
4. Uygulama `https://xxx.streamlit.app` adresinde yayınlanır

> Dikkat: Repo herkese açık olur. Hasta verisi, kurum içi dizi veya ticari bilgi koymayın.
> Bu uygulama yalnızca NCBI'daki **halka açık referans dizileri** kullanır.

### B. Geçici tünel (kendi bilgisayarından, hızlı test)
1. `BASLAT.bat` ile uygulamayı çalıştır
2. Ayrı bir terminalde:

```powershell
winget install --id Cloudflare.cloudflared
cloudflared tunnel --url http://localhost:8501
```

3. Ekranda çıkan `https://xxx.trycloudflare.com` adresini paylaş
4. Terminal kapanınca tünel de kapanır

Not: Tünel açıkken uygulaman bilgisayarın üzerinde çalışmaya devam eder.

---

## 5) Günlük kullanım ipuçları

| İhtiyaç | Yol |
|---|---|
| Yeni virüs eklemek | Kenar çubuğu → "Yeni Virüs Ekle" (sadece o oturum için; kalıcı olsun istersen `DEFAULT_VIRUS_DATABASE` içine ekle) |
| Sonuçları saklamak | Her tasarımdan sonra **Excel raporunu indir** — uygulama geçmişi tutmaz |
| İnternet yok | "Manuel dizi / FASTA" modu (korunmuşluk analizi organizma adı ister) |
| Uygulama açılmıyor | `netstat -ano | findstr :8501` → port doluysa eski süreci kapat |
| Yavaş çalışıyor | Kenar çubuğunda karşılaştırılan dizi sayısını düşür (8 → 5) |

---

## 6) Dosya sorumlulukları

```
realtime_pcr_designer.py   Arayüz (Streamlit, 3 sekme)
pcr_engine.py              Tm, dimer, primer/prob tasarımı
conservation.py            k-mer filtre + hizalama + korunmuşluk
report.py                  7 sayfalık Excel raporu
requirements.txt           Python bağımlılıkları
.streamlit/config.toml     Sunucu/tema ayarları
BASLAT.bat                 Yerel başlatma
BASLAT_AG.bat              Ağdan erişim için başlatma
KURULUM.bat                Paket kurulumu
```

---

## 7) Önemli uyarı

Bu uygulama **in silico ön değerlendirme** aracıdır. Ürettiği primer/prob setleri
laboratuvar ve klinik kullanım öncesi şu doğrulamalardan geçmelidir:

- Primer-BLAST / BLASTn ile özgüllük
- Dimer / hairpin analizi (OligoAnalyzer vb.)
- LOD / LOQ ve klinik örnek validasyonu
