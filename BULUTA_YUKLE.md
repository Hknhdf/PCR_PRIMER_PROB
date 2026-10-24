# Ev bilgisayarı KAPALIYKEN iş yerinden erişim — adım adım

Ev bilgisayarı kapandığında oradaki sunucu da kapanır. Bu yüzden uygulamayı
**ücretsiz bir bulut sunucusuna** koyacağız: Streamlit Community Cloud.

Sonuçta iş yerinde sadece tarayıcıya `https://xxxx.streamlit.app` yazman yeterli olur.
Python kurmak, dosya taşımak, ev bilgisayarını açık bırakmak gerekmez.

---

## Adım 1 — GitHub hesabı

https://github.com/signup adresinden ücretsiz hesap aç (varsa atla).

## Adım 2 — Boş depo (repository) oluştur

1. https://github.com/new
2. **Repository name:** `PCR_PRIMER_PROB`
3. **Private** seçebilirsin (Streamlit Cloud özel depoyu da destekler)
4. **"Add a README file" işaretleme** (klasörde zaten var)
5. **Create repository**

## Adım 3 — Bu klasörü GitHub'a gönder

### Kolay yol (önerilen)

`GITHUB_A_YUKLE.bat` dosyasına **çift tıkla** → GitHub kullanıcı adını yaz → gerisini yapar.

Tarayıcı açılıp GitHub girişi isterse **izin ver**.

### Elle yol

PowerShell aç ve şu 3 satırı çalıştır (`KULLANICI_ADIN` yerine kendi adın):

```powershell
git remote add origin https://github.com/KULLANICI_ADIN/PCR_PRIMER_PROB.git
git branch -M main
git push -u origin main
```

> Not: Bu klasörde git deposu ve commit'ler zaten hazır. Sadece yukarıdaki komutlar gerekli.
> Depoyu oluştururken **"Add a README file" işaretlenmemeli**, yoksa push reddedilir.

## Adım 4 — Streamlit Cloud'da yayınla

1. https://share.streamlit.io → **"Continue with GitHub"**
2. **"New app"** → **"Deploy a public app from GitHub"** (özel depo için GitHub erişim izni istenir)
3. Doldur:
   - **Repository:** `KULLANICI_ADIN/PCR_PRIMER_PROB`
   - **Branch:** `main`
   - **Main file path:** `realtime_pcr_designer.py`
4. **"Advanced settings"** → **Secrets** kutusuna şunu yaz:

```toml
NCBI_EMAIL = "ad.soyad@kurum.gov.tr"
```

5. **Deploy** — 2-3 dakika sonra adres hazır:

```
https://pcr-primer-prob-xxxx.streamlit.app
```

Bu adresi iş yerinde tarayıcıya yaz. Hepsi bu.

---

## Bilmen gerekenler

| Konu | Açıklama |
|---|---|
| Ev bilgisayarı | Kapalı olabilir, hiç gerekmiyor |
| İş yerinde kurulum | Yok — sadece tarayıcı |
| NCBI erişimi | Buluttan yapılır; iş yeri ağı NCBI'yi engellese bile çalışır |
| Uyku modu | Birkaç saat kullanılmazsa uyur; ilk açılışta 20-40 sn uyanır (normal) |
| Gizlilik | Depo private olabilir; uygulama sadece NCBI'nin halka açık dizilerini kullanır |
| Güncelleme | Bilgisayarda değişiklik yapıp `git add -A; git commit -m "guncelleme"; git push` → bulut otomatik yenilenir |
| Excel raporları | Her kullanımda indirilir, sunucuda saklanmaz |

## Alternatif: iş yerindeki bilgisayarda yerel çalıştırma

Bulut istemezsen: `PCR_PRIMER_PROB_tasinabilir.zip` dosyasını iş yerine götür,
Python kur, `KURULUM.bat` → `BASLAT.bat`. Avantajı internet dışına kapalı çalışır;
dezavantajı iş yeri bilgisayarında Python kurma izni gerekmesi.

## Sorun giderme

- **`git push` kimlik doğrulama hatası:** GitHub şifresi yerine *Personal Access Token* ister.
  https://github.com/settings/tokens → Generate new token (classic) → `repo` yetkisi → token'ı şifre yerine yapıştır.
- **Streamlit Cloud "app is sleeping":** sayfayı yenile, kendiliğinden uyanır.
- **Uygulama hata verdi:** share.streamlit.io → uygulama → **Logs** sekmesinde hata mesajı görünür.
- **NCBI e-posta uyarısı:** Secrets alanındaki `NCBI_EMAIL` satırını kontrol et.
