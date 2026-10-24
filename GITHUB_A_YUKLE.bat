@echo off
setlocal
cd /d "%~dp0"
title GitHub'a Yukle - PCR Primer Designer

echo ============================================================
echo   PROJEYI GITHUB'A YUKLE
echo ============================================================
echo.
echo ONCE su adreste BOS bir depo olusturun:
echo    https://github.com/new
echo.
echo    Repository name : PCR_PRIMER_PROB
echo    Private secenegi olabilir
echo    "Add a README file" KUTUSUNU ISARETLEMEYIN
echo.
echo Depoyu olusturduktan sonra bu pencereye donun.
echo.
pause

set /p GHUSER=GitHub kullanici adiniz: 
if "%GHUSER%"=="" (
    echo.
    echo HATA: Kullanici adi bos olamaz.
    pause
    exit /b 1
)

echo.
echo Baglanti kuruluyor: https://github.com/%GHUSER%/PCR_PRIMER_PROB.git
git remote remove origin >nul 2>&1
git remote add origin https://github.com/%GHUSER%/PCR_PRIMER_PROB.git
git branch -M main

echo.
echo Dosyalar gonderiliyor... (tarayici acilip GitHub girisi isterse izin verin)
echo.
git push -u origin main

if errorlevel 1 (
    echo.
    echo ============================================================
    echo  YUKLEME BASARISIZ
    echo ============================================================
    echo  Kontrol listesi:
    echo   1. Depo gercekten olusturuldu mu? (github.com/%GHUSER%/PCR_PRIMER_PROB)
    echo   2. Depoda README dosyasi var mi? Varsa once onu silin.
    echo   3. GitHub sifresi yerine TOKEN isteniyorsa:
    echo      https://github.com/settings/tokens - Generate new token (classic)
    echo      "repo" yetkisini isaretleyip token'i sifre yerine yapistirin.
    echo.
) else (
    echo.
    echo ============================================================
    echo  BASARILI - Kod GitHub'da
    echo ============================================================
    echo  Simdi sirada buluta yayinlamak var:
    echo    https://share.streamlit.io
    echo.
    echo  New app - Repository: %GHUSER%/PCR_PRIMER_PROB
    echo  Main file: realtime_pcr_designer.py
    echo  Advanced settings - Secrets:
    echo      NCBI_EMAIL = "ad.soyad@kurum.gov.tr"
    echo.
    start "" https://share.streamlit.io
)

echo.
pause
