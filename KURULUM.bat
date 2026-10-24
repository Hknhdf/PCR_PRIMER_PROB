@echo off
cd /d "%~dp0"
echo ============================================
echo  Gerekli Python paketleri kuruluyor
echo ============================================
echo.
python --version
if errorlevel 1 (
    echo.
    echo HATA: Python bulunamadi.
    echo https://www.python.org/downloads/ adresinden Python 3.10+ kurun.
    echo Kurulumda "Add Python to PATH" secenegini isaretleyin.
    pause
    exit /b 1
)
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
echo.
if errorlevel 1 (
    echo Paket kurulumunda sorun olustu. Hata mesajini kontrol edin.
) else (
    echo Kurulum tamamlandi. BASLAT.bat dosyasina cift tiklayabilirsiniz.
)
pause
