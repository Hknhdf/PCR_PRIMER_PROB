@echo off
cd /d "%~dp0"
echo ============================================
echo  PCR Tasarim Motoru - Ag Erisimi
echo ============================================
echo.
echo Bu pencereyi KAPATMAYIN.
echo.

for /f "tokens=2 delims=:" %%a in ('ipconfig ^| findstr /i "IPv4"') do set IP=%%a
set IP=%IP: =%

echo Bu bilgisayardan:        http://localhost:8501
if defined IP echo Ayni agdaki cihazlardan:  http://%IP%:8501
echo.
echo Not: Diger cihazdan erisemiyorsaniz Windows Guvenlik Duvari'nda
echo      8501 portuna izin vermeniz gerekir.
echo.

start "" http://localhost:8501
python -m streamlit run realtime_pcr_designer.py --server.port 8501 --server.address 0.0.0.0
pause
