@echo off
cd /d "%~dp0"
echo Real-Time PCR Primer & Probe Designer baslatiliyor...
echo.
echo Tarayicida http://localhost:8501 adresini acin.
echo Bu pencereyi KAPATMAYIN (kapatirsaniz uygulama durur).
echo.
start "" http://localhost:8501
python -m streamlit run realtime_pcr_designer.py --server.port 8501 --server.address localhost
pause
