"""
pdf_report.py
Tek assay icin PDF assay karti (reportlab).
"""

from __future__ import annotations

import io
import os
from datetime import datetime
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak,
)


def _font_name():
    candidates = [
        ("DejaVuSans", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
        ("DejaVuSans", "C:/Windows/Fonts/DejaVuSans.ttf"),
        ("Arial", "C:/Windows/Fonts/arial.ttf"),
    ]
    for name, path in candidates:
        if os.path.exists(path):
            try:
                pdfmetrics.registerFont(TTFont(name, path))
                return name
            except Exception:
                continue
    return "Helvetica"


def _p(text, style):
    return Paragraph(escape(str(text)), style)


def build_assay_pdf(assay: dict, thermo_rows=None, coverage=None, specificity=None,
                    lod=None, notes: str = "", genotype=None):
    """
    assay: {
      pathogen, category, segment, accession, region, forward, reverse, probe,
      forward_start, forward_end, reverse_start, reverse_end, probe_start, probe_end,
      amplicon_size, forward_tm, reverse_tm, probe_tm, forward_gc, reverse_gc, probe_gc
    }
    """
    font = _font_name()
    buf = io.BytesIO()
    doc = SimpleDocTemplate(
        buf, pagesize=A4, leftMargin=16 * mm, rightMargin=16 * mm,
        topMargin=14 * mm, bottomMargin=14 * mm,
        title="qPCR Assay Kartı", author="PCR Primer & Probe Designer",
    )
    styles = getSampleStyleSheet()
    h1 = ParagraphStyle("h1x", parent=styles["Heading1"], fontName=font, fontSize=16,
                        textColor=colors.HexColor("#0d47a1"), spaceAfter=4)
    h2 = ParagraphStyle("h2x", parent=styles["Heading2"], fontName=font, fontSize=11,
                        textColor=colors.HexColor("#1565c0"), spaceBefore=8, spaceAfter=3)
    body = ParagraphStyle("bodyx", parent=styles["BodyText"], fontName=font, fontSize=8.5,
                          leading=11)
    small = ParagraphStyle("smallx", parent=styles["BodyText"], fontName=font, fontSize=7.2,
                           leading=9, textColor=colors.HexColor("#555555"))
    mono = ParagraphStyle("monox", parent=body, fontName="Courier", fontSize=8, leading=10)

    story = []
    story.append(_p("qPCR Primer / Prob Tasarım Kartı", h1))
    story.append(_p(f"Oluşturma: {datetime.now().strftime('%d.%m.%Y %H:%M')}", small))
    story.append(Spacer(1, 4 * mm))

    meta = [
        ["Etken", assay.get("pathogen", "—")],
        ["Kategori", assay.get("category", "—")],
        ["Hedef", assay.get("segment", "—")],
        ["Erişim no", assay.get("accession", "—")],
        ["Bölge", assay.get("region", "—")],
        ["Amplikon", f"{assay.get('amplicon_size', '—')} bp"],
    ]
    t = Table(meta, colWidths=[30 * mm, 145 * mm])
    t.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, -1), font),
        ("FONTSIZE", (0, 0), (-1, -1), 8.5),
        ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#e3f2fd")),
        ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#bbdefb")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(t)

    story.append(_p("Assay dizileri", h2))
    rows = [["Oligo", "Dizi (5'→3')", "Konum", "Uzunluk", "Tm (°C)", "GC (%)"]]
    for name, seq, start, end, tm, gc in (
        ("Forward", assay.get("forward"), assay.get("forward_start"), assay.get("forward_end"),
         assay.get("forward_tm"), assay.get("forward_gc")),
        ("Reverse", assay.get("reverse"), assay.get("reverse_start"), assay.get("reverse_end"),
         assay.get("reverse_tm"), assay.get("reverse_gc")),
        ("Prob", assay.get("probe"), assay.get("probe_start"), assay.get("probe_end"),
         assay.get("probe_tm"), assay.get("probe_gc")),
    ):
        if not seq:
            continue
        pos = f"{start}–{end}" if start and end else "—"
        rows.append([name, Paragraph(escape(str(seq)), mono), pos, str(len(seq)), str(tm or "—"), str(gc or "—")])
    t = Table(rows, colWidths=[18 * mm, 92 * mm, 24 * mm, 15 * mm, 18 * mm, 18 * mm], repeatRows=1)
    t.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, -1), font),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1565c0")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#cfd8dc")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ]))
    story.append(t)

    if thermo_rows:
        story.append(_p("Termodinamik", h2))
        rows = [["Oligo", "Tm p3", "Hairpin ΔG", "Homodimer ΔG", "3' uç ΔG", "Uyarı"]]
        for r in thermo_rows:
            rows.append([
                r.get("name", ""), r.get("tm", "—"), r.get("hairpin_dg", "—"),
                r.get("homodimer_dg", "—"), r.get("three_prime_dg", "—"),
                "; ".join(r.get("warnings", []))[:80] or "—",
            ])
        t = Table(rows, colWidths=[18 * mm, 22 * mm, 30 * mm, 34 * mm, 26 * mm, 60 * mm], repeatRows=1)
        t.setStyle(TableStyle([
            ("FONTNAME", (0, 0), (-1, -1), font), ("FONTSIZE", (0, 0), (-1, -1), 7.5),
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8f5e9")),
            ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#c8e6c9")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ]))
        story.append(t)

    if coverage:
        story.append(_p("Genotip kapsama", h2))
        rows = [["Oligo", "Değerlendirilebilir", "Mükemmel (%)", "%95 GA", "3' risk"]]
        for key, c in coverage.items():
            rows.append([
                key, f"{c.get('n_evaluable')}/{c.get('n_total')}",
                c.get("perfect_pct"), str(c.get("perfect_ci")), c.get("three_prime_risk"),
            ])
        t = Table(rows, colWidths=[18 * mm, 45 * mm, 32 * mm, 45 * mm, 25 * mm], repeatRows=1)
        t.setStyle(TableStyle([
            ("FONTNAME", (0, 0), (-1, -1), font), ("FONTSIZE", (0, 0), (-1, -1), 8),
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#fff8e1")),
            ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#ffe082")),
        ]))
        story.append(t)

    if specificity:
        story.append(_p("Özgüllük", h2))
        story.append(_p(
            f"Skor: {specificity.get('score', '—')}/100 · "
            f"Hedef dışı amplikon: {specificity.get('off_targets', '—')} · "
            f"Not: {specificity.get('notes', '')}", body))

    if lod:
        story.append(_p("LOD / LOQ (simülasyon)", h2))
        story.append(_p(
            f"LOD95: {lod.get('lod95', '—')} kopya "
            f"(%95 GA {lod.get('lod95_ci', '—')}) · LOQ: {lod.get('loq', '—')} kopya · "
            f"Standart eğri eğimi: {lod.get('slope', '—')} · "
            f"Verim: {lod.get('efficiency', '—')}", body))

    if genotype:
        story.append(_p("Genotip özeti", h2))
        story.append(_p(
            f"Amplikon bulunan: {genotype.get('n_found')}/{genotype.get('n_total')} "
            f"(%{genotype.get('amplicon_found_pct')}) · "
            f"Ortalama identite: {genotype.get('identity_mean')} · "
            f"Ortalama korunmuşluk: {genotype.get('conservation_mean')}", body))

    if notes:
        story.append(_p("Notlar", h2))
        story.append(_p(notes, body))

    story.append(Spacer(1, 5 * mm))
    story.append(_p(
        "Yöntem: Tm SantaLucia nearest-neighbour (Biopython); termodinamik primer3-py; "
        "özgüllük ≤3 mismatch toleranslı in-silico PCR; LOD/LOQ Poisson + probit Monte Carlo. "
        "Bu rapor araştırma amaçlıdır; klinik kullanım öncesi Primer-BLAST, dimer analizi, "
        "LOD/LOQ ve ıslak laboratuvar doğrulaması gereklidir.", small))
    doc.build(story)
    return buf.getvalue()
