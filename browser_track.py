"""
browser_track.py
Hafif, Plotly tabanli genom izleyici:
  - gen/ozellik oklari
  - primer/prob konumlari
  - korunmusluk izi
  - mismatch yogunlugu
"""

from __future__ import annotations


def navy_template() -> str:
    """Plotly grafiklerini koyu lacivert temayla uyumlu hale getirir."""
    try:
        import plotly.io as pio

        if "navy" not in pio.templates:
            tpl = pio.templates["plotly_dark"]
            tpl.layout.paper_bgcolor = "#0B1526"
            tpl.layout.plot_bgcolor = "#0B1526"
            tpl.layout.font.color = "#E6EDF7"
            tpl.layout.title.font.color = "#E6EDF7"
            tpl.layout.legend.font.color = "#E6EDF7"
            tpl.layout.xaxis.gridcolor = "#2E425E"
            tpl.layout.yaxis.gridcolor = "#2E425E"
            tpl.layout.xaxis.tickfont.color = "#E6EDF7"
            tpl.layout.yaxis.tickfont.color = "#E6EDF7"
            tpl.layout.xaxis.title.font.color = "#E6EDF7"
            tpl.layout.yaxis.title.font.color = "#E6EDF7"
            pio.templates["navy"] = tpl
    except Exception:
        return "plotly_dark"
    return "navy"


def build_genome_track(region_start: int, region_end: int, features=None, primers=None,
                       conservation=None, mismatch_density=None, title: str = "",
                       ref_len: int | None = None):
    """
    region_start/end: 0-based yari acik.
    features: [{'name','start','end','strand','color'}]  (0-based yari acik)
    primers:  [{'name','start','end','color'}]           (0-based yari acik)
    conservation: 0..1 listesi (region uzunlugunda)
    mismatch_density: sayi listesi (region uzunlugunda)
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    region_start = int(region_start)
    region_end = int(region_end)
    if region_end <= region_start:
        region_end = region_start + 1
    x = list(range(region_start + 1, region_end + 1))
    fig = make_subplots(
        rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.07,
        row_heights=[0.46, 0.30, 0.24],
        subplot_titles=("Gen / özellik ve primer-prob konumları", "Korunmuşluk", "Mismatch yoğunluğu"),
    )

    # Hedef bolge zemini
    fig.add_vrect(x0=region_start + 0.5, x1=region_end + 0.5, fillcolor="#e3f2fd",
                  opacity=0.35, line_width=0, row=1, col=1)

    # Genler
    features = features or []
    if features:
        for i, f in enumerate(features):
            x0 = max(region_start, int(f.get("start", region_start))) + 1
            x1 = min(region_end, int(f.get("end", region_end)))
            if x1 <= x0:
                continue
            color = f.get("color", "#90a4ae")
            y0, y1 = 0.18, 0.52
            fig.add_shape(type="rect", x0=x0, x1=x1, y0=y0, y1=y1,
                          fillcolor=color, opacity=0.45, line=dict(color=color, width=1),
                          row=1, col=1)
            fig.add_annotation(x=(x0 + x1) / 2, y=(y0 + y1) / 2, text=str(f.get("name", ""))[:30],
                               showarrow=False, font=dict(size=9, color="#E6EDF7"), row=1, col=1)
    else:
        fig.add_shape(type="rect", x0=region_start + 1, x1=region_end, y0=0.18, y1=0.52,
                      fillcolor="#cfd8dc", opacity=0.4, line=dict(color="#90a4ae", width=1),
                      row=1, col=1)
        fig.add_annotation(x=(region_start + region_end) / 2, y=0.35, text="Referans bölge",
                           showarrow=False, font=dict(size=9, color="#D6DEEE"), row=1, col=1)

    # Primer / prob
    primers = primers or []
    for p in primers:
        x0 = max(region_start, int(p.get("start", region_start))) + 1
        x1 = min(region_end, int(p.get("end", region_end)))
        if x1 <= x0:
            continue
        color = p.get("color", "#1565c0")
        fig.add_shape(type="rect", x0=x0, x1=x1, y0=0.62, y1=0.92,
                      fillcolor=color, opacity=0.75, line=dict(color=color, width=1),
                      row=1, col=1)
        fig.add_annotation(x=(x0 + x1) / 2, y=0.95, text=str(p.get("name", "")),
                           showarrow=False, font=dict(size=8, color=color), row=1, col=1)

    # Korunmusluk
    if conservation:
        vals = [None if v is None else 100 * float(v) for v in conservation]
        fig.add_trace(go.Scatter(
            x=x[:len(vals)], y=vals, mode="lines", name="Korunmuşluk (%)",
            line=dict(color="#1b5e20", width=2), fill="tozeroy", fillcolor="rgba(27,94,32,0.10)",
            hovertemplate="%{x} bp<br>%{y:.1f}%<extra></extra>",
        ), row=2, col=1)
        fig.update_yaxes(range=[0, 105], title_text="%", row=2, col=1)
    else:
        fig.add_annotation(x=(region_start + region_end) / 2, y=50,
                           text="Korunmuşluk verisi yok (çoklu hizalama çalıştırın)",
                           showarrow=False, font=dict(size=9, color="#B0BEC5"), row=2, col=1)
        fig.update_yaxes(range=[0, 100], title_text="%", row=2, col=1)

    # Mismatch yogunlugu
    if mismatch_density:
        vals = list(mismatch_density)[:len(x)]
        fig.add_trace(go.Bar(
            x=x[:len(vals)], y=vals, name="Mismatch", marker_color="#c62828", opacity=0.75,
            hovertemplate="%{x} bp<br>%{y} mismatch<extra></extra>",
        ), row=3, col=1)
        fig.update_yaxes(title_text="sayı", row=3, col=1)
    else:
        fig.add_annotation(x=(region_start + region_end) / 2, y=1,
                           text="Mismatch verisi yok (genotip analizi çalıştırın)",
                           showarrow=False, font=dict(size=9, color="#B0BEC5"), row=3, col=1)
        fig.update_yaxes(range=[0, 2], title_text="sayı", row=3, col=1)

    fig.update_xaxes(title_text="Genomik pozisyon (bp)", row=3, col=1)
    fig.update_yaxes(showticklabels=False, range=[0, 1], row=1, col=1)
    fig.update_layout(
        height=560, margin=dict(l=20, r=20, t=45, b=40),
        template=navy_template(), showlegend=False,
        title=dict(text=title, font=dict(size=13)),
        hovermode="x unified",
    )
    return fig
