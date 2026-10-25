"""One palette for the whole app: UI tokens (same values as web/src/theme.css), the data colours every Plotly figure
uses, and the soft tints RDKit uses to highlight atoms. Engines take their colours from here, never hard-coded."""
# ---- UI tokens (theme.css: --ink --muted --faint --line --canvas --accent)
INK = "#161A23"
MUTED = "#6A7180"
FAINT = "#9AA1AE"
LINE = "#E6E8ED"
SURFACE = "rgba(0,0,0,0)"   # figures sit on the page's own white
PAPER = "#FFFFFF"           # solid white where something must cover the plot (hover labels)
SURFACE_2 = "#F0F2F5"       # plot gridlines
ACCENT = "#2D5BD0"

# ---- data colours (6, chosen to sit together): series 1 is the accent blue
DATA_BLUE = ACCENT
DATA_AMBER = "#C27C1A"
DATA_TEAL = "#1F8A7E"
DATA_CORAL = "#D0583A"
DATA_VIOLET = "#7457C2"
DATA_GREY = "#A6ADBA"
SERIES = [DATA_BLUE, DATA_AMBER, DATA_TEAL, DATA_CORAL, DATA_VIOLET, DATA_GREY]

BLUE = DATA_BLUE
BLUE_BRIGHT = "#1F45A8"
PLOT_DARK = DATA_BLUE       # matched / explained peaks
PLOT_LIGHT = "#B9C3D6"      # unmatched / unexplained peaks
PLOT_FONT = '"IBM Plex Sans",-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif'
PLOT_MONO = '"IBM Plex Mono",ui-monospace,"SF Mono",Menlo,Consolas,monospace'


def _rgb(hex_color, mix=0.0):
    """'#RRGGBB' -> RDKit (r, g, b) in 0..1, mixed `mix` of the way towards white (a soft highlight tint)."""
    h = hex_color.lstrip("#")
    return tuple(round(int(h[i:i + 2], 16) / 255 * (1 - mix) + mix, 3) for i in (0, 2, 4))


# ---- RDKit atom highlights: the data colours as soft tints, so the structure stays readable on top of them
HL_ACID = _rgb(DATA_CORAL, .55)
HL_BASE = _rgb(DATA_TEAL, .55)
HL_BOTH = _rgb(DATA_VIOLET, .5)
HL_STRONG = _rgb(DATA_BLUE, .5)
HL_WEAK = _rgb(DATA_AMBER, .5)
HL_REACTIVE = _rgb(DATA_TEAL, .5)
