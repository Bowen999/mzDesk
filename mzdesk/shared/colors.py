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
AXIS = "#B9BFCA"            # plot baseline and x ticks (theme.css --line-hover): visible under the sticks, quieter than the data
ACCENT = "#2D5BD0"

# ---- data colours (6, chosen to sit together): series 1 is the accent blue
DATA_BLUE = ACCENT
DATA_AMBER = "#C27C1A"
DATA_TEAL = "#1F8A7E"
DATA_CORAL = "#D0583A"
DATA_VIOLET = "#7457C2"
DATA_GREY = "#A6ADBA"
DATA_SLATE = "#5C6B8A"          # a second line series in the theme family (a dashed line next to a blue one); the warm hues are for chemistry highlights, the alert colour is not a data colour
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


# ---- fragment view (mzdesk/ms2/fragment_view.py): one colour per explained fragment, drawn the same in the structure, the spectrum and the card.
# Equal lightness and chroma in OKLCH (L 0.70, C 0.11), hues spread so neighbours differ: none outweighs another on a white or a dark page.
# FRAGMENT_TEXT is the same hue darker (L 0.52), for the m/z labels on a white page. The whole structure and the unmarked peaks are PURE_INK.
def _oklch_hex(lightness, chroma, hue):
    import math
    a, b = chroma * math.cos(math.radians(hue)), chroma * math.sin(math.radians(hue))
    l, m, s = ((lightness + 0.3963377774 * a + 0.2158037573 * b) ** 3, (lightness - 0.1055613458 * a - 0.0638541728 * b) ** 3,
               (lightness - 0.0894841775 * a - 1.2914855480 * b) ** 3)
    linear = (4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s, -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
              -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s)
    encode = lambda x: 12.92 * x if x <= 0.0031308 else 1.055 * x ** (1 / 2.4) - 0.055
    return "#%02X%02X%02X" % tuple(round(255 * encode(max(0.0, min(1.0, v)))) for v in linear)


FRAGMENT_HUES = (255, 50, 150, 330, 95, 200, 295, 20)          # blue, orange, green, magenta, gold, teal, violet, coral
FRAGMENT_COLORS = tuple(_oklch_hex(0.70, 0.11, h) for h in FRAGMENT_HUES)
FRAGMENT_TEXT = tuple(_oklch_hex(0.52, 0.11, h) for h in FRAGMENT_HUES)
PURE_INK = "#000000"          # the page turns it into pure white in the dark theme (web/src/darkColors.ts)

# ---- RDKit atom highlights: the data colours as soft tints, so the structure stays readable on top of them
HL_ACID = _rgb(DATA_CORAL, .55)
HL_BASE = _rgb(DATA_TEAL, .55)
HL_BOTH = _rgb(DATA_VIOLET, .5)
HL_STRONG = _rgb(DATA_BLUE, .5)
HL_WEAK = _rgb(DATA_AMBER, .5)
HL_REACTIVE = _rgb(DATA_TEAL, .5)
