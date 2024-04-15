from Bio.Data import IUPACData
from colorsys import hls_to_rgb, rgb_to_hsv, hsv_to_rgb
import math
from PyQt5.QtGui import QColor
from typing import cast, List, NamedTuple, Tuple, Union

RgbColor = Tuple[int, int, int]

def hls_to_rgb_int(h : float, l : float, s : float) -> RgbColor:
    r, g, b = hls_to_rgb(h, l, s)
    return (int(r * 255), int(g * 255), int(b * 255))

def generate_colors(n : int) -> List[RgbColor]:
    colors = []
    for i in range(n):
        hue = i / float(n)
        lightness = 0.5
        saturation = 1
        colors.append(hls_to_rgb_int(hue, lightness, saturation))
    return colors

distinct_colors_15 : List[RgbColor] = generate_colors(15)

distinct_colors_5 : List[RgbColor] = generate_colors(5)

residue_letter_colors = dict(
    zip(
        IUPACData.protein_letters,
        generate_colors(len(IUPACData.protein_letters))
    )
)

def to_hex(i : int, desired_length = 2) -> str:
    value = hex(i).replace("0x", "")
    padding = desired_length - len(value)
    return "0"*padding + value

def to_pymol_color(color : Union[QColor, RgbColor]) -> str:

    if isinstance(color, QColor):
        rgb = cast(RgbColor, color.toRgb().getRgb()[0:3])
    else:
        rgb = color        

    return "0x" + "".join(map(to_hex, rgb))

def pymol_color_tuple_to_pymol_color(colors : Tuple[float, float, float]) -> str:
    r,g,b = colors
    return to_pymol_color(
        (
            int(255*r),
            int(255*g),
            int(255*b)
        )
    )

colors_15 = [
    to_pymol_color(rgb)
    for rgb in distinct_colors_15
]

def get_qt_color(item : int, options = distinct_colors_15) -> List[int]:
    return list(options[item % len(options)]) + [100]

def get_color(item : int, options = colors_15) -> str:
    return options[item % len(options)]

def get_color_rgb(item : int, options = distinct_colors_15) -> RgbColor:
    return options[item % len(options)]

def change_intensity(rgb : RgbColor, factor: float) -> RgbColor:

    if factor > 1:
        factor = 1

    if factor < 0:
        factor = 0

    (r,g,b) = rgb
    return (int(r*factor), int(g*factor), int(b*factor))

def int_to_rgb(value : int) -> RgbColor:
    blue = value & 0xFF
    green = (value >> 8) & 0xFF
    red = (value >> 16) & 0xFF
    return (red, green, blue)

def int_to_pymol_color(value : int) -> str:
    return to_pymol_color(int_to_rgb(value))

def norm(color : RgbColor) -> float:
    return math.sqrt(sum(x * x for x in color))

def delta(c1: RgbColor, c2: RgbColor) -> RgbColor:
    return cast(RgbColor, tuple((ca - cb) for (ca, cb) in zip(c2, c1)))

def scale(factor: float, c: RgbColor) -> RgbColor:
    return cast(RgbColor, tuple(factor * i for i in c))

def add(c1 : RgbColor, c2 : RgbColor) -> RgbColor:
    return cast(RgbColor, tuple(ca + cb for (ca, cb) in zip(c1, c2)))

class ColorRange(NamedTuple):

    base_color : RgbColor

    def get_color(self, scale_factor: float) -> RgbColor:

        if scale_factor < 0 or scale_factor > 1:
            raise ValueError("The 'distance' argument must be between 0 and 1.")

        rgb = self.base_color
        r, g, b = [x/255.0 for x in rgb] # Normalize to [0, 1]
        h, s, v = rgb_to_hsv(r, g, b)
        s *= scale_factor # Reduce the value (brightness)
        r, g, b = hsv_to_rgb(h, s, v)
        return (
            int(r * 255),
            int(g * 255),
            int(b * 255)
        )

def color_range_scale(color: RgbColor) -> ColorRange:
    return ColorRange(
        base_color=color
    )

class ColorSpread(NamedTuple):
    low_color_hsv : Tuple[float, float, float]
    high_color_hsv : Tuple[float, float, float]

    def get_color(self, scale_factor: float) -> RgbColor:
        if scale_factor < 0 or scale_factor > 1:
            raise ValueError("The 'distance' argument must be between 0 and 1.")

        h,s,v = [
            start + ((end - start) * scale_factor)
            for end,start in zip(self.high_color_hsv, self.low_color_hsv)
        ]

        r, g, b = hsv_to_rgb(h, s, v)
        return (
            int(r * 255),
            int(g * 255),
            int(b * 255)
        )

def color_spread(rgb1 : RgbColor, rgb2 : RgbColor):
    return ColorSpread(
        low_color_hsv=rgb_to_hsv(*[c / 255 for c in rgb1]),
        high_color_hsv=rgb_to_hsv(*[c / 255 for c in rgb2])
    )