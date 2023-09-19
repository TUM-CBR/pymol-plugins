from colorsys import hls_to_rgb
from typing import List, Tuple

RgbColor = Tuple[int, int, int]

def hls_to_rgb_int(h : float, l : float, s : float) -> RgbColor:
    r, g, b = hls_to_rgb(h, l, s)
    return (int(r * 255), int(g * 255), int(b * 255))

def generate_colors(n : int) -> List[RgbColor]:
    colors = []
    for i in range(n):
        hue = i / float(n)
        lightness = 0.5
        saturation = 0.8
        colors.append(hls_to_rgb_int(hue, lightness, saturation))
    return colors

distinct_colors_15 : List[RgbColor] = generate_colors(15)

def to_hex(i : int, desired_length = 2) -> str:
    value = hex(i).replace("0x", "")
    padding = desired_length - len(value)
    return "0"*padding + value

def to_pymol_color(rgb : RgbColor) -> str:
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