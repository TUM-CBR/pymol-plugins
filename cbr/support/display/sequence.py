from PyQt5.QtGui import QColor

from ...core import color

RESIDUE_COLORS = dict(
    [(resi, QColor(*rgb)) for resi, rgb in color.residue_letter_colors.items()] \
        + [("X", QColor(r=200, g=200, b=200))]
    )
