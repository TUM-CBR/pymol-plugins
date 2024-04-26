from PyQt5.QtGui import QColor

from ...core import color

RESIDUE_COLORS = dict(
    [(resi, QColor(*rgb, 128)) for resi, rgb in color.residue_letter_colors.items()] \
        + [("X", QColor(200, 200, 200, 128))]
    )
