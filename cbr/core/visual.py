from PyQt5.QtWidgets import QComboBox, QPushButton
from pymol import cmd as pymol

def as_structure_selector(
    structures_combo : QComboBox,
    refresh_button : QPushButton
    ) -> None:

    def on_refresh_clicked():
        structures_combo.clear()

        for name in pymol.get_names():
            for chain in pymol.get_chains(name):
                structures_combo.addItem("%s/%s" % (name, chain), (name, chain))

    refresh_button.clicked.connect(on_refresh_clicked)