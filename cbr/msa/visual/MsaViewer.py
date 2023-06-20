from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QFileDialog, QPushButton, QWidget
from typing import Dict

from ...core.Context import Context
from ...core import visual
from ...clustal import msa
from .Ui_MsaViewer import Ui_MsaViewer

class MsaViewer(QWidget):

    def __init__(self, context : Context, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.__ui = Ui_MsaViewer()
        self.__ui.setupUi(self)
        visual.as_structure_selector(self.__ui.structuresCombo, self.__ui.colorConservedButton)
        self.__set_sequences({})

    def __set_sequences(self, sequences : Dict[str, str]):
        self.__sequences = dict(sequences)
        self.__ui.sequenceCombo.clear()
        self.__ui.sequenceCombo.addItems(self.__sequences.keys())

    @pyqtSlot()
    def on_selectMsaButton_clicked(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly  # Optionally set file dialog options

        # Show the file dialog and get the selected file
        file_path, _ = QFileDialog.getOpenFileName(self, "Open File", "", "FASTA alignment Files (*.*)", options=options)

        if not file_path:
            return

        print("Selected file:", file_path)
        self.__set_sequences(msa.parse_alignments(file_path))