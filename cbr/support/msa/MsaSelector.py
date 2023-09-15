from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot
from PyQt5.QtWidgets import QFileDialog, QLabel, QPushButton, QWidget
from os import path
from typing import Callable, Dict

from ...core.stdlib import get_or_raise

from ...clustal import fasta
from ...clustal import msa

Msa = Dict[str, str]

def fasta_parser(in_file : str) -> Msa:

    return dict(
        (k,v)
        for item in fasta.parse_fasta(in_file)
        for k,v in [get_or_raise(item)]
    )

def clustal_parser(in_file : str) -> Msa:
    return msa.parse_alignments(in_file)

MSA_PARSERS : Dict[str, Callable[[str], Msa]] = {
    "fasta" : fasta_parser,
    "fa" : fasta_parser,
    "clustal" : clustal_parser
}

MSA_EXTENSIONS = " ".join(f"*.{extension}" for extension in MSA_PARSERS.keys())

class MsaSelector(QObject):

    msa_file_selected = pyqtSignal(object)

    def __init__(
        self,
        parent : QWidget,
        select_file_button : QPushButton,
        selected_file_label : QLabel
    ):

        super(MsaSelector, self).__init__()
        self.__select_file_button = select_file_button
        self.__selected_file_label = selected_file_label
        self.__parent = parent

        self.__select_file_button.clicked.connect(
            self.__on_select_file
        )

        self.__msa = None

    @property
    def msa(self):
        return self.__msa

    @pyqtSlot()
    def __on_select_file(self):

        msa_file,_ = QFileDialog.getOpenFileName(
            self.__parent,
            "Select MSA File",
            "",
            f"MSA files ({MSA_EXTENSIONS})"
        )

        extension = path.splitext(msa_file)[1][1:]
        self.__selected_file_label.setText(path.basename(msa_file))
        self.__msa = MSA_PARSERS[extension](msa_file)
        self.msa_file_selected.emit(self.__msa)
