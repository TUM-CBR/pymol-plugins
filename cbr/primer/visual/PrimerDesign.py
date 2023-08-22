from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QDoubleValidator
from PyQt5.QtWidgets import QWidget
import re
from typing import Optional

from ...core.Qt.QtWidgets import show_error
from ..Operations import Operations
from .Ui_PrimerDesign import Ui_PrimerDesign

dna_re = re.compile(r"(?P<left>(A|C|T|G)+)\[(?P<design>(A|C|T|G)+)\](?P<right>(A|C|T|G)+)", re.IGNORECASE)

class PrimerDesign(QWidget):

    def __init__(self, operations : Optional[Operations] = None):
        super(PrimerDesign, self).__init__()
        self.__ui = Ui_PrimerDesign()
        self.__ui.setupUi(self)
        self.__operations = operations or Operations()

        self.__ui.designPrimersButton.clicked.connect(self.__on_design_primers)
        self.__ui.desiredTm.setValidator(QDoubleValidator(0, 100, 2))

    @pyqtSlot
    def __on_design_primers(self):

        sequence = self.__ui.sequenceInputText.toPlainText()
        sections = dna_re.match(sequence)

        if sections is None:
            show_error(self, "Invalid sequence", description=invalid_sequence_error)
            return

        try:
            target_tm = float(self.__ui.desiredTm.text())
        except ValueError:
            show_error(self, "Enter a valid Tm")
            return
        


        

invalid_sequence_error =\
"""The sequence must only consists of C,T,G,A and must contain a region surrounded
by brackets ([...CTGA...]) which will be the region for which primers will be
designed.
"""