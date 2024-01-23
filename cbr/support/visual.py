from typing import Any, Iterable, Tuple
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot
from PyQt5.QtWidgets import QComboBox, QPlainTextEdit

from ..core.Qt.QtCore import Throttle
from ..clustal import fasta

class FastaSelector(QObject):

    sequences_changed = pyqtSignal()

    def __init__(
        self,
        fastaTextEdit : QPlainTextEdit,
        fastaCombo : QComboBox,
        *args: Any,
        **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)

        self.__throttle = Throttle(500, self.__on_fasta_text_changed)
        self.__fasta_text_edit = fastaTextEdit
        self.__fasta_text_edit.textChanged.connect(self.__on_text_changed)
        self.__fasta_combo = fastaCombo

    def __on_fasta_text_changed(self):
        self.__fasta_entries = fasta.parse_fasta(self.__fasta_text_edit.toPlainText())
        valid_entries = (item for item in self.__fasta_entries if isinstance(item, tuple))
        self.__fasta_combo.clear()
        self.__fasta_combo.insertItems(0, (key for (key, _) in valid_entries))
        self.sequences_changed.emit()

    @property
    def selected_entry(self):
        return self.__fasta_combo.currentText()

    def selected_sequence(self) -> Tuple[str, str]:
        return next((name, sequence) for (name, sequence) in self.get_items() if name == self.selected_entry)

    def get_items(self) -> Iterable[Tuple[str, str]]:
        for item in self.__fasta_entries:
            if(isinstance(item, Exception)):
                raise item
            yield item

    @pyqtSlot()
    def __on_text_changed(self):
        self.__throttle.trigger()

def as_fasta_selector(fastaTextEdit : QPlainTextEdit, fastaCombo : QComboBox) -> FastaSelector:
    return FastaSelector(fastaTextEdit, fastaCombo)