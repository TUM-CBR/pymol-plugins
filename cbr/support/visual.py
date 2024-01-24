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
        self.__fasta_entries : fasta.FastaSequences = fasta.parse_fast_meta(self.__fasta_text_edit.toPlainText())
        valid_entries = self.__fasta_entries.sequences
        self.__fasta_combo.clear()
        self.__fasta_combo.insertItems(0, (entry.id for entry in valid_entries))
        self.sequences_changed.emit()

    @property
    def selected_entry(self):
        return self.__fasta_combo.currentText()

    def selected_sequence(self) -> Tuple[str, str]:
        return next((name, sequence) for (name, sequence) in self.get_items() if name == self.selected_entry)

    def get_items(self) -> Iterable[Tuple[str, str]]:

        entries = self.__fasta_entries
        for exn in entries.exceptions:
            raise exn

        for item in entries.sequences:
            yield item.as_tuple()

    def get_items_meta(self):
        return self.__fasta_entries

    @pyqtSlot()
    def __on_text_changed(self):
        self.__throttle.trigger()

def as_fasta_selector(fastaTextEdit : QPlainTextEdit, fastaCombo : QComboBox) -> FastaSelector:
    return FastaSelector(fastaTextEdit, fastaCombo)