from io import StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QWidget
from typing import Any, Iterable, List, Optional

from ...display.sequence import RESIDUE_COLORS
from .Ui_FastaViewer import Ui_FastaViewer

class FastaViewerModel(QAbstractTableModel):

    def __init__(
        self,
        sequences: Optional[Iterable[SeqRecord]] = None
    ):
        super().__init__()
        self.set_sequences(sequences, notify=False)

    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.__sequences)
    
    def columnCount(self, parent: Optional[QModelIndex] = None) -> int:
        return self.__length
    
    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        
        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Vertical:
            seq = self.__sequences[section]
            return f"{seq.id} {seq.name}"
        return super().headerData(section, orientation, role)
    
    def __get_sequence(self, index: QModelIndex) -> SeqRecord:
        return self.__sequences[index.row()]
    
    def __get_residue(self, index: QModelIndex) -> Optional[str]:
        seq = self.__get_sequence(index)
        col = index.column()

        if len(seq) > col:
            return seq[col].upper()
        else:
            return None
        
    def __get_data_role(self, index: QModelIndex) -> Optional[str]:
        return self.__get_residue(index)
    
    def set_sequences(
        self,
        sequences: Optional[Iterable[SeqRecord]],
        notify: bool = True
    ):
        if sequences is None:
            self.__sequences = []
            self.__length = 0
        else:
            self.__sequences: List[SeqRecord] = list(sequences)
            self.__length = max(
                len(seq)
                for seq in self.__sequences
            )

        if notify:
            self.modelReset.emit()

    def __get_background(self, index: QModelIndex) -> Optional[QColor]:
        residue = self.__get_residue(index)
        
        if residue is None:
            return None
        
        return RESIDUE_COLORS.get(residue)

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if not index.isValid():
            return None
        
        if role == Qt.ItemDataRole.DisplayRole:
            return self.__get_data_role(index)
        if role == Qt.ItemDataRole.BackgroundColorRole:
            return self.__get_background(index)

        return None

class FastaViewer(QWidget):

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        fasta_file: Optional[str] = None
    ):
        super().__init__(parent)
        self.__ui = Ui_FastaViewer()
        self.__ui.setupUi(self)

        self.__model = FastaViewerModel(None)
        self.__ui.sequencesTable.setModel(self.__model)

        self.set_fasta(fasta_file)

    def set_fasta(self, fasta_file: Optional[str]):

        if fasta_file is None:
            self.__model.set_sequences(fasta_file)
            self.__ui.sequencesTextEdit.setText("")
            return

        sequences: Any = SeqIO.parse(fasta_file, format="fasta")
        self.__model.set_sequences(sequences)

        for col in range(self.__ui.sequencesTable.colorCount()):
            self.__ui.sequencesTable.setColumnWidth(col, 1)

        #with StringIO() as seq_buffer:
        #    SeqIO.write(sequences, seq_buffer, format="fasta")
        #    seq_buffer.seek(0)
        #    self.__ui.sequencesTextEdit.setText(seq_buffer.read())

