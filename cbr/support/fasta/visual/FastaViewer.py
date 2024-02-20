from io import StringIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pymol import cmd
from PyQt5.QtCore import pyqtSlot, QAbstractTableModel, QItemSelection, QModelIndex, Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QDialog, QWidget
from typing import Any, Callable, cast, Dict, Iterable, List, NamedTuple, Optional, Tuple, Union

from ....clustal.Clustal import Clustal
from ....core.pymol.structure import get_pdb_dominant_color_index, get_selection_sequence_index, StructureSelection
from ...msa.visual.MsaStructureSelector import MsaStructureSelector
from ...display.sequence import RESIDUE_COLORS
from .Ui_FastaViewer import Ui_FastaViewer

SequenceStructures = Union[List[Optional[StructureSelection]], Optional[StructureSelection]]

class ResidueEntry(NamedTuple):
    seq_ix: Optional[int]
    seq_residue: Optional[str]
    pdb_resv: Optional[int]
    pdb_residue: Optional[List[str]]

class SequenceData(NamedTuple):
    sequence: SeqRecord
    sequence_structure: Optional[StructureSelection] = None
    structure_alignment: Optional[MultipleSeqAlignment] = None
    pdb_seq: Optional[str] = None
    pdb_ix_to_resv: Optional[List[int]] = None
    pdb_color_ix: Optional[int] = None
    group_header: bool = False

    SEQ_IX : int = 0
    PDB_IX :int = 1

    def set_structure(
        self,
        clustal: Clustal,
        structure: Optional[StructureSelection]
    ) -> 'SequenceData':

        if structure is None:
            alignment = None
            pdb_seq = None
            pdb_ix_to_resv = None
            pdb_color_ix = None
        else:
            seq_index = get_selection_sequence_index(structure)
            alignment = clustal.run_msa_seqs([
                self.sequence,
                SeqRecord(
                    seq = Seq("".join(seq_index.values())),
                    id = structure.show()
                )
            ])

            pdb_ix_to_resv = list(seq_index.keys())
            pdb_ix_to_resv.sort()
            pdb_seq = [seq_index[i] for i in pdb_ix_to_resv]
            pdb_color_ix = get_pdb_dominant_color_index(structure)

        return self._replace(
            sequence_structure = structure,
            structure_alignment = alignment,
            pdb_seq = pdb_seq,
            pdb_ix_to_resv = pdb_ix_to_resv,
            pdb_color_ix = pdb_color_ix
        )
    
    def to_pdb_resv(self, ix: int) -> Optional[int]:
        msa = self.structure_alignment
        pdb_ix_to_resv = self.pdb_ix_to_resv

        if msa is None or pdb_ix_to_resv is None:
            return None

        aligment = msa.alignment
        indexes = aligment.indices
        r_indexes = aligment.inverse_indices

        seq_to_msa = r_indexes[self.SEQ_IX][ix]

        if seq_to_msa < 0:
            return None
        
        msa_to_pdb: int = indexes[self.PDB_IX][seq_to_msa]

        if msa_to_pdb < 0:
            return None

        return pdb_ix_to_resv[msa_to_pdb]

    def iterate_residues(self) -> Iterable[ResidueEntry]:
        msa = self.structure_alignment
        mapping = self.pdb_ix_to_resv
        pdb_seq = self.pdb_seq
        assert msa is not None and mapping is not None and pdb_seq is not None, \
            "This function can only be called if there is a structure mapping"

        aligment = msa.alignment
        indexes = aligment.indices
        seq_seq = self.sequence

        for (seq_ix, pdb_ix) in zip(indexes[self.SEQ_IX], indexes[self.PDB_IX]):

            if seq_ix >= 0:
                seq_residue: Optional[str] = seq_seq[seq_ix]
            else:
                seq_residue = None
                seq_ix = None

            if pdb_ix >= 0:
                pdb_resv: Optional[int] = mapping[pdb_ix]
                pdb_residue : Optional[str] = pdb_seq[pdb_ix]
            else:
                pdb_resv = None
                pdb_residue = None

            yield ResidueEntry(
                seq_ix = seq_ix,
                seq_residue = cast(Any, seq_residue),
                pdb_resv = cast(Any, pdb_resv),
                pdb_residue = cast(Any, pdb_residue)
            )

def get_groups(seq: SeqRecord) -> List[SequenceData]:
    seqs = seq.seq.rsplit('/')

    if len(seqs) == 1:
        return [SequenceData(seq)]

    def make_seq_record(i: int, segment: Seq):
        record = SeqRecord(
            seq=segment,
            id=seq.id,
            description=seq.description
        )

        return SequenceData(
            record,
            # First sequence of the group is the header
            group_header=i == 0
        )
        
    return [make_seq_record(i, new_seq) for i,new_seq in enumerate(seqs)]

def default_render_header(seq: SeqRecord, is_group_header: bool) -> str:
    if is_group_header:
        return f"{seq.id} {seq.name}"
    else:
        return "/"
    
RenderHeaderFunction = Callable[[SeqRecord, bool], str]

class FastaViewerModel(QAbstractTableModel):

    SEQ_STRUCTURE_COL = "Reference Structure"
    MEATA_COLS = [SEQ_STRUCTURE_COL]
    COLOR_PREFIX = "cbr_fasta_"

    def __init__(
        self,
        sequences: Optional[Iterable[SeqRecord]] = None,
        clustal: Optional[Clustal] = None
    ):
        super().__init__()
        self.set_sequences(sequences, notify=False)
        self.__clustal = clustal is not None and clustal or Clustal()
        self.__render_header: RenderHeaderFunction = default_render_header

    def set_render_header(self, render_header: RenderHeaderFunction):
        self.__render_header = render_header
        self.modelReset.emit()

    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.__sequences)
    
    def columnCount(self, parent: Optional[QModelIndex] = None) -> int:
        return self.__length + len(self.MEATA_COLS)
    
    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        
        if role != Qt.ItemDataRole.DisplayRole:
            return super().headerData(section, orientation, role)
        
        if orientation == Qt.Orientation.Vertical:
            entry = self.__sequences[section]
            return self.__render_header(
                entry.sequence,
                entry.group_header
            )
        elif section < len(self.MEATA_COLS):
            return self.MEATA_COLS[section]
        else:
            return str(section - len(self.MEATA_COLS) + 1)
    
    def __get_entry_index(self, index: QModelIndex) -> int:
        return index.row()
    
    def __get_entry(self, index: QModelIndex) -> SequenceData:
        return self.__sequences[self.__get_entry_index(index)]

    def __get_sequence(self, index: QModelIndex) -> SeqRecord:
        return self.__get_entry(index).sequence
    
    def __get_residue_with_ix(self, index: QModelIndex) -> Optional[Tuple[int, str]]:
        seq = self.__get_sequence(index)
        col = index.column()

        if col < len(self.MEATA_COLS):
            return None

        col -= len(self.MEATA_COLS)
        if len(seq) > col:
            return (col, seq[col].upper())

    def __get_residue(self, index: QModelIndex) -> Optional[str]:

        entry = self.__get_residue_with_ix(index)
        if entry is not None:
            return entry[1]
        else:
            return None
        
    def __get_meta(self, index: QModelIndex) -> Optional[str]:

        col = index.column()
        if col >= len(self.MEATA_COLS):
            return None
        
        col_name = self.MEATA_COLS[col]

        if col_name == self.SEQ_STRUCTURE_COL:
            structure = self.__sequences[index.row()].sequence_structure

            return structure is not None and structure.show() or "<double click to select>"
        
    def __get_data_role(self, index: QModelIndex) -> Optional[str]:
        result: Optional[str] = None

        if (result := self.__get_residue(index)) is not None \
            or (result := self.__get_meta(index)) is not None:
            return result
        
        return None
    
    def get_select_structure_index(self, index: QModelIndex) -> Optional[int]:

        col = index.column()

        if col >= len(self.MEATA_COLS) \
            or self.MEATA_COLS[col] != self.SEQ_STRUCTURE_COL:
            return None
        
        return self.__get_entry_index(index)

    def set_structure(
        self,
        group_ix: int,
        new_structure: SequenceStructures
    ):
        group = self.__sequence_groups[group_ix]

        if isinstance(new_structure, list):
            structures = new_structure
        else:
            structures = [new_structure] * len(group)

        for ix, structure in zip(group, structures):
            self.__sequences[ix] = self.__sequences[ix].set_structure(
                self.__clustal,
                structure
            )

            updated_index = self.index(ix, self.MEATA_COLS.index(self.SEQ_STRUCTURE_COL))
            self.dataChanged.emit(updated_index, updated_index)
    
    def set_sequences(
        self,
        sequences: Optional[Iterable[SeqRecord]],
        notify: bool = True
    ):
        if sequences is None:
            self.__sequences: List[SequenceData] = []
            self.__length = 0
            self.__sequence_groups = []
        else:
            new_sequences: List[SequenceData] = []
            groups: List[List[int]] = []

            items = (
                (group_id, member)
                for group_id, sequence in enumerate(sequences)
                for member in get_groups(sequence)
            )

            for i,(gid, item) in enumerate(items):
                new_sequences.append(item)

                while(len(groups) <= gid):
                    groups.append([])

                groups[gid].append(i)

            self.__sequences: List[SequenceData] = new_sequences
            self.__sequence_groups = groups

            self.__length = max(
                len(seq.sequence)
                for seq in self.__sequences
            )

        if notify:
            self.modelReset.emit()

    def __get_background(self, index: QModelIndex) -> Optional[QColor]:
        residue = self.__get_residue(index)
        
        if residue is None:
            return None
        
        return RESIDUE_COLORS.get(residue)
    
    def __select_model(
        self,
        item: SequenceData
    ) -> bool:
        structure = item.sequence_structure

        if structure is None:
            return False
        
        residues_to_color: Dict[Optional[str], List[int]] = {}
        for entry in item.iterate_residues():
            if entry.pdb_resv is not None \
                and entry.pdb_residue != entry.seq_residue:
                residue = entry.seq_residue

                if residue not in residues_to_color:
                    residues_to_color[residue] = []

                residues_to_color[residue].append(entry.pdb_resv)

        cmd.color(
            "white",
            structure.selection
        )

        for residue, indices in residues_to_color.items():

            if residue is None:
                color = f"{self.COLOR_PREFIX}NONE"
                cmd.set_color(color, (0,0,0))
            else:
                color = f"{self.COLOR_PREFIX}{residue}"
                q_color = RESIDUE_COLORS.get(residue)
                if q_color is None:
                    rgb = (255,255,255)
                else:
                    rgb = (q_color.red(), q_color.green(), q_color.blue())
                cmd.set_color(color, rgb)

            cmd.color(
                color,
                structure.residue_selection(indices)
            )

        return True
    
    def update_colors(
        self,
        selected: QItemSelection,
        deselected: QItemSelection
    ):
        selected_selections = set()
        for index in selected.indexes():

            item = self.__get_entry(index)
            self.__select_model(item)

            if item.sequence_structure is not None:
                selected_selections.add(item.sequence_structure.show())

        for index in deselected.indexes():
            item = self.__get_entry(index)
            structure = item.sequence_structure
            pdb_color = item.pdb_color_ix
            if  structure is None \
                or pdb_color is None \
                or structure.show() in selected_selections:
                continue

            cmd.color(
                pdb_color,
                structure.selection
            )

    def __get_selected_resv(self, item: QModelIndex) -> Optional[Tuple[int, int]]:
        ix = self.__get_entry_index(item)
        entry = self.__sequences[ix]
        seq_ix = None
        resv = None

        if ((seq_ix := self.__get_residue_with_ix(item)) is None \
            or (resv := entry.to_pdb_resv(seq_ix[0])) is None):
            return None
        
        return (ix, resv)

    def update_selected_items(
        self,
        items: Iterable[QModelIndex]
    ):
        # mapping from indexes of the sequences to
        # the resv of each sequence to be selected
        residues_to_select: Dict[int, List[int]] = dict()

        for item in items:
            selection = self.__get_selected_resv(item)

            if selection is None:
                continue

            (ix, resv) = selection
            if ix not in residues_to_select:
                residues_to_select[ix] = []
            residues_to_select[ix].append(resv)

        # we can now easily construct a query with
        # all the selections
        selection = " or ".join(
            structure.residue_selection(resv)
            for (ix, resv) in residues_to_select.items()
            for entry in [self.__sequences[ix]]
            for structure in [entry.sequence_structure] if structure is not None
        )

        cmd.select(
            "sele",
            selection
        )

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
        fasta_file: Optional[str] = None,
        clustal: Optional[Clustal] = None
    ):
        super().__init__(parent)
        self.__ui = Ui_FastaViewer()
        self.__ui.setupUi(self)

        self.__model = FastaViewerModel(None)
        sequences_table = self.__ui.sequencesTable
        sequences_table.setModel(self.__model)
        sequences_table.doubleClicked.connect(self.__on_select_structure)
        selection_model = sequences_table.selectionModel()
        assert selection_model is not None, "Selection model is not eexpected to be None"
        selection_model.selectionChanged.connect(self.__on_selection_changed)

        self.__structure_selection_dialog = MsaStructureSelector(self)
        self.set_fasta(fasta_file)

    def set_render_header(self, render_header: RenderHeaderFunction):
        self.__model.set_render_header(render_header)

    @pyqtSlot(QItemSelection, QItemSelection)
    def __on_selection_changed(self, selected: QItemSelection, deselected: QItemSelection):
        self.__model.update_colors(selected, deselected)
        self.__model.update_selected_items(
            self.__ui.sequencesTable.selectedIndexes()
        )

    @pyqtSlot(QModelIndex)
    def __on_select_structure(self, index: QModelIndex):
        self.__select_structure(index)

    def __select_structure(self, index: QModelIndex):

        seq_index = self.__model.get_select_structure_index(index)
        if seq_index is not None \
            and self.__structure_selection_dialog.exec() == QDialog.DialogCode.Accepted:

            self.__model.set_structure(seq_index, self.__structure_selection_dialog.current_selection())            

    def set_fasta(self, fasta_file: Optional[str]):

        if fasta_file is None:
            self.__model.set_sequences(fasta_file)
            self.__ui.sequencesTextEdit.setText("")
            return

        sequences: List[SeqRecord] = list(
            cast(Any, SeqIO.parse(fasta_file, format="fasta"))
        )

        self.set_sequences(sequences)

    def set_sequences(
        self,
        sequences: List[SeqRecord],
        structure_mappings: Optional[Dict[int, SequenceStructures]] = None
    ):
        self.__model.set_sequences(sequences)
        seqs_table = self.__ui.sequencesTable

        for col in range(self.__model.columnCount()):

            if col < len(self.__model.MEATA_COLS):
                seqs_table.resizeColumnToContents(col)
            else:
                seqs_table.setColumnWidth(col, 1)

        with StringIO() as seq_buffer:
            SeqIO.write(sequences, seq_buffer, format="fasta")
            seq_buffer.seek(0)
            self.__ui.sequencesTextEdit.setText(seq_buffer.read())

        if structure_mappings is None:
            return
        
        for ix, structure in structure_mappings.items():
            self.__model.set_structure(ix, structure)

