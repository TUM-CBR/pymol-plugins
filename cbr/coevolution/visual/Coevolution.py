from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import numpy as np
from os import path
from typing import Dict, List, NamedTuple, Optional, Sequence

from PyQt5.QtCore import QModelIndex, QObject, pyqtSlot
from PyQt5.QtWidgets import QAbstractItemView, QDialog, QWidget
from PyQt5.QtGui import QColor

from ...core.Context import Context
from ...core.Qt.QtWidgets import show_error, show_exception
from ...core.pymol.structure import StructureSelection
from ...clustal.Clustal import Clustal, get_clustal_from_context
from ...extra.CbrExtraInteractive import CbrExtraInteractive, CbrProcessExit, MessageQueingPolicy, run_interactive
from ...support.display.sequence import RESIDUE_COLORS
from ...support.msa.MsaSelector import MsaSelector
from ...support.structure.AbstractCompositeTableModel import DEFAULT_PYMOL_ATTRIBUTES, AbstractCompositeTableModel, AbstractRecordView, PymolRecordAttributes, ViewHeaderSpec, ViewRecords, ViewRecordAttributes
from ...support.structure.MultiStructureSelector import MultiStructureSelector
from ...support.structure.StructuresAlignmentMapper import StructuresAlignmentMapper
from ...support.structure.StructurePositionView  import PositionEntry, PositionsWithAttributes, SetStructureArg, StructurePositionView
from ..data import *
from .Ui_coevolution import Ui_Coevolution

class ReferenceStructureEntry(NamedTuple):
    structure: StructureSelection
    reference_sequence_id: str
    msa_to_resv: Dict[int, int]

    @classmethod
    def new(
        cls,
        structure: StructureSelection,
        reference_sequence_id: str,
        msa_mapper: StructuresAlignmentMapper
    ) -> 'ReferenceStructureEntry':
        mapping = msa_mapper.map_structure(structure, reference_sequence_id)
        msa_to_resv = {
            pos: resv
            for resv, pos in mapping.resv_to_msa().items()
            if pos is not None
        }

        return ReferenceStructureEntry(
            structure,
            reference_sequence_id,
            msa_to_resv
        )

SCORE_COLOR_LOW = np.array([128,128,0], dtype=np.uint8)
SCORE_COLOR_HIGH = np.array([0, 255, 0], dtype=np.uint8)
SCORE_COLOR_SPREAD = SCORE_COLOR_HIGH - SCORE_COLOR_LOW

class CoevolutionResultEntry(NamedTuple):
    position: int
    coevolution_result: CoevolutionEntry

    @classmethod
    def color(cls, score: float):
        result = (SCORE_COLOR_LOW + (SCORE_COLOR_SPREAD * score)).astype(np.uint8)
        return QColor(
            result[0],
            result[1],
            result[2]
        )

    def to_position_entry(self) -> PositionEntry:

        return PositionEntry(
            positions_with_attributes = [
                PositionsWithAttributes(
                    positions=[self.position],
                    structure_color=self.color(self.coevolution_result.score)
                )
            ]
        )

class CoevolutionScoreView(AbstractRecordView[CoevolutionResultEntry]):

    K_RESIDUE_1 = "local residue"
    K_RESIDUE_2 = "other residue"
    K_SCORE = "score"
    K_OCCURRENCE = "occurrence"
    K_EXCLUSIVITY = "exclusivity score"
    K_SYMMETRY = "symmetry score"
    K_CONFIDENCE = "confidence score"

    HEADERS = [
        K_RESIDUE_1,
        K_RESIDUE_2,
        K_SCORE,
        K_OCCURRENCE,
        K_EXCLUSIVITY,
        K_SYMMETRY,
        K_CONFIDENCE
    ]

    SCORE_LOW = np.array([255, 0, 0, 128])
    SCORE_HIGH = np.array([0, 255, 0, 128])
    SCORE_SPREAD = SCORE_HIGH - SCORE_LOW

    PYMOL_ATTRIBUTES = [DEFAULT_PYMOL_ATTRIBUTES] * len(HEADERS)

    def __get_score_color(self, score: float) -> QColor:
        score = min(1, max(0, score))
        return QColor(*(self.SCORE_LOW + self.SCORE_SPREAD*score))

    def __to_qt_attributes(self, record: CoevolutionResultEntry) -> List[ViewRecordAttributes]:
        entry = record.coevolution_result
        return [
            ViewRecordAttributes(entry.residue_1, RESIDUE_COLORS[entry.residue_1.upper()]),
            ViewRecordAttributes(entry.residue_2, RESIDUE_COLORS[entry.residue_2.upper()]),
            ViewRecordAttributes(entry.score, self.__get_score_color(entry.score)),
            ViewRecordAttributes(entry.score_occurence, self.__get_score_color(entry.score_occurence)),
            ViewRecordAttributes(entry.score_exclusivity, self.__get_score_color(entry.score_exclusivity)),
            ViewRecordAttributes(entry.score_symmetry, self.__get_score_color(entry.score_symmetry)),
            ViewRecordAttributes(entry.score_confidence, self.__get_score_color(entry.score_confidence))
        ]
    
    def __to_pymol_attributes(self, _record: CoevolutionResultEntry) -> List[PymolRecordAttributes]:
        return self.PYMOL_ATTRIBUTES

    def attributes(self, records: Sequence[CoevolutionResultEntry]) -> ViewRecords:
        return ViewRecords(
            headers=[ViewHeaderSpec(h) for h in self.HEADERS],
            qt_attributes=[self.__to_qt_attributes(r) for r in records],
            pymol_attributes=[self.__to_pymol_attributes(r) for r in records]
        )


class CoevolutionResultTableModel(AbstractCompositeTableModel[CoevolutionResultEntry]):

    def __init__(
        self,
        clustal: Clustal,
        msa: Optional[MultipleSeqAlignment] = None,
        parent: Optional[QObject] = None
    ) -> None:
        super().__init__(parent)

        self.__records: Sequence[CoevolutionResultEntry] = []

        self.__structure_view: StructurePositionView[CoevolutionResultEntry] = StructurePositionView(
            clustal,
            lambda _msa, model: model.to_position_entry(),
            msa=msa,
            parent=self
        )

        self.__score_view = CoevolutionScoreView(self)

    def set_structures(self, structures: Sequence[SetStructureArg]):
        self.__structure_view.set_structures(structures)

    def __views__(self) -> Sequence[AbstractRecordView[CoevolutionResultEntry]]:
        return [self.__structure_view, self.__score_view]
    
    def set_alignment(self, msa: MultipleSeqAlignment):
        self.__structure_view.set_alignment(msa)

    def set_results(self, results: CoevolutionPosition):

        positions = list(results.by_position.keys())
        positions.sort(
            key=lambda pos: results.by_position[pos].score,
            reverse=True
        )
        self.__records = [
            CoevolutionResultEntry(
                position=position,
                coevolution_result=results.by_position[position]
            )
            for position in positions
        ]

        self.records_reset.emit()

    def __records__(self) -> Sequence[CoevolutionResultEntry]:
        return self.__records

class CoevolutionOverviewEntry(NamedTuple):
    position: int

    def to_position_entry(self) -> PositionEntry:

        return PositionEntry(
            positions_with_attributes = [
                PositionsWithAttributes(
                    positions=[self.position]
                )
            ]
        )

class CoevolutionOverviewModel(AbstractCompositeTableModel[CoevolutionOverviewEntry]):

    def __init__(
        self,
        clustal: Clustal,
        msa: Optional[MultipleSeqAlignment] = None,
        parent: Optional[QObject] = None
    ) -> None:
        super().__init__(parent)
        self.__records: Sequence[CoevolutionOverviewEntry] = []
        self.__msa_view: StructurePositionView[CoevolutionOverviewEntry] = StructurePositionView(
            clustal,
            lambda _msa, model: model.to_position_entry(),
            msa=msa
        )

    def set_alignment(self, msa: MultipleSeqAlignment):

        self.__records = [
            CoevolutionOverviewEntry(pos)
            for pos in range(0, msa.get_alignment_length())
        ]

        self.__msa_view.set_alignment(msa)

    def set_structures(self, structures: Sequence[SetStructureArg]):
        self.__msa_view.set_structures(structures)

    def __records__(self) -> Sequence[CoevolutionOverviewEntry]:
        return self.__records

    def __views__(self) -> Sequence[AbstractRecordView[CoevolutionOverviewEntry]]:
        return [self.__msa_view]

CoevolultionHandler = CbrExtraInteractive[InteractiveResponse, InteractiveRequest]

class Coevolution(QWidget):

    MSA_FILENAME = "tmp.aln"

    def __init__(
        self,
        context: Context
    ) -> None:
        super().__init__()
        self.__ui = Ui_Coevolution()
        self.__ui.setupUi(self)

        clustal = get_clustal_from_context(context)

        self.__msa_seletor = MsaSelector(
            parent=self,
            select_file_button=self.__ui.selectAlignmentButton,
            selected_file_label=self.__ui.selectedAlignmentLabel
        )
        self.__msa_seletor.msa_file_selected.connect(self.__on_msa_selected)
        self.__msa: Optional[MultipleSeqAlignment] = None

        self.__alignment_model = CoevolutionOverviewModel(clustal)
        self.__ui.alignmentTable.setModel(self.__alignment_model)
        self.__ui.alignmentTable.doubleClicked.connect(self.__on_position_selected)
        self.__ui.alignmentTable.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)

        self.__results_model = CoevolutionResultTableModel(clustal)
        self.__ui.detailsTable.setModel(self.__results_model)

        self.__working_directory = context.create_temporary_directory()
        self.__coevolution_manager: Optional[CoevolultionHandler] = None

        self.__results_by_position: Dict[int, CoevolutionPosition] = {}
        self.__structure_selector: MultiStructureSelector = MultiStructureSelector(1, None)
        self.__ui.selectStructureButton.clicked.connect(self.__on_slelect_structure)

        self.__set_busy(False)

    @pyqtSlot()
    def __on_slelect_structure(self):

        msa = self.__msa

        if msa is None:
            show_error(self, "Missing Alignment", "You must first select an alignment")
        elif self.__structure_selector.exec_with_sequences([seq.id for seq in msa]) == QDialog.DialogCode.Accepted:

            selected = [
                SetStructureArg(seq.structure, seq.sequence)
                for seq in self.__structure_selector.selected_sequences()
            ]
            self.__alignment_model.set_structures(selected)
            self.__results_model.set_structures(selected)

    def __set_busy(self, busy: bool):

        self.__ui.busyProgress.setVisible(busy)
        self.__ui.alignmentTable.setEnabled(not busy)
        self.__ui.detailsTable.setEnabled(not busy)

    @pyqtSlot(QModelIndex)
    def __on_position_selected(self, index: QModelIndex):

        manager = self.__coevolution_manager
        model = self.__alignment_model
        record = model.get_record(index)

        if manager is None or record is None:
            return
        
        position = record.position

        if position in self.__results_by_position:
            self.__results_model.set_results(self.__results_by_position[position])
        else:
            manager.send_message(
                InteractiveRequest(
                    query=Query(
                        positions=[record.position],
                        max_results=30
                    )
                )
            )

            self.__set_busy(True)

    def __run_interactive_process(self, msa: MultipleSeqAlignment):

        previous = self.__coevolution_manager

        if previous is not None:
            previous.stop()
            previous.dispose_subscriptions()
            self.__set_busy(False)

        self.__alignment_model.set_alignment(msa)
        self.__results_model.set_alignment(msa)

        msa_tmp_file = path.join(self.__working_directory, self.MSA_FILENAME)

        AlignIO.write(msa, msa_tmp_file, 'fasta')
        self.__coevolution_manager = manager = run_interactive(
            [
                "coevolution",
                "interactive",
                "--input-msa", msa_tmp_file
            ],
            interactive_response_parser,
            interactive_request_serializer,
            MessageQueingPolicy.HANDLE_ALL
        )

        manager.observe_status() \
            .for_each(
                action=self.__on_status,
                on_error=self.__on_process_error
            )

        manager.observe_values() \
            .for_each(
                action=self.__on_result,
                on_error=self.__on_error
            )

    def __on_status(self, value: CbrProcessExit):
        self.__set_busy(False)

    def __set_coevolution_results(self, results: CoevolutionResults):

        for position, value in results.positions.items():
            self.__results_by_position[position] = value

        index = self.__ui.alignmentTable.selectedIndexes()

        self.__set_busy(False)

        if len(index) > 0:
            self.__on_position_selected(index[0])

    def __on_result(self, value: InteractiveResponse):

        if value.coevolution is not None:
            self.__set_coevolution_results(value.coevolution)
        else:
            show_error(
                self,
                "Unknown Response",
                "Received an unknown result when running coevolution."
            )

    def __on_process_error(self, error: Exception):
        show_exception(self, error)

    def __on_error(self, error: Exception):
        show_exception(self, error)
        self.__set_busy(False)

    @pyqtSlot(object)
    def __on_msa_selected(self, msa: MultipleSeqAlignment):

        self.__msa = msa
        self.__alignment_model.set_alignment(msa)
        self.__results_model.set_alignment(msa)
        self.__run_interactive_process(msa)
