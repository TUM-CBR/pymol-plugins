from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import numpy as np
from os import path
from typing import Dict, NamedTuple, Optional, Sequence

from PyQt5.QtCore import QModelIndex, QObject, pyqtSlot
from PyQt5.QtWidgets import QWidget
from PyQt5.QtGui import QColor

from ...core.pymol.structure import StructureSelection
from ...core.Context import Context
from ...clustal.Clustal import Clustal, get_clustal_from_context
from ...extra.CbrExtraInteractiveHandler import CbrExtraInteractiveHandler, CbrExtraInteractiveManager, CbrProcessExit, MessageQueingPolicy, run_interactive
from ...support.msa.MsaSelector import MsaSelector
from ...support.structure.AbstractCompositeTableModel import AbstractCompositeTableModel, AbstractRecordView
from ...support.structure.StructuresAlignmentMapper import StructuresAlignmentMapper
from ...support.structure.StructurePositionView  import PositionEntry, PositionsAndColor, StructurePositionView
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

SCORE_COLOR_LOW = np.array([128,128,0], dtype=np.int8)
SCORE_COLOR_HIGH = np.array([0, 255, 0], dtype=np.int8)
SCORE_COLOR_SPREAD = SCORE_COLOR_HIGH - SCORE_COLOR_LOW

class CoevolutionResultEntry(NamedTuple):
    position: int
    coevolution_result: CoevolutionEntry

    @classmethod
    def color(cls, score: float):
        result = (SCORE_COLOR_LOW + (SCORE_COLOR_SPREAD * score)).astype(np.int8)
        return QColor(
            r = result[0],
            g = result[1],
            b = result[2]
        )

    def to_position_entry(self) -> PositionEntry:

        return PositionEntry(
            positions_and_colors = [
                PositionsAndColor(
                    positions=[self.position],
                    structure_color=self.color(self.coevolution_result.score)
                )
            ]
        )

class CoevolutionResultTable(AbstractCompositeTableModel[CoevolutionResultEntry]):

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
            msa=msa
        )

    def __views__(self) -> Sequence[AbstractRecordView[CoevolutionResultEntry]]:
        return [self.__structure_view]
    
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

    def __records__(self) -> Sequence[CoevolutionResultEntry]:
        return self.__records

class CoevolutionOverviewEntry(NamedTuple):
    position: int

    def to_position_entry(self) -> PositionEntry:

        return PositionEntry(
            positions_and_colors = [
                PositionsAndColor(
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

    def __records__(self) -> Sequence[CoevolutionOverviewEntry]:
        return self.__records

    def __views__(self) -> Sequence[AbstractRecordView[CoevolutionOverviewEntry]]:
        return [self.__msa_view]

CoevolultionHandler = CbrExtraInteractiveHandler[InteractiveResponse, InteractiveRequest]

class CoevolutionProcessContext(NamedTuple):
    manager: CbrExtraInteractiveManager
    handler: CoevolultionHandler

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

        self.__alignment_model = CoevolutionOverviewModel(clustal)
        self.__ui.alignmentTable.setModel(self.__alignment_model)
        self.__ui.alignmentTable.doubleClicked.connect(self.__on_position_selected)

        self.__results_model = CoevolutionResultTable(clustal)
        self.__ui.detailsTable.setModel(self.__results_model)

        self.__working_directory = context.create_temporary_directory()

        self.__coevolution_manager: Optional[CoevolutionProcessContext] = None

        self.__set_busy(False)

    def __set_busy(self, busy: bool):

        self.__ui.busyProgress.setVisible(not busy)
        self.__ui.alignmentTable.setEnabled(not busy)
        self.__ui.detailsTable.setEnabled(not busy)

    @pyqtSlot(QModelIndex)
    def __on_position_selected(self, index: QModelIndex):

        manager = self.__coevolution_manager

        if manager is None:
            return
        
        manager.handler.send_message(
            InteractiveRequest(
                query=Query(
                    
                )
            )
        )

    def __run_interactive_process(self, msa: MultipleSeqAlignment):

        previous = self.__coevolution_manager

        if previous is not None:
            previous.manager.stop()

        msa_tmp_file = path.join(self.__working_directory, self.MSA_FILENAME)

        AlignIO.write(msa, msa_tmp_file, 'fasta')
        manager = run_interactive([
            "coevolution",
            "interactive",
            "--input-msa", msa_tmp_file
        ])

        handler = CbrExtraInteractiveHandler(
            manager,
            interactive_response_parser,
            interactive_request_serializer,
            MessageQueingPolicy.HANDLE_ALL
        )

        handler.observe_values() \
            .for_each(
                action=self.__on_result,
                on_error=self.__on_error
            )

        self.__coevolution_manager = CoevolutionProcessContext(manager, handler)

    def __on_result(self, value: InteractiveResponse):
        pass

    def __on_error(self, error: Exception):
        pass

    @pyqtSlot(object)
    def __on_msa_selected(self, msa: MultipleSeqAlignment):

        self.__alignment_model.set_alignment(msa)
        self.__results_model.set_alignment(msa)
        self.__run_interactive_process(msa)
