from typing import Callable, Optional
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt, pyqtSlot
from PyQt5.QtGui import QColor, QDoubleValidator
from PyQt5.QtWidgets import QAbstractItemDelegate, QComboBox, QWidget
from typing import Any, cast, TypeVar

from ...core import color
from ..data import *
from ..operations import query_organisms

from .Ui_CascadesViewer import Ui_CascadesViewer

TOut = TypeVar('TOut')

identity_color_spread = color.color_spread(
    (100,0,0),
    (0,100,0)
)

white = QColor(255,255,255)

class OrganismsTableModel(QAbstractTableModel):

    POLICY_ROW_IX = 0
    META_ROWS = [POLICY_ROW_IX]

    ORGANISM_NAME_COL_IX = 0
    META_COLS = [ORGANISM_NAME_COL_IX]
    META_COLS_HEADERS = ["Organism"]

    def __init__(
        self,
        results : QueryCascadeResult
    ) -> None:
        super().__init__()
        self.__organism_results = results
        self.__step_policies = {}

    def __get_policy(self, step_id : int):
        step = self.steps[step_id]
        return self.__step_policies.get(step.step_id) or QueryStepPolicy.any

    @property
    def __organisms(self):
        return self.__organism_results.organisms

    def enumerate_organisms_rows(self):
        organisms = self.__organisms
        offset = len(self.META_ROWS)
        return zip(
            range(offset, offset + len(organisms)),
            organisms
        )

    @property
    def steps(self) -> List[QueryCascadeResultStepEntry]:
        if len(self.__organisms) > 0:
            return self.__organisms[0].steps
        else:
            return []

    def __get_headers(self):
        return self.META_COLS_HEADERS + [step.step_name for step in self.steps]

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            headers = self.__get_headers()
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)

    def columnCount(self, parent = None) -> int:
        return len(self.META_COLS) + len(self.steps)

    def rowCount(self, parent = None) -> int:
        return len(self.__organisms) + len(self.META_ROWS)

    def __get_organism(self, index: QModelIndex) -> QueryCascadeResultOrganismEntry:
        ix = index.row() - len(self.META_ROWS)
        return self.__organisms[ix]

    def __get_organism_record(
        self,
        col_ix : int,
        entry : QueryCascadeResultOrganismEntry,
        with_organism : Callable[[OrganismResultEntry], TOut],
        with_step : Callable[[QueryCascadeResultStepEntry], TOut]
    ) -> TOut:
        if col_ix == 0:
            return with_organism(entry.organism)
        else:
            step = entry.steps[col_ix - 1]
            return with_step(step)

    def __render_organism(
        self,
        col_ix : int,
        entry : QueryCascadeResultOrganismEntry
    ) -> str:
        return self.__get_organism_record(
            col_ix,
            entry,
            with_organism=lambda o: f"{o.name} (taxid:{o.tax_id})",
            with_step=lambda s: f"{round(s.identity*100,2)}%"
        )

    def __display_role_data(self, index : QModelIndex) -> Any:
        row_ix = index.row()
        col_ix = index.column()

        if row_ix == self.POLICY_ROW_IX:

            meta_cols_count = len(self.META_COLS)
            if col_ix < meta_cols_count:
                return ""
            
            return self.__get_policy(col_ix - meta_cols_count)

        return self.__render_organism(col_ix, self.__get_organism(index))

    def __color_role_data(self, index : QModelIndex) -> Any:

        if index.row() < len(self.META_ROWS) or index.column() < len(self.META_COLS):
            return None

        organism = self.__get_organism(index)
        return self.__get_organism_record(
            index.column(),
            organism,
            with_organism=lambda _: None,
            with_step=lambda step: QColor(*identity_color_spread.get_color(step.identity))
        )

    def flags(self, index : QModelIndex) -> Qt.ItemFlags:
        flags = super().flags(index)

        if index.row() == self.POLICY_ROW_IX and index.column() > 1:
            return flags | Qt.ItemIsEditable

        return flags

    def data(self, index: QModelIndex, role: int = Qt.DisplayRole) -> Any:

        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            return self.__display_role_data(index)
        elif role == Qt.BackgroundColorRole:
            return self.__color_role_data(index)
        elif role == Qt.TextColorRole:
            return self.__color_role_data(index) and white

        return None

class CascadesViewer(QWidget):

    def __init__(
        self,
        db_path : str,
        create_cascade_args : Optional[CreateCascadeDatabaseArgs] = None
    ):
        super().__init__()

        self.__ui = Ui_CascadesViewer()
        self.__ui.setupUi(self)
        self.__model = OrganismsTableModel(QueryCascadeResult(organisms=[]))

        self.__ui.equivalenceLineEdit.setValidator(QDoubleValidator(0,1,2))
        self.__ui.equivalenceLineEdit.setText("0.8")

        self.__db_path = db_path
        self.__create_initial_table()

        self.__policy_combos = {}

    def __create_initial_table(self):
        organisms = query_organisms(self.__db_path)
        self.__set_organisms(organisms)

    def __set_organisms(self, organisms : QueryCascadeResult):
        self.__model = model = OrganismsTableModel(organisms)
        self.__ui.organismsTable.setModel(model)
        self.__ui.organismsTable.resizeColumnsToContents()
        self.__policy_combos = {}

        for i,step in enumerate(model.steps):
            combo = QComboBox()
            combo.setFocusPolicy(Qt.StrongFocus)
            combo.addItems([item.value for item in QueryStepPolicy])
            combo.setCurrentText(cast(str, QueryStepPolicy.any))
            self.__ui.organismsTable.setIndexWidget(
                model.index(0, i + len(model.META_COLS)),
                combo
            )
            self.__policy_combos[step.step_id] = combo