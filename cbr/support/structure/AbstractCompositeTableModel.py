import pymol
from PyQt5.QtCore import pyqtSignal, pyqtSlot, QAbstractTableModel, QModelIndex, QObject, Qt
from PyQt5.QtGui import QColor
from typing import Any, Dict, Generic, Iterator, NamedTuple, Optional, Sequence, TypeVar, Union

from ...core.color import to_pymol_color
from ...core.pymol.structure import StructureSelection

TModelRecord = TypeVar('TModelRecord')

class ViewHeaderSpec(NamedTuple):
    name: str

DisplayValue = Union[str, float, int]

class ViewRecordAttributes(NamedTuple):
    display: Optional[DisplayValue] = None
    background_color: Optional[QColor] = None

    @classmethod
    def default(cls) -> 'ViewRecordAttributes':
        return DEFAULT_QT_ATTRIBUTES
    
class PymolRecordAttributes(NamedTuple):
    selections: Optional[Sequence[StructureSelection]] = None
    colors: Optional[Sequence[Optional[QColor]]] = None

    @classmethod
    def default(cls) -> 'PymolRecordAttributes':
        return DEFAULT_PYMOL_ATTRIBUTES

class ViewRecords(NamedTuple):
    headers: Sequence[ViewHeaderSpec]
    qt_attributes: Sequence[Sequence[ViewRecordAttributes]]
    pymol_attributes: Sequence[Sequence[PymolRecordAttributes]]

DEFAULT_PYMOL_ATTRIBUTES = PymolRecordAttributes()

DEFAULT_QT_ATTRIBUTES = ViewRecordAttributes()

class PymolStructureState(NamedTuple):
    structure: StructureSelection
    original_colors: Dict[int, int]

    def restore(self):
        self.structure.set_colors(self.original_colors)

    @classmethod
    def create(cls, model_name: str) -> 'PymolStructureState':
        selection = StructureSelection(model_name, None, None)
        return PymolStructureState(
            selection,
            selection.get_color_indexes()
        )

class AbstractRecordView(QObject, Generic[TModelRecord]):

    DEFAULT_PYMOL_ATTRIBUTES = PymolRecordAttributes.default()

    DEFAULT_QT_ATTRIBUTES = ViewRecordAttributes.default()

    EMPTY_HEADERS: Sequence[ViewHeaderSpec] = []
    EMPTY_VIEW_RECORDS: Sequence[ViewRecordAttributes] = []
    EMPTY_PYMOL_RECORDS: Sequence[PymolRecordAttributes] = []

    content_changed = pyqtSignal(QObject)
    
    def attributes(self, records: Sequence[TModelRecord]) -> ViewRecords:
        return ViewRecords(
            headers=self.EMPTY_HEADERS,
            qt_attributes=[self.EMPTY_VIEW_RECORDS for _record in records],
            pymol_attributes=[self.EMPTY_PYMOL_RECORDS for _record in records]
        )

class IndexEntry(NamedTuple):
    entry_index: QModelIndex
    view_attributes: ViewRecordAttributes
    pymol_attribtues: PymolRecordAttributes
    record: Any

class AbstractCompositeTableModel(QAbstractTableModel, Generic[TModelRecord]):

    def __init__(self, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)
        self.__views: Sequence[AbstractRecordView[TModelRecord]] = []
        self.__records: Dict[QModelIndex, TModelRecord] = {}
        self.__headers: Sequence[ViewHeaderSpec] = []
        self.__orientation: Qt.Orientation = Qt.Orientation.Vertical
        self.__qt_attributes: Dict[QModelIndex, ViewRecordAttributes] = {}
        self.__structures_state: Dict[str, PymolStructureState] = {}

        self.__view_data_version = 1
        self.__current_data_version = 0

    def __records__(self) -> Sequence[TModelRecord]:
        raise NotImplementedError("The function __records__ must have an implementation.")
    
    def __views__(self) -> Sequence[AbstractRecordView[TModelRecord]]:
        raise NotImplementedError("The function __views__ must  habe an implementation.")
    
    def __orientation__(self) -> Qt.Orientation:
        """
        This function describes in what orientation should the views of this table
        be displayed. Vertical means that records will be stacked vertically meaning
        each header of the view will be mapped to a column.
        """

        return Qt.Orientation.Vertical

    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:

        if self.__orientation == Qt.Orientation.Vertical:
            return len(self.__records)
        else:
            return len(self.__headers)
        
    def columnCount(self, parent: Optional[QModelIndex] = None) -> int:
        if self.__orientation == Qt.Orientation.Vertical:
            return len(self.__headers)
        else:
            return len(self.__records)
        
    def get_record(self, index: QModelIndex) -> Optional[TModelRecord]:
        return self.__records.get(index)
        
    def __enumerate_attributes(
        self,
        records: Sequence[TModelRecord],
        views_attributes: Sequence[ViewRecords]
    ) -> Iterator[IndexEntry]:
        """
        This function pre-computes all the data related to the rendering of the
        records in the table as well as any changes that will happen in the
        pymol user interface.
        """
        
        i_header_offset = 0
        for attributes in views_attributes:
            qt_attributes = attributes.qt_attributes
            pymol_attributes = attributes.pymol_attributes

            headers_len = len(attributes.headers)

            if len(pymol_attributes) != len(records) or len(qt_attributes) != len(records):
                raise Exception("There must be an attribute for every record!")

            for i_record,(record, qt_attribute_list, pymol_attribute_list) in enumerate(zip(records, qt_attributes, pymol_attributes)):

                if len(qt_attribute_list) != headers_len or len(pymol_attribute_list) != headers_len:
                    raise Exception("There must be an attribute for every header")

                for i_header,(qt_attributes, pymol_attributes) in enumerate(zip(qt_attribute_list, pymol_attribute_list)):

                    if self.__orientation == Qt.Orientation.Vertical:
                        index = self.index(i_record, i_header + i_header_offset)
                    else:
                        index = self.index(i_header + i_header_offset, i_record)

                    yield IndexEntry(index, qt_attributes, pymol_attributes, record)

            i_header_offset += headers_len

    @pyqtSlot(QObject)
    def __on_content_changed(self, view: AbstractRecordView[TModelRecord]):
        self.__setup()
        self.modelReset.emit()
    
    def __setup(self):

        for view in self.__views:
            view.content_changed.disconnect(self.__on_content_changed)

        self.__views = self.__views__()

        for view in self.__views:
            view.content_changed.connect(self.__on_content_changed)

        self.__get_view_data()

    def __get_view_data(self):
        records = self.__records__()
        views_attributes = [view.attributes(records) for view in self.__views]
        self.__headers = [
            header
            for attributes in views_attributes
            for header in attributes.headers
        ]
        self.__orientation = self.__orientation__()
        self.__records = records_dict = {} 
        self.__qt_attributes = qt_dict = {}

        for state in self.__structures_state.values():
            state.restore()

        self.__structures_state = {}

        for entry in self.__enumerate_attributes(records, views_attributes):
            index = entry.index
            qt_dict[index] = entry.view_attributes
            records_dict[index] = entry.record
            pymol_attributes = entry.pymol_attribtues
            pymol_selections = pymol_attributes.selections
            pymol_colors = pymol_attributes.colors

            if pymol_selections is None or pymol_colors is None:
                continue

            for (selection, color) in zip(pymol_selections, pymol_colors):
                if color is None:
                    continue

                if selection.structure_name not in  self.__structures_state:
                    self.__structures_state[selection.structure_name] = PymolStructureState.create(selection.structure_name)

                pymol.cmd.color(
                    to_pymol_color(color),
                    selection.selection
                )

    def __ensure_data_version(self):

        if self.__current_data_version != self.__view_data_version:
            self.__setup()
            self.__current_data_version = self.__view_data_version

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        self.__ensure_data_version()

        record = self.__qt_attributes.get(index)

        if record is None:
            return None
        elif role == Qt.ItemDataRole.DisplayRole:
            return record.display
        elif role  == Qt.ItemDataRole.BackgroundColorRole:
            return record.background_color
        else:
            return  None