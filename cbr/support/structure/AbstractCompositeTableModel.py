import pymol
from PyQt5.QtCore import pyqtSignal, QAbstractTableModel, QModelIndex, QObject, Qt
from PyQt5.QtGui import QColor
from typing import Dict, Generic, Iterator, NamedTuple, Optional, Sequence, Tuple, TypeVar, Union

from ...core.color import to_pymol_color
from ...core.pymol.structure import StructureSelection

TModelRecord = TypeVar('TModelRecord')

class ViewHeaderSpec(NamedTuple):
    name: str

DisplayValue = Union[str, float, int]

class ViewRecordAttributes(NamedTuple):
    display: Optional[DisplayValue] = None
    background_color: Optional[QColor] = None

class PymolRecordAttributes(NamedTuple):
    selections: Optional[Sequence[StructureSelection]] = None
    colors: Optional[Sequence[Optional[QColor]]] = None

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

    DEFAULT_PYMOL_ATTRIBUTES = PymolRecordAttributes()

    DEFAULT_QT_ATTRIBUTES = ViewRecordAttributes()

    headers_changed = pyqtSignal()

    def headers(self) -> Sequence[ViewHeaderSpec]:
        raise NotImplementedError("The function headers must have an implementation.")
    
    def pymol_attributes(self, records: Sequence[TModelRecord]) -> Sequence[Sequence[PymolRecordAttributes]]:
        return [
            [self.DEFAULT_PYMOL_ATTRIBUTES for _header in self.headers()]
            for _record in records
        ]
    
    def qt_attributes(self, records: Sequence[TModelRecord]) -> Sequence[Sequence[ViewRecordAttributes]]:
        return [
            [self.DEFAULT_QT_ATTRIBUTES for _header in self.headers()]
            for _record in records
        ]

class AbstractCompositeTableModel(QAbstractTableModel, Generic[TModelRecord]):

    def __init__(self, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)
        self.__views: Sequence[AbstractRecordView[TModelRecord]] = []
        self.__records: Sequence[TModelRecord] = []
        self.__headers: Sequence[ViewHeaderSpec] = []
        self.__orientation: Qt.Orientation = Qt.Orientation.Vertical
        self.__qt_attributes: Dict[QModelIndex, ViewRecordAttributes] = {}
        self.__structures_state: Dict[str, PymolStructureState] = {}

        self.modelReset.connect(self.__on_model_reset)

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
    
    def __on_model_reset(self):
        self.__prepare()

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
        
        
    def __enumerate_attributes(self, records: Sequence[TModelRecord]) -> Iterator[Tuple[QModelIndex, ViewRecordAttributes, PymolRecordAttributes]]:
        """
        This function pre-computes all the data related to the rendering of the
        records in the table as well as any changes that will happen in the
        pymol user interface.
        """
        
        i_header_offset = 0
        for view in self.__views:
            pymol_attributes = view.pymol_attributes(records)
            qt_attributes = view.qt_attributes(records)
            headers_len = len(view.headers())

            if len(pymol_attributes) != len(records) or len(qt_attributes) != len(records):
                raise Exception("There must be an attribute for every record!")

            for i_record,(qt_attribute_list, pymol_attribute_list) in enumerate(zip(qt_attributes, pymol_attributes)):

                if len(qt_attribute_list) != headers_len or len(pymol_attribute_list) != headers_len:
                    raise Exception("There must be an attribute for every header")

                for i_header,(qt_attributes, pymol_attributes) in enumerate(zip(qt_attribute_list, pymol_attribute_list)):

                    if self.__orientation == Qt.Orientation.Vertical:
                        index = self.index(i_record, i_header + i_header_offset)
                    else:
                        index = self.index(i_header + i_header_offset, i_record)

                    yield (index, qt_attributes, pymol_attributes)

            i_header_offset += headers_len
    
    def __prepare(self):
        self.__views = self.__views__()
        self.__records = self.__records__()
        self.__headers = [
            header
            for view in self.__views
            for header in view.headers()
        ]
        self.__orientation = self.__orientation__()
        self.__qt_attributes = qt_dict = {}

        for state in self.__structures_state.values():
            state.restore()

        self.__structures_state = {}

        for (index, qt_attributes, pymol_attribtues) in self.__enumerate_attributes(self.__records):
            qt_dict[index] = qt_attributes
            pymol_selections = pymol_attribtues.selections
            pymol_colors = pymol_attribtues.colors

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