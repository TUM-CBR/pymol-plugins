from PyQt5.QtWidgets import QDialog, QWidget
from typing import Generic, List, Optional, Type

from .NamedTupleEditor import FieldOrientation, MetaFieldOverridesDict, TTuple, namedtuple_eidtor
from .Ui_NamedtupleEditorDialog import Ui_NamedtupleEditorDialog

class NamedtupleEditorDialog(QDialog, Generic[TTuple]):

    def __init__(
        self,
        *values: Optional[TTuple],
        tuple_type : Optional[Type[TTuple]] = None,
        tuple_field_overrides: Optional[MetaFieldOverridesDict] = None,
        default_item : Optional[TTuple] = None,
        fields_orientation : FieldOrientation = FieldOrientation.Vertical,
        parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(parent)

        self.__ui = Ui_NamedtupleEditorDialog()
        self.__ui.setupUi(self)

        self.__editor_model = namedtuple_eidtor(
            self.__ui.settingsTable,
            *values,
            tuple_type=tuple_type,
            tuple_field_overrides=tuple_field_overrides,
            default_item=default_item,
            fields_orientation=fields_orientation
        )
        self.__previous_value: Optional[List[Optional[TTuple]]] = None

    def exec(self) -> int:

        self.__previous_value = list(self.__editor_model.iterate_values())

        return super().exec()
    
    def __getitem__(self, ix: int) -> Optional[TTuple]:
        return self.__editor_model[ix]
    
    def reject(self) -> None:

        previous = self.__previous_value

        if previous is not None:
            for i,value in enumerate(previous):
                self.__editor_model[i] = value

        return super().reject()
    
def namedtuple_dialog(
    *values: Optional[TTuple],
    tuple_type : Optional[Type[TTuple]] = None,
    tuple_field_overrides: Optional[MetaFieldOverridesDict] = None,
    default_item : Optional[TTuple] = None,
    fields_orientation : FieldOrientation = FieldOrientation.Vertical,
    parent: Optional[QWidget] = None
):
    return NamedtupleEditorDialog(
        *values,
        tuple_type=tuple_type,
        tuple_field_overrides=tuple_field_overrides,
        default_item=default_item,
        fields_orientation=fields_orientation,
        parent=parent
    )
