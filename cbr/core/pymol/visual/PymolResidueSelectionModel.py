from typing import Callable, Optional
import pymol
from PyQt5.QtCore import Qt
from PyQt5.QtCore import QItemSelection, QItemSelectionModel, pyqtSlot
from PyQt5.QtWidgets import QTableView

from ....control import viter

class PymolResidueSelectionModel(QItemSelectionModel):

    def __init__(
        self,
        view : QTableView,
        get_selection_scope: Callable[[], Optional[str]],
        resi_role : Qt.ItemDataRole,
        selection_name : str = "sele"
    ):
        super().__init__()

        self.__selection_name = selection_name
        self.__role = resi_role
        self.__get_selection_scope = get_selection_scope

        view.setSelectionModel(self)
        self.selectionChanged.connect(self.__on_selection_changed)


    @pyqtSlot(QItemSelection,QItemSelection)
    def __on_selection_changed(self, new: QItemSelection, old: QItemSelection):

        if self.__get_selection_scope() is None:
            return

        pymol.cmd.select(self.__selection_name)

        resi = [
            f"resi {position}"
            for item_range in new
            for item in item_range.indexes()
            for position in viter(item.data(self.__role))
        ]

        selection = " or ".join(resi)
        pymol.cmd.select(
            self.__selection_name,
            f"{self.__get_selection_scope()} and ({selection})"
        )

pymol_residue_selection_model = PymolResidueSelectionModel