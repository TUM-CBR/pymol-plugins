import itertools
from pymol import cmd
from PyQt5.QtCore import QObject, pyqtSlot
from PyQt5.QtWidgets import QTableWidget
from typing import Dict, List, Optional

class PymolResidueSelector(QObject):

    __counter = itertools.count()

    def __init__(
        self,
        qtable : QTableWidget,
        selection : Optional[str] = None,
        selection_name : Optional[str] = None,
        *args,
        **kwargs
    ):
        super(PymolResidueSelector, self).__init__(*args, **kwargs)
        self.__default_selection_name = "TableSelection_%i" % next(PymolResidueSelector.__counter)
        self.__selection : Optional[str] = selection
        self.__mappings : Dict[int, List[int]] = {}
        self.__qtable = qtable
        qtable.itemSelectionChanged.connect(self.__on_item_selection_changed)
        self.__selection_name = selection_name

    def reset(self):
        self.__selection = None
        self.__mappings = {}

    def set_selection(self, selection : Optional[str]=None):
        self.__selection = selection

    def set_item_riesidues(self, item : int, residue : List[int]):
        self.__mappings[item] = residue

    @pyqtSlot()
    def __on_item_selection_changed(self):
        if self.__selection is None:
            return

        selection = self.__qtable.selectedIndexes()
        residues = [
            "resi %i" % residue
            for index in selection
            for residues in [self.__mappings.get(index.row())] if residues
            for residue in residues
        ]

        selection = self.__selection_name or self.__default_selection_name

        if(len(residues) == 0):
            cmd.select(
                selection
            )
        else:
            selection_str = " & ".join([
                self.__selection,
                " | ".join(residues)
            ])
            cmd.select(
                selection,
                selection=selection_str
            )
