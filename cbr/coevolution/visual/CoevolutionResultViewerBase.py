from abc import abstractmethod
from PyQt5.QtWidgets import QWidget
from typing import List

from ...core.pymol.structure import StructureSelection

class CoevolutionResultViewerBase(QWidget):

    @abstractmethod
    def set_visualization_structure(
        self,
        structure_sele : StructureSelection,
        link : List[int]
    ):
        pass