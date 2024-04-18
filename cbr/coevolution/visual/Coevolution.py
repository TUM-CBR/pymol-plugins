
from ...core.Context import Context
from ...support.structure.AbstractStructureTableModel import AbstractStructureTableModel
from .Ui_coevolution import Ui_Coevolution

class AlignmentTableModel(AbstractStructureTableModel):
    pass

class Coevolution(Ui_Coevolution):

    def __init__(
        self,
        context: Context
    ) -> None:
        super().__init__()