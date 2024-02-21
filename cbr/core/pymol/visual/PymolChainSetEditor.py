from typing import Iterable

from .. import objects
from ...Qt.visual.SetEditorModel import SetsEditorModel
from ..structure import StructureSelection

def list_available_selections() -> Iterable[StructureSelection]:

    for (model, chain) in objects.iter_chains():
        yield StructureSelection(
            structure_name=model,
            chain_name=chain,
            segment_identifier=None
        )

class PymolChainSetEditorModel(SetsEditorModel[StructureSelection]):

    def __init__(
        self,
        exclusive: bool = False
    ):
        super().__init__(
            set(list_available_selections()),
            [set() for _ in range(5)],
            exclusive=exclusive,
            render=lambda ms: ms.show()
        )

    def refresh(self):
        self.set_universe(set(list_available_selections()))