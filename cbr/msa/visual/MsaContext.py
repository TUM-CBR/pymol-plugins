from abc import abstractproperty
from typing import Optional

from ...core.pymol import structure
from ...core.pymol.structure import StructureSelection
from ...support.msa import Msa


class MsaContext:

    @abstractproperty
    def selected_msa_sequence_name(self) -> str:
        raise NotImplementedError()

    @abstractproperty
    def selected_structure(self) -> Optional[StructureSelection]:
        raise NotImplementedError()

    @abstractproperty
    def sequences(self) -> Msa:
        raise NotImplementedError()

    def get_structure_sequence(self) -> str:

        if self.selected_structure:
            return structure.get_pdb_sequence(self.selected_structure)
        else:
            raise Exception("A structure must be selected to get its sequence!")