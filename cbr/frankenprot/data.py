from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from pymol import cmd
from typing import Dict, List, NamedTuple, Optional

from ..core.pymol.structure import StructureSelection

class FrankenProt(NamedTuple):
    base_structure: StructureSelection
    fragments_structure: StructureSelection
    replacements: Dict[int, List[int]]
    base_structure_state: int = 1
    fragments_structure_state: int = 1

    def show_preview(self):
        base_resv = list(self.replacements.keys())
        fragment_resv = list(r for resvs in self.replacements.values() for r in resvs)

        base_model = self.base_structure
        base_sel = base_model.selection

        fragment_model = self.fragments_structure
        fragment_sel = fragment_model.selection

        cmd.show(representation='cartoon', selection=base_sel)
        cmd.color("blue", base_sel)

        if len(base_resv) > 0:
            cmd.hide(selection=base_model.residue_selection(base_resv))

        cmd.hide(selection=fragment_sel)

        if len(fragment_resv) > 0:
            start = min(fragment_resv)
            end = max(fragment_resv)
            cmd.show(representation='cartoon', selection=fragment_model.residue_selection([start - 1, end + 1] + fragment_resv))
            cmd.color("green", fragment_sel)

    def get_structure(self, structure_name: str, chain_id: str = "A") -> Structure:

        scope = list(self.base_structure.get_sequence().keys())
        scope.sort()
        sub_positions = list(self.replacements.keys())
        sub_positions.sort()
        residues: List[Residue] = []
        atom_serial: int = 1

        last_resv: Optional[int] = None
        def add_atom(
            resv: int,
            resn: str,
            name: str,
            x: float,
            y: float,
            z: float,
            bfactor: float,
            occupancy: float,
            altloc: str,
            full_name: str,
            element: str
        ):
            nonlocal atom_serial
            nonlocal last_resv
            nonlocal residues
            if len(residues) == 0 \
                or last_resv != resv:
                residue = Residue(
                    (' ', len(residues) + 1, ' '),
                    resn,
                    '    '
                )
                residues.append(residue)
            else:
                residue = residues[-1]

            atom = Atom(
                name,
                [x,y,z],
                bfactor,
                occupancy,
                ' ',
                full_name,
                atom_serial,
                element
            )

            try:
                residue.add(atom)
            except PDBConstructionException:
                # Atoms have alternative conformations,
                # We only need the first
                pass

            atom_serial += 1
            last_resv = resv

        while(len(sub_positions) > 0):

            # get the position of the next substitution
            sub_start = sub_positions.pop(0)

            # Add the atoms of the base structure upto the
            # position to be substituted
            sub_start_ix = scope.index(sub_start)
            positions_to_add = scope[:sub_start_ix]
            scope = scope[sub_start_ix+1:]

            # Add atoms from the base structure up to the residue
            # that will be substituted
            if len(positions_to_add) > 0:
                last_resv: Optional[int] = None
                cmd.iterate_state(
                    self.base_structure_state,
                    self.base_structure.residue_selection(positions_to_add),
                    'add_atom(resv, resn, name, x, y, z, b, q, alt, name, elem)',
                    space={'add_atom': add_atom}
                )

            # add the atoms from the fragments structure to be placed
            # instead of a residue from the base structure
            new_atoms = self.replacements[sub_start]
            if len(new_atoms) > 0:
                last_resv: Optional[int] = None
                cmd.iterate_state(
                    self.fragments_structure_state,
                    self.fragments_structure.residue_selection(new_atoms),
                    'add_atom(resv, resn, name, x, y, z, b, q, alt, name, elem)',
                    space={'add_atom': add_atom}
                )

        # Add the remaining atoms
        if len(scope) > 0:
            last_resv: Optional[int] = None
            cmd.iterate_state(
                self.base_structure_state,
                self.base_structure.residue_selection(scope),
                'add_atom(resv, resn, name, x, y, z, b, q, alt, name, elem)',
                space={'add_atom': add_atom}
            )

        chain = Chain(chain_id)
        for residue in residues:
            chain.add(residue)

        model = Model(structure_name)
        model.add(chain)
        structure = Structure(structure_name)
        structure.add(model)

        return structure
    
    def save_pdb(self, name: str, location: str, chain_id: str = "A"):

        structure = self.get_structure(name, chain_id)
        pdb = PDBIO()
        pdb.set_structure(structure)
        pdb.save(location)
