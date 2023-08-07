import json
import os
from os import path
import pymol
from pymol.wizard.mutagenesis import Mutagenesis, obj_name
import random
from typing import Iterable, List, NamedTuple, Optional, TextIO

from ..core.pymol.structure import get_selection_sequece, StructureSelection
from ..core.sequence import ensure_abbreviation, ensure_oneletter

from . import ProThermDB

class MutationResult(NamedTuple):
    structure_name : str
    mutated_structure_name : str
    original_residue : str
    mutation_position : int
    new_reside : str

    def to_json_dict(self):
        return {
            'structure_name': self.structure_name,
            'mutated_structure_name': self.mutated_structure_name,
            'original_residue': self.original_residue,
            'mutation_position': self.mutation_position,
            'new_reside': self.new_reside
        }

    @staticmethod
    def from_json_dict(json_dict : dict) -> 'MutationResult':
        return MutationResult(**json_dict)

    @staticmethod
    def from_stream(stream : TextIO) -> 'MutationResult':
        return MutationResult.from_json_dict(json.load(stream))

    def get_energy_log_file(self, directory : str):
        return path.join(directory, "%s.log" % self.structure_name)

class Mutation(NamedTuple):
    structure_position : int
    new_residue : str

class MutationSpace(NamedTuple):
    mutations : Iterable[Mutation]
    structure_selection : StructureSelection
    max_rotamers_per_mutation : Optional[int]

class MutationContext(NamedTuple):
    results_folder : str
    wizard : Mutagenesis

def apply_mutation(
    context : MutationContext,
    space : MutationSpace,
    mutation : Mutation
) -> Iterable[MutationResult]:

    structure = space.structure_selection
    sequence = get_selection_sequece(structure.selection)
    wizard = context.wizard

    # Strings start at index 0, mutation positions are given relative to index 1
    original_residue = sequence[mutation.structure_position - 1]
    new_reside = ensure_oneletter(mutation.new_residue)
    candidates = None

    while(True):
        rotamer = candidates and candidates.pop() or 1
        mutated_structure_name = \
            "{structure}.{original_residue}{position}{new_residue}.rot{rotamer}".format(
                structure = structure.structure_name,
                original_residue = original_residue,
                new_residue = new_reside,
                position = mutation.structure_position,
                rotamer = rotamer
            )

        mutated_structure_selection = structure._replace(
            structure_name = mutated_structure_name
        )
        max_rotamers_per_mutation = space.max_rotamers_per_mutation or 999

        try:
            selection_name = mutated_structure_name + ".selection"
            wizard.cleanup()
            pymol.cmd.copy(mutated_structure_name, structure.structure_name)
            pymol.cmd.select(
                selection_name,
                "%s & resi %i" % (mutated_structure_selection.selection, mutation.structure_position)
            )
            wizard.do_select(selection_name)
            wizard.set_mode(ensure_abbreviation(mutation.new_residue))
            wizard.set_hyd("none")
            rotamer_count = pymol.cmd.count_states(obj_name)

            if candidates is None:
                candidates = \
                    random.sample(
                        range(
                            2,
                            rotamer_count + 1
                        ),
                        min(max_rotamers_per_mutation, rotamer_count) - 1
                    )
            wizard.do_state(rotamer)
            wizard.apply()
            yield MutationResult(
                structure_name = structure.structure_name,
                mutated_structure_name = mutated_structure_name,
                original_residue = original_residue,
                mutation_position = mutation.structure_position,
                new_reside = new_reside
            )

            if not any(candidates):
                return
        except Exception:
            pymol.cmd.delete(mutated_structure_name)
            raise
        
def structure_directory(context : MutationContext, structure_name : str):
    directory = path.join(context.results_folder, structure_name)
    if not os.path.exists(directory):
        os.makedirs(directory)

    return directory

def generate_mutations(
    context : MutationContext,
    space : MutationSpace
) -> Iterable[str]:

    for mutation in space.mutations:
        for mutation_result in apply_mutation(context, space, mutation):

            mutated_structure = mutation_result.mutated_structure_name
            result_directory = structure_directory(context, mutated_structure)
            try:
                result_file = path.join(
                    result_directory,
                    "%s.pdb" % mutated_structure
                )

                metadata_file = path.join(
                    result_directory,
                    "%s.pdb.metadata.json" % mutated_structure
                )

                with open(metadata_file, 'w') as metadata:
                    json.dump(mutation_result.to_json_dict(), metadata)

                pymol.cmd.save(
                    result_file,
                    "model %s" % mutated_structure
                )

                yield result_file
            finally:
                pymol.cmd.delete(mutated_structure)

def run_prothermdb(
    context : MutationContext,
    structure : StructureSelection,
    location : str,
    max_rotamers : Optional[int] = None
) -> List[str]:

    entries = [entry for entry in ProThermDB.parse(location) if entry.pdb_code.upper() == structure.structure_name.upper()]
    mutated = set()
    mutations = [
        Mutation(structure_position=entry.mutation.position, new_residue=entry.mutation.new_residue) \
        for entry in entries if entry.mutation
        for item in [(entry.mutation.position, entry.mutation.new_residue)] if item not in mutated and (mutated.add(item) or True)
    ]

    return list(
        generate_mutations(
            context,
            MutationSpace(
                structure_selection=structure,
                mutations=mutations,
                max_rotamers_per_mutation = max_rotamers
            )
        )
    )