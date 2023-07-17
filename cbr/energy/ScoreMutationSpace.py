import json
import os
from os import path
import pymol
from pymol.wizard.mutagenesis import Mutagenesis
import tempfile
from typing import Dict, List, NamedTuple, Tuple

from ..gromacs import gmx_box, gmx_configure_emin, gmx_emin_script, gmx_neutralize, gmx_solvate, gmx_topology
from ..packmol import pack_structure

batch_all_script = \
    """
    ls | xargs 
    """

class Mutation(NamedTuple):
    selection : str
    mutation : str

    def to_json_dict(self):
        return {
            'selection' : self.selection,
            'mutation' : self.mutation
        }

class MutationContext(NamedTuple):

    name : str
    directory : str
    mutations : List[Mutation]

    def to_json_dict(self):
        return {
            'name' : self.name,
            'directory' : self.directory,
            'mutations': [m.to_json_dict() for m in self.mutations]
        }

    @property
    def structure_selection(self):
        return 'model %s' % self.name

    def mutation(self, i : int) -> str:
        return self.mutations[i].mutation

    @property
    def __pdb_base_filename(self) -> str:
        return "%s.pdb" % self.name

    def __full_path(self, file_name):
        return path.join(self.directory, file_name)

    @property
    def structure_file(self) -> str:
        return self.__full_path(self.__pdb_base_filename)

    @property
    def packed_structure_file(self) -> str:
        return self.__full_path("packed.%s" % self.__pdb_base_filename)

    @property
    def __gmx_base_filename(self) -> str:
        return "%s.gro" % self.name

    @property
    def topology_structure_file(self) -> str:
        return self.__full_path("topology.%s" % self.__gmx_base_filename)

    @property
    def topology_file(self) -> str:
        return self.__full_path("%s.top" % self.name)

    @property
    def box_file(self) -> str:
        return self.__full_path("boxed.%s" % self.__gmx_base_filename)

    @property
    def solvated_file(self) -> str:
        return self.__full_path("solvated.%s" % self.__gmx_base_filename)

    @property
    def neutralized_file(self) -> str:
        return self.__full_path("neutralized.%s" % self.__gmx_base_filename)

    @property
    def emin_config_file(self) -> str:
        return self.__full_path("emin.%s.tpr" % self.name)

    @property
    def emin_log_file(self) -> str:
        return self.__full_path("emin.%s.log" % self.name)

    @property
    def energy_file(self) -> str:
        return self.__full_path("%s.edr" % self.name)

    @property
    def md_script_emin(self) -> str:
        return self.__full_path("run_em.sh")

class ScoreMutationSpace():

    def __init__(
        self,
        structure : str,
        mutations : Dict[str, List[str]]
    ):

        self.__structure = structure
        self.__mutations = mutations
        self.__working_directory = tempfile.TemporaryDirectory()
        self.__wizzard = Mutagenesis()

    def __del__(self):
        self.__working_directory.cleanup()

    def __apply_muation(self, ctx : MutationContext):
        try:
            # Create a copy of the structure
            pymol.cmd.copy(ctx.name, self.__structure)
            sele_name = ctx.name + "_mut"

            for mutation in ctx.mutations:
                self.__wizzard.cleanup()
                pymol.cmd.select(
                    sele_name,
                    "%s and %s" % (ctx.structure_selection, mutation.selection)
                )
                self.__wizzard.do_select(sele_name)
                self.__wizzard.set_mode(mutation.mutation)
                self.__wizzard.apply()
            pymol.cmd.save(
                ctx.structure_file,
                ctx.structure_selection
            )
        finally:
            pass
            pymol.cmd.delete(ctx.name)

    def __repackage_structure(self, ctx : MutationContext):
        pack_structure(
            ctx.structure_file,
            ctx.packed_structure_file
        )

    def __score_mutation(self, mutations : List[Tuple[str, str]]):

        suffix = str(hash("".join([str(hash(m)) for m in mutations])))[-8:]
        name = "%s_%s" % (self.__structure, suffix)
        directory = path.join(self.__working_directory.name, name)

        if not path.exists(directory):
            os.mkdir(directory)

        context = MutationContext(
            name = name,
            directory = directory,
            mutations = [Mutation(*m) for m in mutations]
        )

        self.__apply_muation(context)
        # self.__repackage_structure(context)
        gmx_topology(context.directory, context.structure_file, context.topology_structure_file, context.topology_file)
        gmx_box(context.directory, context.topology_structure_file, context.box_file)
        gmx_solvate(context.directory, context.box_file, context.solvated_file, context.topology_file)
        gmx_neutralize(context.directory, context.solvated_file, context.neutralized_file, context.topology_file)
        gmx_configure_emin(context.directory, context.neutralized_file, context.emin_config_file, context.topology_file)
        gmx_emin_script(context.emin_config_file, context.md_script_emin, context.energy_file, context.emin_log_file)

        with open(path.join(directory, 'metadata.json'), 'w') as metadata:
            json.dump(context.to_json_dict(), metadata)

    def write_batch_script(self):
        pass

    def scoring_test(self, mutations : List[Tuple[str, str]]):
        self.__score_mutation(mutations)
        return self.__working_directory

    def score_all(self):
        for (src, mutations) in CANDIDATES.items():
            for mutation in mutations:
                self.__score_mutation([(src, mutation)])

        self.__score_mutation([])
        return self.__working_directory

CANDIDATES = {
    "resi 168": [
        "ARG",
        "GLN",
        "THR"
    ],
    "resi 297": [
        "GLN",
        "GLU"
    ],
    "resi 329": [
        "LYS",
        "GLN",
        "GLU",
        "ASP"
    ],
    "resi 390": [
        "GLN",
        "ASP",
        "GLU",
        "THR"
    ]
}