import os
from os import path
import pymol
import tempfile
from typing import Dict, List, NamedTuple

from ..gromacs import gmx_box, gmx_configure_emin, gmx_emin_script, gmx_neutralize, gmx_solvate, gmx_topology
from ..packmol import pack_structure

class MutationContext(NamedTuple):

    name : str
    directory : str
    selection : str
    mutation : str

    @property
    def model_copy_selection(self) -> str:
        return 'model %s and %s' % (self.name, self.selection)

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
        return self.__full_path("emin.%s.tpr")

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

    def __del__(self):
        self.__working_directory.cleanup()

    def __apply_muation(self, ctx : MutationContext):
        try:
            # Create a copy of the structure
            pymol.cmd.copy(ctx.name, self.__structure)
            pymol.cmd.alter(
                ctx.model_copy_selection,
                ctx.mutation
            )
            pymol.cmd.rebuild()
            pymol.cmd.save(
                ctx.directory,
                ctx.model_copy_selection
            )
        finally:
            pymol.cmd.delete(ctx.name)

    def __repackage_structure(self, ctx : MutationContext):
        pack_structure(
            ctx.structure_file,
            ctx.packed_structure_file
        )

    def __score_mutation(self, selection : str, mutation : str):

        name = "%s_%s" % (selection[0:8], str(hash(selection))[-4:])
        directory = path.join(self.__working_directory.name, name)

        if not path.exists(directory):
            os.mkdir(directory)

        context = MutationContext(
            name = name,
            directory = directory,
            selection = selection,
            mutation = mutation
        )

        self.__apply_muation(context)
        self.__repackage_structure(context)
        gmx_topology(context.packed_structure_file, context.topology_file, context.topology_file)
        gmx_box(context.topology_structure_file, context.box_file)
        gmx_solvate(context.box_file, context.solvated_file, context.topology_file)
        gmx_neutralize(context.solvated_file, context.neutralized_file, context.topology_file)
        gmx_configure_emin(context.neutralized_file, context.emin_config_file, context.topology_file)
        gmx_emin_script(context.emin_config_file, context.md_script_emin, context.energy_file)

    def scoring_test(self, selection, mutation):
        self.__score_mutation(selection, mutation)
        return self.__working_directory
