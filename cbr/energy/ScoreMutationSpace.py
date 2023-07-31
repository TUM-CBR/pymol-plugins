from io import StringIO
import json
import os
from os import path
import pymol
from pymol.wizard.mutagenesis import Mutagenesis, obj_name
import tempfile
from typing import Dict, List, NamedTuple, Optional, TextIO, Tuple

from ..core import sequence
from ..gromacs import ForceField, gmx_box, gmx_configure_emin, gmx_emin_script, gmx_neutralize, gmx_solvate, gmx_topology

from .Foundations import F_METADATA, Mutation, MutationContextBase
from . import ProThermDB

RUN_ITEM_SCRIPT = \
"""
if ! [[ -f {output_log} ]]; then
    srun -D $PWD/{directory} --export GMX_BIN=$GMX_BIN /bin/sh $PWD/{directory}/run_em.sh
#    sbatch -D $PWD --time 02:00:00 --export GMX_BIN=$GMX_BIN $PWD/run.sh
fi

"""

def hash_mutations(mutations : List[Tuple[str, str]]) -> int:
    return hash("".join([str(hash(m)) for m in mutations]))


class MutationContext(MutationContextBase):
    runs : List['MutationContextRun']
    force_field : ForceField

    @property
    def structure_selection(self):
        return 'model %s' % self.name

class MutationContextRun(NamedTuple):
    name : str
    context : 'MutationContext'

    @property
    def force_field(self):
        return self.context.force_field

    @property
    def directory(self):
        return path.join(self.context.directory, self.name)

    @property
    def mutations(self):
        return self.context.mutations

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

    def __apply_muation(self, ctx : MutationContext, reps = 5) -> MutationContext:
        try:
            # Create a copy of the structure
            pymol.cmd.copy(ctx.name, self.__structure)
            sele_name = ctx.name + "_mut"

            #for mutation in ctx.mutations:
            mutation = next(ctx.mutations.__iter__(), None)
            runs = []

            if mutation:

                rot_state = 1
                done = False
                while(not done):
                    run = MutationContextRun(
                        name="rot_%i" % rot_state,
                        context = ctx
                    )
                    self.__wizzard.cleanup()
                    pymol.cmd.select(
                        sele_name,
                        "%s and %s" % (ctx.structure_selection, mutation.selection)
                    )
                    self.__wizzard.do_select(sele_name)
                    self.__wizzard.set_mode(mutation.mutation)
                    self.__wizzard.set_hyd("none")
                    total_states = pymol.cmd.count_states(obj_name)
                    self.__wizzard.do_state(rot_state)
                    self.__wizzard.apply()
                    pymol.cmd.save(
                        run.structure_file,
                        ctx.structure_selection
                    )
                    runs.append(run)
                    rot_state = rot_state + 1
                    done = rot_state <= total_states
            else:
                runs.append(MutationContextRun(name="run_1", context = ctx))

            return ctx._replace(runs = ctx.runs + runs)

        finally:
            pass
            pymol.cmd.delete(ctx.name)

    #def __repackage_structure(self, ctx : MutationContext):
    #    pack_structure(
    #        ctx.structure_file,
    #        ctx.packed_structure_file
    #    )

    def __create_run(self, context : MutationContextRun, run_all : TextIO):
        force_field = context.force_field
        gmx_topology(context.directory, context.structure_file, context.topology_structure_file, context.topology_file, force_field)
        gmx_box(context.directory, context.topology_structure_file, context.box_file)
        gmx_solvate(context.directory, context.box_file, context.solvated_file, context.topology_file)
        gmx_neutralize(context.directory, context.solvated_file, context.neutralized_file, context.topology_file)
        gmx_configure_emin(context.directory, context.neutralized_file, context.emin_config_file, context.topology_file, force_field)
        gmx_emin_script(context.emin_config_file, context.md_script_emin, context.energy_file, context.emin_log_file)
        directory = path.basename(context.directory)
        output_log = path.join("$PWD", directory, path.basename(context.emin_log_file))
        run_all.write(
            RUN_ITEM_SCRIPT.format(
                directory=directory,
                output_log=output_log
            )
        )

    def __score_mutation(self, mutations : List[Tuple[str, str]], force_field : ForceField) -> MutationContext:

        suffix = str(hash_mutations(mutations))[-8:]
        name = "%s_%s" % (self.__structure, suffix)
        directory = path.join(self.__working_directory.name, name)

        if not path.exists(directory):
            os.mkdir(directory)

        context = MutationContext(
            name = name,
            directory = directory,
            mutations = [Mutation(*m) for m in mutations]
        )

        context = self.__apply_muation(context)
        # self.__repackage_structure(context)

        with open(path.join(directory, F_METADATA), 'w') as metadata:
            json.dump(context.to_json_dict(), metadata)

        return context

    def scoring_test(self, mutations : List[Tuple[str, str]], force_field : ForceField):
        self.__score_mutation(mutations, force_field)
        return self.__working_directory

    def score_all(self, force_field : ForceField, m_candidates : Optional[dict] = None):

        candidates = m_candidates or CANDIDATES
        run_all = StringIO()
        run_all.write("#!/bin/sh\n\n")
        done = set()

        def write_run_item(run_args):
            nonlocal run_all
            context = self.__score_mutation(run_args, force_field)

            for run in context.runs:
                self.__create_run(run, run_all)


        for (src, mutations) in candidates.items():
            for mutation in mutations:
                space = [(src, mutation)]
                id = hash_mutations(space)

                if id not in done:
                    write_run_item([(src, mutation)])
                    done.add(id)

        write_run_item([])

        with open(path.join(self.__working_directory.name, "run.sh"), 'w') as run_script:
            run_all.seek(0)
            run_script.write(run_all.read())

        return self.__working_directory

    def score_pro_therm_db(self, db: str, force_field : ForceField):
        entries = ProThermDB.parse(db)
        candidates = {}
        for entry in entries:
            mutation = entry.mutation
            if entry.pdb_code.lower() != self.__structure.lower() \
                or not mutation:
                continue

            key = "resi %i" % mutation.position

            if key not in candidates:
                candidates[key] = []
            mutations = candidates[key]
            mutations.append(sequence.aa_from_letter(mutation.new_residue))

        self.score_all(force_field, m_candidates = candidates)


CANDIDATES = {
    "resi 40": [
        "GLU"
    ],
    "resi 67": [
        "VAL"
    ],
    "resi 89": [
        "ASP"
    ],
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
    ],
    "resi 465": [
        "ILE"
    ],
    "resi 469": [
        "MET"
    ]
}