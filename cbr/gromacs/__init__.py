from os import path
import pathlib
import tempfile

from ..core.process import simple_execute

def gmx_executable():
    return "gmx_mpi"

def gmx_topology(input_pdb : str, output_file : str, topology : str):

    #!gmx pdb2gmx -f 1fjs_protein.pdb -o 1fjs_processed.gro -water tip3p -ff "charmm27"
    simple_execute(
        [ gmx_executable()
        , "-f"
        , input_pdb
        , "-o"
        , output_file
        , "-water"
        , "tip3p"
        , "-ff"
        , "amber99sb-ildn"
        , "-p"
        , topology
    ])

def gmx_box(input_gr : str, output_gr : str):
    # gmx editconf -f 1fjs_processed.gro -o 1fjs_newbox.gro -c -d 1.0 -bt dodecahedron
    simple_execute(
        [ gmx_executable()
        , "editconf"
        , "-f", input_gr
        , "-o", output_gr
        , "-c"
        , "-d", "1.0"
        , "-bt", "cube"
        , "-box", "10.0" #todo: we need to be smarter here
        ]
    )

def gmx_solvate(input_gr : str, output_gr : str, topl : str):
    # !gmx solvate -cp 1fjs_newbox.gro -cs spc216.gro -o 1fjs_solv.gro -p topol.top
    simple_execute(
        [ gmx_executable()
        , "solvate"
        , "-cp", input_gr
        , "-cs", "spc216.gro"
        , "-o", output_gr
        , "-p", topl
        ]
    )

def gmx_neutralize(input_gr : str, output_gr : str, topl : str):
    # !gmx grompp -f ions.mdp -c 1fjs_solv.gro -p topol.top -o ions.tpr

    with tempfile.TemporaryDirectory() as workdir:
        ions = path.join(workdir, "ions.mdp")
        ions_tpr = path.join(workdir, "ions.tpr")
        pathlib.Path(ions).touch()
        simple_execute(
            [ gmx_executable()
            , "grompp"
            , "-f", ions
            , "-c", input_gr
            , "-p", topl
            , "-o", ions_tpr
            ]
        )

        # !printf "SOL\n" | gmx genion -s ions.tpr
        # -o 1fjs_solv_ions.gro -conc 0.15 -p topol.top -pname NA -nname CL -neutral
        simple_execute(
            [ gmx_executable()
            , "genion"
            , "-s", ions_tpr
            , "-conc", "0.15"
            , "-p", topl
            , "-pname", "NA"
            , "-nname", "CL"
            , "-neutral"
            , "-o", output_gr
            ],
            "SOL"
        )

amber_emin_cfg = \
    """
    title       = Amber steepest descent enrgy minimisation

    ; Parameters describing what to do, when to stop and what to save
    integrator  = steep  ; Algorithm (steep = steepest descent minimization)
    emtol       = 1000.0 ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
    emstep      = 0.01   ; Minimization step size
    nstenergy   = 500    ; save energies every 1.0 ps, so we can observe if we are successful
    nsteps      = -1     ; run as long as we need
    ; Settings that make sure we run with parameters in harmony with the selected force-field
    constraints             = h-bonds   ; bonds involving H are constrained
    rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
    rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
    vdw-modifier            = Potential-shift ; Amber specific
    DispCorr                = EnerPres  ; account for cut-off vdW scheme
    coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
    fourierspacing          = 0.125     ; grid spacing for FFT
    """

def gmx_configure_emin(input_gr : str, output_gr : str, topl : str, maxwarns = 1):
# !gmx grompp -f input/emin-charmm.mdp -c 1fjs_solv_ions.gro -p topol.top -o em.tpr
    with tempfile.TemporaryDirectory() as workdir:
        mdp_config = path.join(workdir, "emin_amber.mdp")
        with open(mdp_config, 'w') as mdp:
            mdp.write(amber_emin_cfg)

        simple_execute(
            [ gmx_executable()
            , "grompp"
            , "-f", mdp_config
            , "-c", input_gr
            , "-p", topl
            , "-o", output_gr
            , "-maxwarn", str(maxwarns)
            ]
        )

run_em_sh_template = \
    """
    #!/bin/bash
    $GMX_BIN mdrun %s
    """

def gmx_emin_script(input_tpr : str, output_script : str, output_energy : str):
    args = " ".join(
        [ "-s", input_tpr
        , "-e", output_energy
        , "-v"
        ]
    )

    with open(output_script, 'w') as output_script_file:
        output_script_file.write(run_em_sh_template % args)