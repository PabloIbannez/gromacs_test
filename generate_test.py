#!/usr/bin/env python3

from pathlib import Path
import argparse
import textwrap
import subprocess


def fmt(x: float) -> str:
    """Format number for folder names."""
    return f"{x:g}"


def write_file(path: Path, content: str):
    path.write_text(textwrap.dedent(content).lstrip())


def main():
    parser = argparse.ArgumentParser(
        description="Create a GROMACS pure-water box test system."
    )

    parser.add_argument("x", type=float, help="Box size in X direction, in nm")
    parser.add_argument("y", type=float, help="Box size in Y direction, in nm")
    parser.add_argument("z", type=float, help="Box size in Z direction, in nm")

    parser.add_argument(
        "--gmx",
        default="gmx",
        help="Path to GROMACS binary. Default: gmx",
    )

    parser.add_argument(
        "--run",
        action="store_true",
        help="Immediately run make after creating the folder.",
    )

    args = parser.parse_args()

    x, y, z = fmt(args.x), fmt(args.y), fmt(args.z)
    folder = Path(f"water_{x}_{y}_{z}")
    folder.mkdir(exist_ok=True)

    write_file(
        folder / "topol.template.top",
        """
        #include "oplsaa.ff/forcefield.itp"
        #include "oplsaa.ff/spc.itp"

        [ system ]
        Pure water box

        [ molecules ]
        """,
    )

    write_file(
        folder / "em.mdp",
        """
        integrator              = steep
        nsteps                  = 5000
        emtol                   = 100.0
        emstep                  = 0.01

        cutoff-scheme           = Verlet
        nstlist                 = 20

        coulombtype             = PME
        rlist                   = 1.0
        rvdw                    = 1.0
        rcoulomb                = 1.0

        constraints             = h-bonds
        constraint-algorithm    = lincs

        fourier-spacing         = 0.125
        pbc                     = xyz
        """,
    )

    write_file(
        folder / "test.mdp",
        """
        integrator              = md
        nsteps                  = -1
        dt                      = 0.002

        nstfout                 = 0
        nstxout                 = 0
        nstvout                 = 0
        nstxout-compressed      = 0
        nstcalcenergy           = 5000
        nstcomm                 = 5000
        nstenergy               = 0
        nstlog                  = 0

        constraints             = h-bonds
        constraint-algorithm    = lincs

        cutoff-scheme           = Verlet
        nstlist                 = 20

        coulombtype             = PME
        rlist                   = 1.0
        rvdw                    = 1.0
        rcoulomb                = 1.0

        tcoupl                  = v-rescale
        tc-grps                 = System
        tau-t                   = 1.0
        ref-t                   = 300
        nsttcouple              = 100

        fourier-spacing         = 0.125
        pbc                     = xyz
        """,
    )

    makefile = f"""
GMX ?= {args.gmx}

X = {args.x:g}
Y = {args.y:g}
Z = {args.z:g}

MDP = test.mdp
EM_MDP = em.mdp

TOP = topol.top
TEMPLATE_TOP = topol.template.top

CONF = conf.gro
EM_TPR = em.tpr
EM_CONF = em.gro
TPR = topol.tpr

EM_DEFFNM = em
DEFFNM = md

.PHONY: all info solvate minimize prepare run clean distclean

all: run

info:
\t@echo "GMX = $(GMX)"
\t@echo "Box = $(X) $(Y) $(Z) nm"

solvate:
\tcp $(TEMPLATE_TOP) $(TOP)
\t$(GMX) solvate \\
\t\t-cs spc216.gro \\
\t\t-box $(X) $(Y) $(Z) \\
\t\t-o $(CONF) \\
\t\t-p $(TOP)

minimize: solvate
\t$(GMX) grompp \\
\t\t-f $(EM_MDP) \\
\t\t-c $(CONF) \\
\t\t-p $(TOP) \\
\t\t-o $(EM_TPR)
\t$(GMX) mdrun \\
\t\t-s $(EM_TPR) \\
\t\t-deffnm $(EM_DEFFNM) \\
\t\t-nb cpu \\
\t\t-v

prepare: minimize
\t$(GMX) grompp \\
\t\t-f $(MDP) \\
\t\t-c $(EM_CONF) \\
\t\t-p $(TOP) \\
\t\t-o $(TPR)

run: prepare
\t$(GMX) mdrun \\
\t\t-s $(TPR) \\
\t\t-deffnm $(DEFFNM) \\
\t\t-nb gpu \\
\t\t-update gpu \\
\t\t-v

clean:
\trm -f $(CONF) $(TOP)
\trm -f $(EM_TPR) $(EM_CONF) em.log em.edr em.trr em.xtc em.cpt em_prev.cpt
\trm -f $(TPR)
\trm -f md.log md.edr md.trr md.xtc md.gro md.cpt md_prev.cpt

distclean: clean
\trm -f \\#* *.log
"""

    write_file(folder / "Makefile", makefile)

    print(f"Created folder: {folder}")
    print()
    print("To run:")
    print(f"  make -C {folder}")
    print()
    print("Or with another GROMACS binary:")
    print(f"  make -C {folder} GMX=/path/to/gmx")
    print()
    print("Available targets:")
    print(f"  make -C {folder} solvate")
    print(f"  make -C {folder} minimize")
    print(f"  make -C {folder} prepare")
    print(f"  make -C {folder} run")
    print(f"  make -C {folder} clean")

    if args.run:
        subprocess.run(["make", "-C", str(folder)], check=True)


if __name__ == "__main__":
    main()
