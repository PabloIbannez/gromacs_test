#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import subprocess
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class GmxBuild:
    name: str
    binary: str
    execution: dict[str, Any]


@dataclass(frozen=True)
class ProteinCase:
    protein_name: str
    copy_count: int
    copy_axis: str
    copy_spacing_nm: float
    box_x_nm: float
    box_y_nm: float
    box_z_nm: float

    @property
    def folder_name(self) -> str:
        return f"protein_{self.protein_name}_N{self.copy_count}"


def fmt(value: float) -> str:
    return f"{value:g}"


def write_file(path: Path, content: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(content).lstrip())
    if executable:
        path.chmod(path.stat().st_mode | 0o111)


def require_mapping(config: dict[str, Any], key: str) -> dict[str, Any]:
    value = config.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"Configuration key '{key}' must be an object")
    return value


def require_list(config: dict[str, Any], key: str) -> list[Any]:
    value = config.get(key)
    if not isinstance(value, list):
        raise ValueError(f"Configuration key '{key}' must be a list")
    return value


def load_config(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        config = json.load(handle)
    if not isinstance(config, dict):
        raise ValueError("Top-level configuration must be a JSON object")
    return config


def parse_builds(config: dict[str, Any]) -> list[GmxBuild]:
    builds = []
    for entry in require_list(config, "builds"):
        if not isinstance(entry, dict):
            raise ValueError("Each build entry must be an object")
        execution = entry.get("execution", {})
        if not isinstance(execution, dict):
            raise ValueError("Build entry 'execution' must be an object")
        builds.append(GmxBuild(str(entry["name"]), str(entry["binary"]), execution))
    if not builds:
        raise ValueError("At least one build must be configured")
    return builds


def make_mdp(parameters: dict[str, Any]) -> str:
    return "".join(f"{key:<24}= {value}\n" for key, value in parameters.items())


def make_cases(campaign: dict[str, Any], explicit_copy_count: int | None) -> list[ProteinCase]:
    protein_name = str(campaign.get("protein_name", "protein"))
    copy_axis = str(campaign.get("copy_axis", "x"))
    if copy_axis not in {"x", "y", "z"}:
        raise ValueError("'campaign.copy_axis' must be one of: x, y, z")

    if explicit_copy_count is None:
        copy_counts = [int(value) for value in require_list(campaign, "copy_counts")]
    else:
        copy_counts = [explicit_copy_count]
    if not copy_counts or any(value < 1 for value in copy_counts):
        raise ValueError("Protein copy counts must be positive integers")

    spacing = float(campaign.get("copy_spacing_nm", 4.0))
    margin = float(campaign.get("box_margin_nm", 2.0))
    box_y = float(campaign.get("box_y_nm", 6.0))
    box_z = float(campaign.get("box_z_nm", 6.0))

    cases = []
    for copy_count in copy_counts:
        box_x = float(campaign.get("box_x_nm", copy_count * spacing + margin))
        if copy_axis == "y":
            box_y = float(campaign.get("box_y_nm", copy_count * spacing + margin))
        elif copy_axis == "z":
            box_z = float(campaign.get("box_z_nm", copy_count * spacing + margin))
        cases.append(
            ProteinCase(
                protein_name=protein_name,
                copy_count=copy_count,
                copy_axis=copy_axis,
                copy_spacing_nm=spacing,
                box_x_nm=box_x,
                box_y_nm=box_y,
                box_z_nm=box_z,
            )
        )
    return cases


def make_replicate_gro_script() -> str:
    return r"""
    #!/usr/bin/env python3

    from __future__ import annotations

    import argparse
    from pathlib import Path


    def parse_args() -> argparse.Namespace:
        parser = argparse.ArgumentParser(description="Replicate a GROMACS .gro file along one axis.")
        parser.add_argument("--input", required=True, type=Path)
        parser.add_argument("--output", required=True, type=Path)
        parser.add_argument("--copies", required=True, type=int)
        parser.add_argument("--axis", default="x", choices=("x", "y", "z"))
        parser.add_argument("--spacing", required=True, type=float)
        parser.add_argument("--box", nargs=3, required=True, type=float, metavar=("X", "Y", "Z"))
        return parser.parse_args()


    def parse_atom_line(line: str) -> tuple[int, str, str, int, float, float, float, str]:
        residue_number = int(line[0:5])
        residue_name = line[5:10]
        atom_name = line[10:15]
        atom_number = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        rest = line[44:].rstrip("\n")
        return residue_number, residue_name, atom_name, atom_number, x, y, z, rest


    def format_atom_line(
        residue_number: int,
        residue_name: str,
        atom_name: str,
        atom_number: int,
        x: float,
        y: float,
        z: float,
        rest: str,
    ) -> str:
        return f"{residue_number % 100000:5d}{residue_name}{atom_name}{atom_number % 100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}{rest}\n"


    def main() -> None:
        args = parse_args()
        if args.copies < 1:
            raise SystemExit("--copies must be >= 1")

        lines = args.input.read_text().splitlines()
        title = lines[0] if lines else "Replicated protein"
        atom_count = int(lines[1].strip())
        atom_lines = lines[2 : 2 + atom_count]
        atoms = [parse_atom_line(line) for line in atom_lines]

        axis_index = {"x": 0, "y": 1, "z": 2}[args.axis]
        coordinates = [(atom[4], atom[5], atom[6]) for atom in atoms]
        mins = [min(values[index] for values in coordinates) for index in range(3)]
        maxs = [max(values[index] for values in coordinates) for index in range(3)]
        span = [maxs[index] - mins[index] for index in range(3)]
        origin = [
            max((args.box[0] - (span[0] if axis_index != 0 else span[0] + args.spacing * (args.copies - 1))) / 2.0, 0.0),
            max((args.box[1] - (span[1] if axis_index != 1 else span[1] + args.spacing * (args.copies - 1))) / 2.0, 0.0),
            max((args.box[2] - (span[2] if axis_index != 2 else span[2] + args.spacing * (args.copies - 1))) / 2.0, 0.0),
        ]

        output_lines = [f"{title} x {args.copies}", f"{len(atoms) * args.copies:5d}"]
        atom_number = 1
        residue_offset = 0
        max_residue = max(atom[0] for atom in atoms)
        for copy_index in range(args.copies):
            shift = [origin[index] - mins[index] for index in range(3)]
            shift[axis_index] += args.spacing * copy_index
            for residue_number, residue_name, atom_name, _old_atom_number, x, y, z, rest in atoms:
                output_lines.append(
                    format_atom_line(
                        residue_number + residue_offset,
                        residue_name,
                        atom_name,
                        atom_number,
                        x + shift[0],
                        y + shift[1],
                        z + shift[2],
                        rest,
                    ).rstrip("\n")
                )
                atom_number += 1
            residue_offset += max_residue
        output_lines.append(f"{args.box[0]:10.5f}{args.box[1]:10.5f}{args.box[2]:10.5f}")
        args.output.write_text("\n".join(output_lines) + "\n")


    if __name__ == "__main__":
        main()
    """


def make_filter_pdb_script() -> str:
    return r"""
    #!/usr/bin/env python3

    from __future__ import annotations

    import argparse
    from pathlib import Path


    def parse_args() -> argparse.Namespace:
        parser = argparse.ArgumentParser(description="Keep protein ATOM/TER records from a PDB file.")
        parser.add_argument("--input", required=True, type=Path)
        parser.add_argument("--output", required=True, type=Path)
        return parser.parse_args()


    def main() -> None:
        args = parse_args()
        output_lines = []
        for line in args.input.read_text().splitlines():
            record = line[:6].strip()
            if record in {"ATOM", "TER"}:
                output_lines.append(line)
        output_lines.append("END")
        args.output.write_text("\n".join(output_lines) + "\n")


    if __name__ == "__main__":
        main()
    """


def make_topology_helper_script() -> str:
    return r"""
    #!/usr/bin/env python3

    from __future__ import annotations

    import argparse
    from pathlib import Path


    def parse_args() -> argparse.Namespace:
        parser = argparse.ArgumentParser(description="Set protein molecule count in a GROMACS topology.")
        parser.add_argument("--topology", required=True, type=Path)
        parser.add_argument("--copies", required=True, type=int)
        return parser.parse_args()


    def main() -> None:
        args = parse_args()
        lines = args.topology.read_text().splitlines()
        in_molecules = False
        replaced = False
        output = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("[") and stripped.endswith("]"):
                in_molecules = stripped.lower() == "[ molecules ]"
                output.append(line)
                continue
            if in_molecules and stripped and not stripped.startswith(";"):
                parts = stripped.split()
                if len(parts) >= 2 and parts[-1].isdigit():
                    parts[-1] = str(args.copies)
                    output.append(f"{parts[0]:<20}{parts[-1]}")
                    replaced = True
                    continue
            output.append(line)
        if not replaced:
            raise SystemExit(f"Could not find protein molecule entry in {args.topology}")
        args.topology.write_text("\n".join(output) + "\n")


    if __name__ == "__main__":
        main()
    """


def write_common_inputs(folder: Path, mdp_config: dict[str, Any]) -> None:
    write_file(folder / "em.mdp", make_mdp(require_mapping(mdp_config, "energy_minimization")))
    write_file(folder / "nvt.mdp", make_mdp(require_mapping(mdp_config, "equilibration")))
    write_file(folder / "md.mdp", make_mdp(require_mapping(mdp_config, "production")))
    write_file(folder / "filter_protein_pdb.py", make_filter_pdb_script(), executable=True)
    write_file(folder / "replicate_gro.py", make_replicate_gro_script(), executable=True)
    write_file(folder / "set_topology_copies.py", make_topology_helper_script(), executable=True)


def make_case_makefile(
    case: ProteinCase,
    build: GmxBuild,
    execution: dict[str, Any],
    campaign: dict[str, Any],
) -> str:
    execution = {**execution, **build.execution}
    gmxflags = execution.get("gmxflags", "")
    solvate_flags = execution.get("solvate_flags", "-cs spc216.gro")
    minimization_mdrun_flags = execution.get("minimization_mdrun_flags", "-nb cpu -v")
    production_mdrun_flags = execution.get("production_mdrun_flags", "-nb gpu -bonded gpu -update gpu -v")
    equilibration_mdrun_flags = execution.get("equilibration_mdrun_flags", production_mdrun_flags)
    energy_term = execution.get("energy_term", "Total-Energy")
    protein_pdb = Path(str(campaign["protein_pdb"])).resolve()
    force_field = campaign.get("force_field", "oplsaa")
    water_model = campaign.get("water_model", "spc")
    pdb2gmx_flags = campaign.get("pdb2gmx_flags", "-ignh")
    grompp_maxwarn = int(campaign.get("grompp_maxwarn", 1))

    return f"""
    GMX ?= {build.binary}
    GMXFLAGS ?= {gmxflags}
    SOLVATE_FLAGS ?= {solvate_flags}
    MINIMIZATION_MDRUN_FLAGS ?= {minimization_mdrun_flags}
    EQUILIBRATION_MDRUN_FLAGS ?= {equilibration_mdrun_flags}
    PRODUCTION_MDRUN_FLAGS ?= {production_mdrun_flags}
    ENERGY_TERM ?= {energy_term}

    PROTEIN_PDB = {protein_pdb}
    FORCE_FIELD = {force_field}
    WATER_MODEL = {water_model}
    PDB2GMX_FLAGS = {pdb2gmx_flags}
    GROMPP_MAXWARN = {grompp_maxwarn}

    COPIES = {case.copy_count}
    COPY_AXIS = {case.copy_axis}
    COPY_SPACING = {fmt(case.copy_spacing_nm)}
    X = {fmt(case.box_x_nm)}
    Y = {fmt(case.box_y_nm)}
    Z = {fmt(case.box_z_nm)}

    PROTEIN_INPUT = protein_input.pdb
    SINGLE_PROTEIN = protein_single.gro
    SOLUTE_CONF = protein_replicated.gro
    CONF = conf.gro
    TOP = topol.top

    EM_MDP = em.mdp
    EM_DEFFNM = em
    EM_TPR = $(EM_DEFFNM).tpr
    EM_CONF = $(EM_DEFFNM).gro

    NVT_MDP = nvt.mdp
    NVT_DEFFNM = nvt
    NVT_TPR = $(NVT_DEFFNM).tpr
    NVT_CONF = $(NVT_DEFFNM).gro

    MDP = md.mdp
    DEFFNM = md
    TPR = topol.tpr
    PBC_TRAJ = $(DEFFNM)_pbc.xtc
    ENERGY_XVG = total_energy.xvg
    ENERGY_TXT = total_energy.dat
    DRIFT_TXT = energy_drift.txt

    .PHONY: all info topology solvate minimize equilibrate prepare run pbc energy clean distclean

    all: energy

    info:
    \t@echo "Build = {build.name}"
    \t@echo "GMX = $(GMX)"
    \t@echo "Protein copies = $(COPIES)"
    \t@echo "Box = $(X) $(Y) $(Z) nm"
    \t@echo "PRODUCTION_MDRUN_FLAGS = $(PRODUCTION_MDRUN_FLAGS)"

    topology:
    \t./filter_protein_pdb.py --input $(PROTEIN_PDB) --output $(PROTEIN_INPUT)
    \t$(GMX) pdb2gmx $(GMXFLAGS) \\
    \t\t-f $(PROTEIN_INPUT) \\
    \t\t-o $(SINGLE_PROTEIN) \\
    \t\t-p $(TOP) \\
    \t\t-ff $(FORCE_FIELD) \\
    \t\t-water $(WATER_MODEL) \\
    \t\t$(PDB2GMX_FLAGS)
    \t./replicate_gro.py \\
    \t\t--input $(SINGLE_PROTEIN) \\
    \t\t--output $(SOLUTE_CONF) \\
    \t\t--copies $(COPIES) \\
    \t\t--axis $(COPY_AXIS) \\
    \t\t--spacing $(COPY_SPACING) \\
    \t\t--box $(X) $(Y) $(Z)
    \t./set_topology_copies.py --topology $(TOP) --copies $(COPIES)

    solvate: topology
    \t$(GMX) solvate $(GMXFLAGS) $(SOLVATE_FLAGS) \\
    \t\t-cp $(SOLUTE_CONF) \\
    \t\t-box $(X) $(Y) $(Z) \\
    \t\t-o $(CONF) \\
    \t\t-p $(TOP)

    minimize: solvate
    \t$(GMX) grompp $(GMXFLAGS) \\
    \t\t-f $(EM_MDP) \\
    \t\t-c $(CONF) \\
    \t\t-p $(TOP) \\
    \t\t-o $(EM_TPR) \\
    \t\t-maxwarn $(GROMPP_MAXWARN)
    \t$(GMX) mdrun $(GMXFLAGS) \\
    \t\t-s $(EM_TPR) \\
    \t\t-deffnm $(EM_DEFFNM) \\
    \t\t$(MINIMIZATION_MDRUN_FLAGS)

    equilibrate: minimize
    \t$(GMX) grompp $(GMXFLAGS) \\
    \t\t-f $(NVT_MDP) \\
    \t\t-c $(EM_CONF) \\
    \t\t-p $(TOP) \\
    \t\t-o $(NVT_TPR) \\
    \t\t-maxwarn $(GROMPP_MAXWARN)
    \t$(GMX) mdrun $(GMXFLAGS) \\
    \t\t-s $(NVT_TPR) \\
    \t\t-deffnm $(NVT_DEFFNM) \\
    \t\t$(EQUILIBRATION_MDRUN_FLAGS)

    prepare: equilibrate
    \t$(GMX) grompp $(GMXFLAGS) \\
    \t\t-f $(MDP) \\
    \t\t-c $(NVT_CONF) \\
    \t\t-p $(TOP) \\
    \t\t-o $(TPR) \\
    \t\t-maxwarn $(GROMPP_MAXWARN)

    run: prepare
    \t$(GMX) mdrun $(GMXFLAGS) \\
    \t\t-s $(TPR) \\
    \t\t-deffnm $(DEFFNM) \\
    \t\t$(PRODUCTION_MDRUN_FLAGS)

    pbc: run
    \tprintf "System\\nSystem\\n" | $(GMX) trjconv $(GMXFLAGS) \\
    \t\t-s $(TPR) \\
    \t\t-f $(DEFFNM).trr \\
    \t\t-o $(PBC_TRAJ) \\
    \t\t-pbc mol \\
    \t\t-ur compact \\
    \t\t-center

    energy: pbc
    \tprintf "$(ENERGY_TERM)\\n0\\n" | $(GMX) energy $(GMXFLAGS) \\
    \t\t-f $(DEFFNM).edr \\
    \t\t-o $(ENERGY_XVG) \\
    \t\t-xvg none
    \tawk 'BEGIN {{ first_set=0 }} /^[@#]/ {{ next }} NF >= 2 {{ if (!first_set) {{ first_t=$$1; first_e=$$2; first_set=1 }} last_t=$$1; last_e=$$2 }} END {{ if (!first_set) {{ print "No energy samples found" > "/dev/stderr"; exit 1 }} drift=last_e-first_e; duration=last_t-first_t; printf "first_time_ps %.8g\\nfirst_energy_kj_mol %.12g\\nlast_time_ps %.8g\\nlast_energy_kj_mol %.12g\\ndrift_kj_mol %.12g\\ndrift_per_ps_kj_mol %.12g\\n", first_t, first_e, last_t, last_e, drift, (duration == 0 ? 0 : drift / duration) }}' $(ENERGY_XVG) > $(DRIFT_TXT)
    \tcp $(ENERGY_XVG) $(ENERGY_TXT)

    clean:
    \trm -f $(PROTEIN_INPUT) $(SINGLE_PROTEIN) $(SOLUTE_CONF) $(CONF) $(TOP) posre.itp
    \trm -f $(EM_TPR) $(EM_CONF) em.log em.edr em.trr em.xtc em.cpt em_prev.cpt
    \trm -f $(NVT_TPR) $(NVT_CONF) nvt.log nvt.edr nvt.trr nvt.xtc nvt.cpt nvt_prev.cpt
    \trm -f $(TPR)
    \trm -f md.log md.edr md.trr md.xtc md.gro md.cpt md_prev.cpt $(PBC_TRAJ)
    \trm -f $(ENERGY_XVG) $(ENERGY_TXT) $(DRIFT_TXT) mdout.mdp

    distclean: clean
    \trm -f \\#*
    """


def write_case(
    folder: Path,
    case: ProteinCase,
    build: GmxBuild,
    mdp_config: dict[str, Any],
    execution: dict[str, Any],
    campaign: dict[str, Any],
) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    write_common_inputs(folder, mdp_config)
    write_file(folder / "Makefile", make_case_makefile(case, build, execution, campaign))
    write_file(
        folder / "run_local.sh",
        """
        #!/usr/bin/env bash
        set -euo pipefail
        cd "$(dirname "$0")"
        make all
        """,
        executable=True,
    )
    write_file(
        folder / "run_case.sh",
        """
        #!/usr/bin/env bash
        set -euo pipefail
        exec "$(dirname "$0")/run_local.sh"
        """,
        executable=True,
    )


def write_manifest(output: Path, case_dirs: list[Path]) -> None:
    lines = ["# Generated GROMACS protein-water test cases", ""]
    lines.extend(str(path) for path in case_dirs)
    write_file(output / "cases.txt", "\n".join(lines) + "\n")


def write_slurm_script(output: Path, case_dirs: list[Path], slurm_config: dict[str, Any]) -> None:
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={slurm_config.get('job_name', 'gmx-protein-water')}",
        f"#SBATCH --array=0-{len(case_dirs) - 1}",
        f"#SBATCH --nodes={slurm_config.get('nodes', 1)}",
        f"#SBATCH --ntasks={slurm_config.get('ntasks', 1)}",
        f"#SBATCH --cpus-per-task={slurm_config.get('cpus_per_task', 8)}",
        f"#SBATCH --gres={slurm_config.get('gres', 'gpu:1')}",
        f"#SBATCH --time={slurm_config.get('time', '02:00:00')}",
        f"#SBATCH --output={slurm_config.get('output', 'slurm-%A_%a.out')}",
        f"#SBATCH --error={slurm_config.get('error', 'slurm-%A_%a.err')}",
    ]
    partition = slurm_config.get("partition")
    if partition:
        lines.append(f"#SBATCH --partition={partition}")

    lines.extend(["", "set -euo pipefail", "", "CASE_DIRS=("])
    lines.extend(f'    "{path}"' for path in case_dirs)
    lines.extend(
        [
            ")",
            "",
            'CASE_DIR="${CASE_DIRS[$SLURM_ARRAY_TASK_ID]}"',
            "",
            'echo "Running case $SLURM_ARRAY_TASK_ID in $CASE_DIR"',
            'make -C "$CASE_DIR" all',
            "",
        ]
    )
    write_file(output / "submit_slurm.sh", "\n".join(lines), executable=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate a GROMACS protein-water test campaign from split JSON configuration files."
    )
    parser.add_argument(
        "--general-config",
        default="config/general_config.json",
        help="Path to the general JSON configuration. Default: config/general_config.json",
    )
    parser.add_argument(
        "--system-config",
        default="config/protein_water_system_config.json",
        help="Path to the protein-water system JSON configuration. Default: config/protein_water_system_config.json",
    )
    parser.add_argument(
        "--copy-count",
        type=int,
        default=None,
        help="Override the configured copy_counts list with a single N value.",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Run all generated cases locally after creating the campaign.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    general_config_path = Path(args.general_config)
    system_config_path = Path(args.system_config)
    general_config = load_config(general_config_path)
    system_config = load_config(system_config_path)
    campaign = require_mapping(system_config, "campaign")
    mdp_config = require_mapping(system_config, "mdp")
    execution = require_mapping(general_config, "execution")
    slurm_config = require_mapping(general_config, "slurm")

    protein_pdb = Path(str(campaign["protein_pdb"]))
    if not protein_pdb.exists():
        raise FileNotFoundError(f"Protein PDB not found: {protein_pdb}")

    output = Path(campaign.get("output", "protein_water_fixed_precision_campaign")).resolve()
    cases = make_cases(campaign, args.copy_count)
    builds = parse_builds(general_config)

    case_dirs: list[Path] = []
    for case in cases:
        for build in builds:
            case_dir = output / build.name / case.folder_name
            write_case(case_dir, case, build, mdp_config, execution, campaign)
            case_dirs.append(case_dir)

    write_manifest(output, case_dirs)
    write_slurm_script(output, case_dirs, slurm_config)

    print(f"Loaded general configuration: {general_config_path}")
    print(f"Loaded system configuration: {system_config_path}")
    print(f"Created campaign: {output}")
    print(f"Generated {len(cases)} copy counts x {len(builds)} builds = {len(case_dirs)} case folders")
    print()
    print("Run locally:")
    print(f"  {case_dirs[0] / 'run_local.sh'}")
    print(f"  make -C {case_dirs[0]} all")
    print()
    print("Run on SLURM:")
    print(f"  sbatch {output / 'submit_slurm.sh'}")

    if args.run:
        for case_dir in case_dirs:
            subprocess.run(["make", "-C", str(case_dir), "all"], check=True)


if __name__ == "__main__":
    main()
