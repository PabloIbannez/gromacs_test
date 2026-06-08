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


@dataclass(frozen=True)
class WaterCase:
    x: float
    y: float
    z: float

    @property
    def folder_name(self) -> str:
        return f"water_{fmt(self.x)}_{fmt(self.y)}_{fmt(self.z)}"


def fmt(value: float) -> str:
    """Format a number for command-line values and folder names."""
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


def make_cases(base_box: list[Any], axis: str, axis_lengths: list[Any]) -> list[WaterCase]:
    if len(base_box) != 3:
        raise ValueError("'campaign.base_box_nm' must contain exactly three values")
    if axis not in {"x", "y", "z"}:
        raise ValueError("'campaign.sweep_axis' must be one of: x, y, z")

    axis_index = {"x": 0, "y": 1, "z": 2}[axis]
    cases = []
    for length in axis_lengths:
        box = [float(value) for value in base_box]
        box[axis_index] = float(length)
        cases.append(WaterCase(*box))
    return cases


def make_mdp(parameters: dict[str, Any]) -> str:
    lines = []
    for key, value in parameters.items():
        lines.append(f"{key:<24}= {value}")
    return "\n".join(lines) + "\n"


def write_common_inputs(folder: Path, mdp_config: dict[str, Any]) -> None:
    topology = mdp_config.get("topology", {})
    if not isinstance(topology, dict):
        raise ValueError("'mdp.topology' must be an object")

    force_field = topology.get("force_field", "oplsaa.ff/forcefield.itp")
    water_model = topology.get("water_model", "oplsaa.ff/spc.itp")
    system_name = topology.get("system_name", "Pure water fixed-precision coordinate test")

    write_file(
        folder / "topol.template.top",
        f"""
        #include "{force_field}"
        #include "{water_model}"

        [ system ]
        {system_name}

        [ molecules ]
        """,
    )

    write_file(folder / "em.mdp", make_mdp(require_mapping(mdp_config, "energy_minimization")))
    write_file(folder / "md.mdp", make_mdp(require_mapping(mdp_config, "production")))


def make_case_makefile(
    case: WaterCase,
    build: GmxBuild,
    execution: dict[str, Any],
) -> str:
    gmxflags = execution.get("gmxflags", "")
    solvate_flags = execution.get("solvate_flags", "-cs spc216.gro")
    minimization_mdrun_flags = execution.get("minimization_mdrun_flags", "-nb cpu -v")
    production_mdrun_flags = execution.get("production_mdrun_flags", "-nb gpu -bonded gpu -update gpu -v")
    energy_term = execution.get("energy_term", "Total-Energy")

    return f"""
    GMX ?= {build.binary}
    GMXFLAGS ?= {gmxflags}
    SOLVATE_FLAGS ?= {solvate_flags}
    MINIMIZATION_MDRUN_FLAGS ?= {minimization_mdrun_flags}
    PRODUCTION_MDRUN_FLAGS ?= {production_mdrun_flags}
    ENERGY_TERM ?= {energy_term}

    X = {fmt(case.x)}
    Y = {fmt(case.y)}
    Z = {fmt(case.z)}

    TOP = topol.top
    TEMPLATE_TOP = topol.template.top
    CONF = conf.gro

    EM_MDP = em.mdp
    EM_DEFFNM = em
    EM_TPR = $(EM_DEFFNM).tpr
    EM_CONF = $(EM_DEFFNM).gro

    MDP = md.mdp
    DEFFNM = md
    TPR = topol.tpr
    ENERGY_XVG = total_energy.xvg
    ENERGY_TXT = total_energy.dat
    DRIFT_TXT = energy_drift.txt

    .PHONY: all info solvate minimize prepare run energy clean distclean

    all: energy

    info:
    \t@echo "Build = {build.name}"
    \t@echo "GMX = $(GMX)"
    \t@echo "Box = $(X) $(Y) $(Z) nm"
    \t@echo "PRODUCTION_MDRUN_FLAGS = $(PRODUCTION_MDRUN_FLAGS)"

    solvate:
    \tcp $(TEMPLATE_TOP) $(TOP)
    \t$(GMX) solvate $(GMXFLAGS) $(SOLVATE_FLAGS) \\
    \t\t-box $(X) $(Y) $(Z) \\
    \t\t-o $(CONF) \\
    \t\t-p $(TOP)

    minimize: solvate
    \t$(GMX) grompp $(GMXFLAGS) \\
    \t\t-f $(EM_MDP) \\
    \t\t-c $(CONF) \\
    \t\t-p $(TOP) \\
    \t\t-o $(EM_TPR)
    \t$(GMX) mdrun $(GMXFLAGS) \\
    \t\t-s $(EM_TPR) \\
    \t\t-deffnm $(EM_DEFFNM) \\
    \t\t$(MINIMIZATION_MDRUN_FLAGS)

    prepare: minimize
    \t$(GMX) grompp $(GMXFLAGS) \\
    \t\t-f $(MDP) \\
    \t\t-c $(EM_CONF) \\
    \t\t-p $(TOP) \\
    \t\t-o $(TPR)

    run: prepare
    \t$(GMX) mdrun $(GMXFLAGS) \\
    \t\t-s $(TPR) \\
    \t\t-deffnm $(DEFFNM) \\
    \t\t$(PRODUCTION_MDRUN_FLAGS)

    energy: run
    \tprintf "$(ENERGY_TERM)\\n0\\n" | $(GMX) energy $(GMXFLAGS) \\
    \t\t-f $(DEFFNM).edr \\
    \t\t-o $(ENERGY_XVG) \\
    \t\t-xvg none
    \tawk 'BEGIN {{ first_set=0 }} /^[@#]/ {{ next }} NF >= 2 {{ if (!first_set) {{ first_t=$$1; first_e=$$2; first_set=1 }} last_t=$$1; last_e=$$2 }} END {{ if (!first_set) {{ print "No energy samples found" > "/dev/stderr"; exit 1 }} drift=last_e-first_e; duration=last_t-first_t; printf "first_time_ps %.8g\\nfirst_energy_kj_mol %.12g\\nlast_time_ps %.8g\\nlast_energy_kj_mol %.12g\\ndrift_kj_mol %.12g\\ndrift_per_ps_kj_mol %.12g\\n", first_t, first_e, last_t, last_e, drift, (duration == 0 ? 0 : drift / duration) }}' $(ENERGY_XVG) > $(DRIFT_TXT)
    \tcp $(ENERGY_XVG) $(ENERGY_TXT)

    clean:
    \trm -f $(CONF) $(TOP)
    \trm -f $(EM_TPR) $(EM_CONF) em.log em.edr em.trr em.xtc em.cpt em_prev.cpt
    \trm -f $(TPR)
    \trm -f md.log md.edr md.trr md.xtc md.gro md.cpt md_prev.cpt
    \trm -f $(ENERGY_XVG) $(ENERGY_TXT) $(DRIFT_TXT)

    distclean: clean
    \trm -f \\#*
    """


def write_case(
    folder: Path,
    case: WaterCase,
    build: GmxBuild,
    mdp_config: dict[str, Any],
    execution: dict[str, Any],
) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    write_common_inputs(folder, mdp_config)
    write_file(folder / "Makefile", make_case_makefile(case, build, execution))
    write_file(
        folder / "run_case.sh",
        """
        #!/usr/bin/env bash
        set -euo pipefail
        make -C "$(dirname "$0")" all
        """,
        executable=True,
    )


def write_manifest(output: Path, case_dirs: list[Path]) -> None:
    lines = ["# Generated GROMACS water test cases", ""]
    lines.extend(str(path) for path in case_dirs)
    write_file(output / "cases.txt", "\n".join(lines) + "\n")


def write_slurm_script(output: Path, case_dirs: list[Path], slurm_config: dict[str, Any]) -> None:
    array_max = len(case_dirs) - 1
    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={slurm_config.get('job_name', 'gmx-water-fixed')}",
        f"#SBATCH --array=0-{array_max}",
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

    lines.extend(
        [
            "",
            "set -euo pipefail",
            "",
            "CASE_DIRS=(",
        ]
    )
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


def parse_builds(config: dict[str, Any]) -> list[GmxBuild]:
    raw_builds = require_list(config, "builds")
    builds = []
    for entry in raw_builds:
        if not isinstance(entry, dict):
            raise ValueError("Each build entry must be an object")
        builds.append(GmxBuild(str(entry["name"]), str(entry["binary"])))
    if not builds:
        raise ValueError("At least one build must be configured")
    return builds


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate a GROMACS pure-water GPU test campaign from a JSON configuration."
    )
    parser.add_argument(
        "--config",
        default="campaign_config.json",
        help="Path to the campaign JSON configuration. Default: campaign_config.json",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Run all generated cases locally after creating the campaign.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    config_path = Path(args.config)
    config = load_config(config_path)
    campaign = require_mapping(config, "campaign")
    mdp_config = require_mapping(config, "mdp")
    execution = require_mapping(config, "execution")
    slurm_config = require_mapping(config, "slurm")

    output = Path(campaign.get("output", "water_fixed_precision_campaign")).resolve()
    cases = make_cases(
        require_list(campaign, "base_box_nm"),
        str(campaign.get("sweep_axis", "x")),
        require_list(campaign, "axis_lengths_nm"),
    )
    builds = parse_builds(config)

    case_dirs: list[Path] = []
    for case in cases:
        for build in builds:
            case_dir = output / build.name / case.folder_name
            write_case(case_dir, case, build, mdp_config, execution)
            case_dirs.append(case_dir)

    write_manifest(output, case_dirs)
    write_slurm_script(output, case_dirs, slurm_config)

    print(f"Loaded configuration: {config_path}")
    print(f"Created campaign: {output}")
    print(f"Generated {len(cases)} box sizes x {len(builds)} builds = {len(case_dirs)} case folders")
    print()
    print("Run locally:")
    print(f"  make -C {case_dirs[0]} all")
    print()
    print("Run on SLURM:")
    print(f"  sbatch {output / 'submit_slurm.sh'}")
    print()
    print("Energy outputs per case:")
    print("  total_energy.xvg")
    print("  total_energy.dat")
    print("  energy_drift.txt")

    if args.run:
        for case_dir in case_dirs:
            subprocess.run(["make", "-C", str(case_dir), "all"], check=True)


if __name__ == "__main__":
    main()
