#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-gromacs-energy")

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    print(
        "Missing required Python package: matplotlib. Install it before running this script.",
        file=sys.stderr,
    )
    raise SystemExit(1)


@dataclass(frozen=True)
class EnergyTrace:
    case_dir: Path
    relative_name: str
    build_name: str
    box_name: str
    natoms: int
    times: list[float]
    energies: list[float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot production-phase total-energy traces for completed campaign cases."
    )
    parser.add_argument(
        "campaign",
        nargs="?",
        default="water_fixed_precision_campaign",
        type=Path,
        help="Campaign directory to scan.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for plots and summary files. Defaults to CAMPAIGN/energy_plots.",
    )
    parser.add_argument(
        "--energy-term",
        default="Total-Energy",
        help="GROMACS energy term to extract when total_energy.xvg is missing.",
    )
    parser.add_argument(
        "--gmx",
        default=None,
        help="GROMACS binary to use for extraction. Defaults to the case Makefile GMX value.",
    )
    parser.add_argument(
        "--extract-missing",
        action="store_true",
        help="Extract total_energy.xvg from md.edr for completed cases where it is missing.",
    )
    return parser.parse_args()


def discover_case_dirs(campaign_dir: Path) -> list[Path]:
    manifest = campaign_dir / "cases.txt"
    if manifest.exists():
        case_dirs = []
        for line in manifest.read_text().splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            path = Path(stripped)
            if not path.is_absolute():
                path = campaign_dir.parent / path
            if (path / "Makefile").exists():
                case_dirs.append(path)
        return sorted(case_dirs)

    return sorted(path.parent for path in campaign_dir.rglob("Makefile"))


def production_finished(case_dir: Path) -> bool:
    if (case_dir / "total_energy.xvg").exists():
        return True

    md_log = case_dir / "md.log"
    if not md_log.exists() or not (case_dir / "md.edr").exists():
        return False

    try:
        log_text = md_log.read_text(errors="replace")
    except OSError:
        return False
    return "Finished mdrun" in log_text and "Segmentation fault" not in log_text


def make_relative_name(campaign_dir: Path, case_dir: Path) -> str:
    try:
        return str(case_dir.relative_to(campaign_dir))
    except ValueError:
        return str(case_dir)


def case_parts(campaign_dir: Path, case_dir: Path) -> tuple[str, str]:
    relative = Path(make_relative_name(campaign_dir, case_dir))
    if len(relative.parts) >= 2:
        return relative.parts[0], relative.parts[1]
    if len(relative.parts) == 1:
        return "unknown", relative.parts[0]
    return "unknown", case_dir.name


def read_gro_atom_count(path: Path) -> int:
    with path.open() as handle:
        next(handle, None)
        atom_count_line = next(handle, None)
    if atom_count_line is None:
        raise ValueError(f"Could not read atom count from {path}")
    try:
        return int(atom_count_line.strip())
    except ValueError as error:
        raise ValueError(f"Invalid atom count in {path}: {atom_count_line!r}") from error


def read_atom_count(case_dir: Path) -> int:
    for filename in ("md.gro", "nvt.gro", "em.gro", "conf.gro"):
        path = case_dir / filename
        if path.exists():
            return read_gro_atom_count(path)
    raise ValueError(f"No .gro file found to determine atom count in {case_dir}")


def read_makefile_gmx(case_dir: Path) -> str | None:
    makefile = case_dir / "Makefile"
    if not makefile.exists():
        return None

    for line in makefile.read_text().splitlines():
        match = re.match(r"\s*GMX\s*\?=\s*(.+?)\s*$", line)
        if match:
            return str(Path(match.group(1)).expanduser())
    return None


def extract_energy_xvg(case_dir: Path, gmx: str, energy_term: str) -> Path:
    output = case_dir / "total_energy.xvg"
    command = [gmx, "energy", "-f", "md.edr", "-o", output.name, "-xvg", "none"]
    subprocess.run(
        command,
        cwd=case_dir,
        input=f"{energy_term}\n0\n",
        text=True,
        check=True,
    )
    return output


def read_xvg(path: Path) -> tuple[list[float], list[float]]:
    times = []
    energies = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith(("@", "#")):
            continue
        columns = stripped.split()
        if len(columns) < 2:
            continue
        try:
            times.append(float(columns[0]))
            energies.append(float(columns[1]))
        except ValueError as error:
            raise ValueError(f"Could not parse {path}:{line_number}: {line}") from error

    if not times:
        raise ValueError(f"No energy samples found in {path}")
    return times, energies


def load_trace(
    campaign_dir: Path,
    case_dir: Path,
    energy_term: str,
    explicit_gmx: str | None,
    extract_missing: bool,
) -> EnergyTrace | None:
    energy_xvg = case_dir / "total_energy.xvg"
    if not energy_xvg.exists():
        if not extract_missing:
            return None
        gmx = explicit_gmx or read_makefile_gmx(case_dir)
        if gmx is None:
            raise ValueError(f"No GMX binary found for {case_dir}")
        energy_xvg = extract_energy_xvg(case_dir, gmx, energy_term)

    times, energies = read_xvg(energy_xvg)
    build_name, box_name = case_parts(campaign_dir, case_dir)
    return EnergyTrace(
        case_dir=case_dir,
        relative_name=make_relative_name(campaign_dir, case_dir),
        build_name=build_name,
        box_name=box_name,
        natoms=read_atom_count(case_dir),
        times=times,
        energies=energies,
    )


def safe_plot_name(relative_name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", relative_name)


def build_sort_key(trace: EnergyTrace) -> tuple[int, str]:
    order = {"reference": 0, "candidate": 1, "ref_dbl": 2}
    return order.get(trace.build_name, 100), trace.build_name


def energy_per_atom(trace: EnergyTrace) -> list[float]:
    return [energy / trace.natoms for energy in trace.energies]


def drift_per_atom(trace: EnergyTrace) -> list[float]:
    values = energy_per_atom(trace)
    initial = values[0]
    return [value - initial for value in values]


def plot_trace(trace: EnergyTrace, output_dir: Path) -> Path:
    output = output_dir / f"{safe_plot_name(trace.relative_name)}_total_energy.png"
    fig, axis = plt.subplots(figsize=(9, 5))
    axis.plot(trace.times, energy_per_atom(trace), linewidth=1.5)
    axis.set_title(trace.relative_name)
    axis.set_xlabel("Production time (ps)")
    axis.set_ylabel("Total energy per atom (kJ/mol/atom)")
    axis.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return output


def plot_combined(traces: list[EnergyTrace], output_dir: Path) -> Path | None:
    if not traces:
        return None

    output = output_dir / "all_completed_total_energy.png"
    fig, axis = plt.subplots(figsize=(10, 6))
    for trace in traces:
        axis.plot(trace.times, energy_per_atom(trace), linewidth=1.2, label=trace.relative_name)
    axis.set_xlabel("Production time (ps)")
    axis.set_ylabel("Total energy per atom (kJ/mol/atom)")
    axis.grid(True, alpha=0.3)
    axis.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return output


def plot_box_group(box_name: str, traces: list[EnergyTrace], output_dir: Path) -> Path:
    output = output_dir / f"{safe_plot_name(box_name)}_total_energy_by_build.png"
    ordered = sorted(traces, key=build_sort_key)
    fig, axis = plt.subplots(figsize=(10, 6))
    for trace in ordered:
        total_drift_per_atom = (trace.energies[-1] - trace.energies[0]) / trace.natoms
        label = f"{trace.build_name} ({total_drift_per_atom:.4g} kJ/mol/atom drift)"
        axis.plot(trace.times, drift_per_atom(trace), linewidth=1.4, label=label)
    axis.set_title(box_name)
    axis.set_xlabel("Production time (ps)")
    axis.set_ylabel("Energy drift per atom from t=0 (kJ/mol/atom)")
    axis.grid(True, alpha=0.3)
    axis.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return output


def plot_box_groups(traces: list[EnergyTrace], output_dir: Path) -> dict[str, Path]:
    traces_by_box: dict[str, list[EnergyTrace]] = {}
    for trace in traces:
        traces_by_box.setdefault(trace.box_name, []).append(trace)

    plot_paths = {}
    for box_name, box_traces in sorted(traces_by_box.items()):
        plot_paths[box_name] = plot_box_group(box_name, box_traces, output_dir)
    return plot_paths


def write_summary(traces: list[EnergyTrace], output_dir: Path, plot_paths: dict[str, Path]) -> Path:
    output = output_dir / "energy_summary.csv"
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "case",
                "build",
                "box",
                "natoms",
                "samples",
                "first_time_ps",
                "first_energy_kj_mol",
                "last_time_ps",
                "last_energy_kj_mol",
                "drift_kj_mol",
                "drift_per_atom_kj_mol",
                "drift_per_ps_kj_mol",
                "drift_per_atom_per_ps_kj_mol",
                "plot",
                "box_plot",
            ]
        )
        for trace in traces:
            duration = trace.times[-1] - trace.times[0]
            drift = trace.energies[-1] - trace.energies[0]
            drift_per_ps = 0.0 if duration == 0 else drift / duration
            drift_per_atom = drift / trace.natoms
            drift_per_atom_per_ps = 0.0 if duration == 0 else drift_per_atom / duration
            writer.writerow(
                [
                    trace.relative_name,
                    trace.build_name,
                    trace.box_name,
                    trace.natoms,
                    len(trace.times),
                    trace.times[0],
                    trace.energies[0],
                    trace.times[-1],
                    trace.energies[-1],
                    drift,
                    drift_per_atom,
                    drift_per_ps,
                    drift_per_atom_per_ps,
                    plot_paths[trace.relative_name].name,
                    plot_paths[trace.box_name].name,
                ]
            )
    return output


def main() -> None:
    args = parse_args()
    campaign_dir = args.campaign.resolve()
    output_dir = args.output_dir or campaign_dir / "energy_plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    case_dirs = discover_case_dirs(campaign_dir)
    completed = [case_dir for case_dir in case_dirs if production_finished(case_dir)]
    traces = []
    skipped = []
    plot_paths = {}

    for case_dir in completed:
        trace = load_trace(
            campaign_dir=campaign_dir,
            case_dir=case_dir,
            energy_term=args.energy_term,
            explicit_gmx=args.gmx,
            extract_missing=args.extract_missing,
        )
        if trace is None:
            skipped.append(make_relative_name(campaign_dir, case_dir))
            continue
        traces.append(trace)
        plot_paths[trace.relative_name] = plot_trace(trace, output_dir)

    if not traces:
        print(f"Scanned cases: {len(case_dirs)}")
        print(f"Completed production cases: {len(completed)}")
        print("Plotted cases: 0")
        print(f"Output directory: {output_dir}")
        if skipped:
            print("Completed cases skipped because total_energy.xvg is missing:")
            for name in skipped:
                print(f"  {name}")
            print("Rerun with --extract-missing to extract total_energy.xvg from md.edr.")
        return

    plot_combined(traces, output_dir)
    plot_paths.update(plot_box_groups(traces, output_dir))
    summary = write_summary(traces, output_dir, plot_paths)

    print(f"Scanned cases: {len(case_dirs)}")
    print(f"Completed production cases: {len(completed)}")
    print(f"Plotted cases: {len(traces)}")
    print(f"Output directory: {output_dir}")
    print(f"Summary: {summary}")
    if skipped:
        print("Completed cases skipped because total_energy.xvg is missing:")
        for name in skipped:
            print(f"  {name}")
        print("Rerun with --extract-missing to extract total_energy.xvg from md.edr.")


if __name__ == "__main__":
    main()
