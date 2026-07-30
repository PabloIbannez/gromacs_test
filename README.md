# GROMACS Fixed-Precision Test Campaigns

This repository generates GROMACS test campaigns for comparing a reference
binary against a candidate binary and a double-precision reference binary. It
currently supports pure-water systems and protein-in-water systems with repeated
copies of ubiquitin.

General execution and build parameters live in `config/general_config.json`.
System-specific parameters live in `config/water_system_config.json` and
`config/protein_water_system_config.json`.

## Shared Build Configuration

Edit the GROMACS binary paths in `config/general_config.json`:

```json
"builds": [
  {
    "name": "reference",
    "binary": "/path/to/reference/gmx",
    "execution": {
      "equilibration_mdrun_flags": "-nb gpu -bonded cpu -update gpu -v",
      "production_mdrun_flags": "-nb gpu -bonded cpu -update gpu -v"
    }
  },
  {
    "name": "candidate",
    "binary": "/path/to/candidate/gmx",
    "execution": {
      "equilibration_mdrun_flags": "-nb gpu -bonded cpu -update gpu -v",
      "production_mdrun_flags": "-nb gpu -bonded cpu -update gpu -v"
    }
  },
  {
    "name": "reference_double",
    "binary": "/path/to/double/reference/gmx_d",
    "execution": {
      "equilibration_mdrun_flags": "-nb cpu -bonded cpu -update cpu -v",
      "production_mdrun_flags": "-nb cpu -bonded cpu -update cpu -v"
    }
  }
]
```

## Generate a Water Campaign

The default water campaign sweeps the `x` box length with `5x3x3`, `10x3x3`,
`50x3x3`, `100x3x3`, and `150x3x3` nm water boxes, creating one folder per
build and box size.

Generate it with:

```bash
python3 generate/generate_water_test.py \
  --general-config config/general_config.json \
  --system-config config/water_system_config.json
```

The default output directory is `water_fixed_precision_campaign`. It contains:

- `reference/water_*`: runs using the reference binary.
- `candidate/water_*`: runs using the candidate binary.
- `reference_double/water_*`: runs using the double-precision reference binary.
- `cases.txt`: the generated case list.
- `submit_slurm.sh`: a SLURM array job that launches every case.
- `run_local.sh` in each case folder: launches that case from the folder.
- `run_case.sh` in each case folder: compatibility wrapper for `run_local.sh`.

Each case writes:

- `total_energy.xvg`: total energy extracted from `md.edr`.
- `total_energy.dat`: a copy of the same numeric energy series.
- `energy_drift.txt`: first and last energy values, total drift, and drift per
  ps.

Run one generated case from its folder:

```bash
cd water_fixed_precision_campaign/reference/water_100_3_3
./run_local.sh
```

Or from the repository root:

```bash
make -C water_fixed_precision_campaign/reference/water_100_3_3 all
```

## Simulation Parameters

The box sweep is configured under `campaign` in `config/water_system_config.json`:

```json
"campaign": {
  "output": "water_fixed_precision_campaign",
  "base_box_nm": [3.0, 3.0, 3.0],
  "sweep_axis": "x",
  "axis_lengths_nm": [5.0, 10.0, 50.0, 100.0, 150.0]
}
```

The generated `em.mdp`, `nvt.mdp`, and `md.mdp` files are written from the
`mdp` section. Each case minimizes, equilibrates with NVT, then runs production
and writes total energy every 100 steps.

Electrostatics use PME with the Verlet neighbor list:

```json
"coulombtype": "PME",
"fourierspacing": 0.12,
"pme-order": 4,
"rlist": 1.0,
"rcoulomb": 1.0
```

## GPU Path

The default `mdrun` flags are configured under `execution` in
`config/general_config.json` and can be overridden per build:

```json
"equilibration_mdrun_flags": "-nb gpu -bonded cpu -update gpu -v",
"production_mdrun_flags": "-nb gpu -bonded cpu -update gpu -v"
```

This keeps nonbonded and update work on the GPU for the single-precision builds,
while the double-precision reference defaults to CPU execution.

## Energy Plots

After cases complete, plot available production total-energy traces:

```bash
python3 analysis/analyze_energy.py water_fixed_precision_campaign
```

Individual plots show absolute total energy per atom. Combined plots subtract
the fitted line intercept from each trace, so traces with different absolute
energy offsets start from the same fitted origin.

Use `--extract-missing` to extract `total_energy.xvg` from completed cases that
only have `md.edr`.

Use `--extract-partial` to also include unfinished cases that already have
`md.edr`:

```bash
python3 analysis/analyze_energy.py water_fixed_precision_campaign --extract-partial
```

`--extract-missing` only operates on runs whose `md.log` says `Finished mdrun`.
`--extract-partial` also extracts from running or interrupted production runs
when `md.edr` exists, so the plot covers the available production time instead
of waiting for the final `md.gro`.

## SLURM

SLURM parameters are configured under `slurm` in `config/general_config.json`:

```json
"slurm": {
  "job_name": "gmx-water-fixed",
  "nodes": 1,
  "ntasks": 1,
  "cpus_per_task": 8,
  "gres": "gpu:1",
  "time": "02:00:00",
  "partition": "",
  "output": "slurm-%A_%a.out",
  "error": "slurm-%A_%a.err"
}
```

Submit all generated cases with:

```bash
sbatch water_fixed_precision_campaign/submit_slurm.sh
```

## Generate a Protein-Water Campaign

The protein-water campaign uses ubiquitin (`1UBQ`) as a compact structured
protein. The input PDB lives at `data/1ubq.pdb`; regenerate it from RCSB with:

```bash
python3 data/fetch_1ubq.py
```

Configure the number of horizontal copies in
`config/protein_water_system_config.json`:

```json
"copy_counts": [1, 2, 4],
"copy_axis": "x",
"copy_spacing_nm": 4.0
```

Generate all configured copy counts:

```bash
python3 generate/generate_protein_water_test.py \
  --general-config config/general_config.json \
  --system-config config/protein_water_system_config.json
```

Or generate one `N` explicitly:

```bash
python3 generate/generate_protein_water_test.py --copy-count 4
```

The default output directory is `protein_water_fixed_precision_campaign`. It
contains case folders such as:

- `reference/protein_1ubq_N4`
- `candidate/protein_1ubq_N4`
- `reference_double/protein_1ubq_N4`

Each protein case runs `pdb2gmx`, replicates the processed protein `N` times in
the configured horizontal axis, solvates, minimizes, equilibrates with NVT, runs
production, writes `md.trr`, postprocesses `md_pbc.xtc`, and extracts total
energy.

Run one generated protein case from the repository root:

```bash
make -C protein_water_fixed_precision_campaign/reference/protein_1ubq_N4 all
```

View a completed protein trajectory with VMD:

```bash
vmd protein_water_fixed_precision_campaign/reference/protein_1ubq_N4/md.gro \
    protein_water_fixed_precision_campaign/reference/protein_1ubq_N4/md_pbc.xtc
```

For a partial run that has `md.trr` but not `md.gro`, use:

```bash
vmd protein_water_fixed_precision_campaign/reference/protein_1ubq_N4/topol.tpr \
    protein_water_fixed_precision_campaign/reference/protein_1ubq_N4/md.trr
```

Plot partial protein-water energy traces with:

```bash
python3 analysis/analyze_energy.py protein_water_fixed_precision_campaign --extract-partial
```
