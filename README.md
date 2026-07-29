# GROMACS Water Fixed-Precision Test Campaigns

This repository generates pure-water GROMACS test campaigns for comparing a
reference binary against a candidate binary and a double-precision reference
binary. The default campaign sweeps the `x` box length with `5x3x3`,
`10x3x3`, `50x3x3`, `100x3x3`, and `150x3x3` nm water boxes, creating one
folder per build and box size.

General execution and build parameters live in `config/general_config.json`.
Water system and MDP parameters live in `config/water_system_config.json`.

## Generate a Campaign

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

Then generate the campaign:

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

Use `--extract-missing` to extract `total_energy.xvg` from completed cases that
only have `md.edr`.

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
