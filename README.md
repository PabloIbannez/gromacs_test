# GROMACS Water Fixed-Precision Test Campaigns

This repository generates pure-water GROMACS test campaigns for comparing a
reference binary against a candidate binary. The default campaign sweeps the
`x` box length with `5x5x5`, `10x5x5`, `50x5x5`, and `100x5x5` nm water boxes,
creating one folder per build and box size.

All campaign, simulation, execution, and SLURM parameters live in
`campaign_config.json`.

## Generate a Campaign

Edit the GROMACS binary paths in `campaign_config.json`:

```json
"builds": [
  {
    "name": "reference",
    "binary": "/path/to/reference/gmx"
  },
  {
    "name": "candidate",
    "binary": "/path/to/candidate/gmx"
  }
]
```

Then generate the campaign:

```bash
python3 generate_test.py --config campaign_config.json
```

The default output directory is `water_fixed_precision_campaign`. It contains:

- `reference/water_*`: runs using the reference binary.
- `candidate/water_*`: runs using the candidate binary.
- `cases.txt`: the generated case list.
- `submit_slurm.sh`: a SLURM array job that launches every case.

Each case writes:

- `total_energy.xvg`: total energy extracted from `md.edr`.
- `total_energy.dat`: a copy of the same numeric energy series.
- `energy_drift.txt`: first and last energy values, total drift, and drift per
  ps.

## Simulation Parameters

The box sweep is configured under `campaign`:

```json
"campaign": {
  "output": "water_fixed_precision_campaign",
  "base_box_nm": [5.0, 5.0, 5.0],
  "sweep_axis": "x",
  "axis_lengths_nm": [5.0, 10.0, 50.0, 100.0]
}
```

The generated `em.mdp` and `md.mdp` files are written from the `mdp` section.
The default production run is NVE and writes total energy every 100 steps.

Electrostatics use the Verlet neighbor list with a cutoff instead of PME:

```json
"coulombtype": "Cut-off",
"coulomb-modifier": "Potential-shift",
"rlist": 1.0,
"rcoulomb": 1.0
```

## GPU Path

The production `mdrun` flags are configured under `execution`:

```json
"production_mdrun_flags": "-nb gpu -bonded gpu -update gpu -v"
```

This keeps nonbonded, bonded, and update work on the GPU while avoiding PME.

## SLURM

SLURM parameters are configured under `slurm`:

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
