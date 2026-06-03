# Spine_biomech

MATLAB workflow for replicating spine kinematics methods described in the references under `ref/`, using C3D marker trajectories and spline-based local coordinate systems.

## Goal

This repository is intended to reproduce key processing steps reported in the reference publications/thesis files in `ref/`:

- `ref/s10439-017-1970-x.pdf`
- `ref/Maksimovic_Mina_2021_thesis.pdf`

Current implementation focus:

- Build spine marker sets from C3D labels (`L`, midline spinous, `R`)
- Filter marker trajectories (4th-order low-pass Butterworth)
- Fit left/mid/right spline curves along the spine
- Construct vertebral local coordinate systems (LCS)
- Compute thoracic and lumbar 3D angles over time
- Visualize static frames and dynamic animation

## Repository Structure

```text
LICENSE
README.md
ref/                         # Reference publications used for replication targets
sample_data/                 # Example C3D files
src/
	Spine_biomech_test.m       # Static/single-frame pipeline and sanity checks
	Spine_biomech_dyn.m        # Dynamic trial processing and animation
	funs/                      # Helper functions (legacy + plotting/utilities)
```

## Requirements

- MATLAB (recommended R2021b or newer)
- Signal Processing Toolbox (`butter`, `filtfilt`)
- Curve Fitting Toolbox / Spline functions (`spap2`, `fnval`, `fnder`)
- `ezc3dRead` available on MATLAB path (from EZC3D MATLAB bindings)

## Data Expectations

The scripts assume marker labels follow this convention:

- Left markers end with `L` (example: `T12L`)
- Right markers end with `R` (example: `T12R`)
- Midline spinous markers have no side suffix (example: `T12`)

Expected vertebral ordering is:

`C6, C7, T1..T12, L1..L5, S1..S5`

Only levels present in the C3D file are used.

## Quick Start

1. Open MATLAB in the repository root.
2. Add source folders and EZC3D to path.
3. Run the static script first to validate labels and coordinate systems.
4. Run the dynamic script to compute time-series angles.

Example session:

```matlab
cd('.../Spine_biomech');
addpath(genpath('src'));
% addpath('.../path/to/ezc3d/matlab');

% Static quality check
run('src/Spine_biomech_test.m');

% Dynamic processing
run('src/Spine_biomech_dyn.m');
```

## Important Configuration Notes

Before running on your machine, update hard-coded paths and file names.

- In `src/Spine_biomech_test.m`:
	- `path = ('E:\XiongSpine');`
	- static trial name (`SpineVST_Buiding_ST.c3d`)
- In `src/Spine_biomech_dyn.m`:
	- dynamic trial name (`Spine_Extention01.c3d`)

If you use files from `sample_data/`, either:

- change the script filenames to match, or
- copy/rename your data files to the expected names.

## Processing Overview

### 1. Marker import and filtering

- Read C3D marker trajectories with `ezc3dRead`
- Convert to per-marker `nFrames x 3`
- Apply low-pass filter (`fc = 6 Hz`, order `4`)

### 2. Spine marker organization

- Parse labels into left, midline, right columns
- Reorder by anatomical level
- Build `SpineCell{level, side}`

### 3. Spline reconstruction

- Fit cubic splines to left/mid/right columns (`spap2`)
- Evaluate dense points along the spine
- Match each vertebral level to nearest midline spline location

### 4. Local coordinate systems per level

- SI axis: tangent of midline spline
- LR direction: right spline point minus left spline point
- AP axis: `cross(SI, LR)`
- ML axis: `cross(AP, SI)`
- Enforce right-handedness (`det(R) > 0`)

### 5. Angle computation

- Thoracic: C7 relative to T12
- Lumbar: T12 relative to S1
- Cardan XYZ extraction used for `[FE, LB, AT]`

## Outputs

From `src/Spine_biomech_test.m`:

- Single-frame marker plot with labels
- Spline reconstruction + coordinate system visualization
- Console printout of thoracic/lumbar angles for the selected frame

From `src/Spine_biomech_dyn.m`:

- Frame-wise thoracic/lumbar angle arrays (`thoracic_deg`, `lumbar_deg`)
- Intersegmental angle arrays (`FE_int`, `LB_int`, `AT_int`)
- Animation of spline geometry and selected LCS axes

## Legacy Helper Functions

Files in `src/funs/` include Robertson-style teaching/demo utilities for lower-limb examples (`CoordSystemsRobertsonGood`, `Find3DAnglesRobertson`, etc.).

For spine replication, the primary entry points are:

- `src/Spine_biomech_test.m`
- `src/Spine_biomech_dyn.m`

## Current Limitations

- Some paths and filenames are still hard-coded.
- `useMask` selection in dynamic script may need explicit level lists depending on MATLAB string behavior.
- No automated unit tests yet.
- No batch runner for multiple trials yet.

## Suggested Next Steps

- Move all hard-coded parameters into a config section or separate config file.
- Add a function-based API (instead of script-only workflow).
- Save outputs to `.mat`/`.csv` for direct figure/table replication.
- Add validation plots against target publication metrics.

## License

MIT. See `LICENSE`.