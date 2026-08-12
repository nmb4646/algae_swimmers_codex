Batch folder workflow:

```matlab
standardize_trajectory_folder( ...
    './new_raw_folder', ...
    './data_standardized/pico_water', ...
    3.2e-7, ...
    1/240, ...
    0.01, ...
    'DoSmoothing', true, ...
    'SgolayOrder', 3, ...
    'SgolayWindow', 11)
```

This takes every `*.csv` in `./new_raw_folder`, converts `X` and `Y` from pixels to meters, resamples them to `0.01 s`, and writes them into `./data_standardized/pico_water`.

Output:
- standardized trajectory CSVs in your chosen output folder
- a `standardization_summary.csv` in that same folder

Standardized files contain:
- `FrameNumber`
- `Time_s`
- `X`
- `Y`

Units:
- `X`, `Y` are meters.
- `Time_s` is seconds.

If you want a script to read standardized files safely, use:

```matlab
[x, y, t, dt, meta] = load_standardized_trajectory("./data_standardized/pico_water/1.csv");
```

That returns metric coordinates plus the per-file timestep inferred from `Time_s`.

The older manifest workflow is still available through [standardize_trajectories.m](/home/nickbroussinos/Code/algae_swimmers_codex/standardize_trajectories.m:1), but for a new batch from one camera setting you should usually use [standardize_trajectory_folder.m](/home/nickbroussinos/Code/algae_swimmers_codex/standardize_trajectory_folder.m:1).
