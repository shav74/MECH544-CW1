
---

## Requirements
- MATLAB (R2021a or newer recommended)
- Base MATLAB only (no special toolboxes required)

---

## Quick Start
1. **Add data**: copy all Motive files into `data_raw/`.  
   - Assumes headers on row 8 and data from cell **A8**.
2. **Open MATLAB** and set the **current folder** to the project root.
3. **Run** `main.m` (run all sections).
4. **Outputs**:
   - Figures in `results/figs/`
   - Tables in `results/` (e.g., `avg_speed_and_duration_visible.csv`)
   - Console prints fastest/farthest flights and summary stats.

---

## What the Script Does
1. **Batch import** from `data_raw/` → coerce to numeric → convert **mm → m** (`scale = 0.001`).
2. **Axis remap** so ground plane = **X–Y**, up = **Z**.
3. **Clean** with `fillmissing(...,'nearest')` and `smoothdata(...,'movmean',50)`.
4. **Keep flights only** using simple rules:
   - peak height **Z > 1.5 m**
   - span along throw axis **> 4.5 m**
5. **Trim arm swing** using a speed gate on the ground plane:
   - enter when speed ≥ **20%** of vmax
   - exit when speed ≥ **12%** of vmax
6. **Metrics** (computed on the trimmed in-flight segment):
   - **distance (m)**, **duration (s)**, **average speed (m/s)**
7. **Visuals**:
   - 3D trajectories (filtered & trimmed)
   - Bar charts: average speed, distance, flight time
   - Scatter plots: distance vs speed, distance vs time
   - Histograms: speed, distance, time
   - Landing positions map (ground X–Y)

---

## Key Parameters (tweak if needed)
```matlab
fs_default  = 120;    % sampling rate (Hz)
scale       = 0.001;  % mm -> m
windowSize  = 50;     % smoothing window (samples)
floor_height= 1.5;    % min peak Z (m) to count as a flight
min_length  = 4.5;    % min span along throw axis (m)
enterFrac   = 0.20;   % speed-gate entry (fraction of vmax)
exitFrac    = 0.12;   % speed-gate exit  (fraction of vmax)
