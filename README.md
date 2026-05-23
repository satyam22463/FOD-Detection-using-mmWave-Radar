# FOD Detection on Airport Runways — 77 GHz mmWave Radar

A real-time Foreign Object Debris (FOD) detection system using a **TI AWR1843BOOST 77 GHz mmWave MIMO radar** and MATLAB. Implements the full radar signal processing pipeline to detect small runway hazards, with a comparative evaluation of CFAR detection algorithms in high-clutter environments.

---

## Pipeline

```
[ 77 GHz Radar Front-End ]  -->  Transmits FMCW Chirps
            │
            ▼
  [ IF Signal / ADC Buffer ]  -->  Captures Raw ADC Data
            │
            ▼
  [ 1D Range-FFT Pipeline ]   -->  Resolves Target Distance
            │
            ▼
  [ CFAR Detection Engine  ]  -->  CA / OS / CMLD Thresholding
            │
            ▼
  [ Binary Hazard Map Out  ]  -->  Identifies Runway FOD Locations
```

---

## Technical Specifications

| Parameter | Value |
|---|---|
| Radar Hardware | TI AWR1843BOOST (MIMO) |
| Frequency Band | 77 – 81 GHz |
| Waveform | FMCW (Frequency Modulated Continuous Wave) |
| Core Algorithm | 1D Range-FFT + Adaptive CFAR |
| CFAR Variants | CA-CFAR, OS-CFAR, CMLD-CFAR |
| Processing Tool | MATLAB |

---

## Signal Processing

### Range Estimation

The radar transmits a linear chirp sweeping over bandwidth $B$. The reflected beat signal frequency $f_b$ maps directly to target distance:

$$R = \frac{c \cdot f_b \cdot T_c}{2B}$$

A 1D Range-FFT converts raw ADC samples into frequency-domain peaks, where each spike corresponds to a physical object on the runway.

### Adaptive CFAR Thresholding

A fixed threshold fails in clutter-heavy runway environments. The adaptive threshold is computed as:

$$\tau = \alpha \cdot P_{\text{clutter}}$$

where $\alpha$ controls the false alarm rate $P_{fa}$ and $P_{\text{clutter}}$ is estimated dynamically from neighboring range cells.

```
[ Training Cells ] [ Guard Cells ] [ CUT ] [ Guard Cells ] [ Training Cells ]
        └──────────────────────> [ Noise Estimate ] <──────────────────────┘
```

| CFAR Variant | Method | Best For |
|---|---|---|
| **CA-CFAR** | Arithmetic mean of training cells | Homogeneous clutter |
| **OS-CFAR** | Selects $k$-th ranked training cell | Dense multi-target scenes |
| **CMLD-CFAR** | Censors highest-power cells before averaging | Clutter-edge transitions |

---

## Engineering Challenges

### Target Masking in Dense Debris Clusters
**Problem:** CA-CFAR averages all training cells — a large piece of debris raises the local noise floor and blinds the radar to smaller nearby objects.

**Fix:** Switched to OS-CFAR. Sorting cell magnitudes and selecting the $k$-th percentile isolates large outliers, cleanly resolving closely grouped hazards.

### Clutter-Edge False Alarms
**Problem:** The sharp power contrast at asphalt–grass runway borders causes standard cell-averaging to flag false debris alerts at the edges.

**Fix:** Implemented a split-window configuration with CMLD censoring. Tracking leading and trailing training windows separately lets the threshold follow the step-change smoothly without false spikes.

---

## Key Results

- **Sub-centimetre resolution** — 77 GHz sweeps resolve small metallic and non-metallic debris within targeted range constraints
- **Reduced false alarms** — OS-CFAR and CMLD-CFAR significantly outperform CA-CFAR at runway borders while keeping detection probability $P_d$ stable in multi-target scenarios
---

## Hardware

- **TI AWR1843BOOST** — 77 GHz single-chip MIMO radar (3 TX / 4 RX)
- **MATLAB** — Signal processing, FFT, CFAR engine, visualisation
