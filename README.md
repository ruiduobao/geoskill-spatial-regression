# Spatial Regression (geoskill-spatial-regression)

> OLS spatial diagnostics + SLM/SEM spatial regression models

---

## 1. Overview

Full spatial regression pipeline: OLS fitting → residual spatial autocorrelation diagnostics (Moran's I + Lagrange Multiplier) → maximum likelihood estimation of the spatial lag model (SLM) and the spatial error model (SEM) (grid search + exact log determinant).

## 3. Quick Start

```bash
pip install -r requirements.txt
python geoskill-spatial-regression.py --bbox 116 39 117 40 --synthetic --output-dir ./out
```

## 4. CLI Parameters

Run `python <skill>.py --help` for the full list. Common parameters:

| Parameter | Type | Description |
|---|---|---|
| `--bbox` | `float[4]` | WGS84 bounding box `min_lon min_lat max_lon max_lat` |
| `--input` | `path` | Local input file (GeoJSON/GeoTIFF/etc.) |
| `--output-dir` | `path` | Output directory (default `./output`) |
| `--synthetic` | `flag` | Use synthetic data instead of real input |
| `--quiet` | `flag` | Suppress non-essential stdout |

## 5. Input / Output

| File | Format | Description |
|---|---|---|
| `regression_stats.json` | GeoTIFF/GeoJSON/JSON | Primary output |
| `output-manifest.json` | JSON | Run manifest |

## 6. Technical Principle

- OLS and residual Moran's I diagnostics
- LM-lag / LM-error tests
- SLM/SEM maximum likelihood estimation

## 8. License

MIT License. See [`LICENSE`](./LICENSE) for full text.

---

> 中文完整版：[references/SKILL.zh-CN.md](references/SKILL.zh-CN.md)
