---
name: geoskill-spatial-regression
description: 'OLS空间诊断+SLM/SEM空间回归模型'
---

# 空间回归分析 | Spatial Regression

End-to-end spatial regression workflow: OLS fitting → residual spatial autocorrelation diagnostics (Moran's I + Lagrange Multiplier) → maximum likelihood estimation of the spatial lag model (SLM) and spatial error model (SEM) (grid search + exact log-determinant).

## Core Algorithm / 核心算法

- OLS and residual Moran's I diagnostics
- LM-lag / LM-error tests
- SLM/SEM maximum likelihood estimation

## Dependencies / 依赖

```bash
pip install numpy rasterio scipy geopandas shapely scikit-learn
```

## Usage / 使用方法

### Example 1 (synthetic data, offline)

```bash
python geoskill-spatial-regression.py --bbox 116.0 39.0 117.0 40.0 --synthetic --output-dir ./out
```

### Example 2 (custom region + quiet mode)

```bash
python geoskill-spatial-regression.py --bbox 121.0 31.0 122.0 32.0 --synthetic --output-dir ./out2 --quiet
```

### Example 3 (real input)

```bash
python geoskill-spatial-regression.py --input <your data file> --output-dir ./out3
```

### Example 4 (tiny-area boundary test)

```bash
python geoskill-spatial-regression.py --bbox 116.39 39.90 116.40 39.91 --synthetic --output-dir ./out4 --quiet
```

## Output / 输出

| File | Format | Description |
|---|---|---|
| `regression_stats.json` | GeoTIFF/GeoJSON/JSON | Primary output |
| `output-manifest.json` | JSON | Run manifest |

## Data Source / 数据源 / Source

- Synthetic mode: locally generates physically consistent simulated data; no external data sources.
- Real-data mode: reads local input files; no network requests.

## Privacy / 隐私声明 / Privacy

- Runs fully offline by default; makes no network requests.
- `--synthetic` mode reads no external data.
- All computation is done locally; user data is never uploaded.

## License / License

MIT

---

> 中文完整版：需要时再读 [references/SKILL.zh-CN.md](references/SKILL.zh-CN.md)，不要默认全量加载。
