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

<!-- ===== 中文原文 (Chinese Original) ===== -->

---
name: geoskill-spatial-regression
description: 'OLS空间诊断+SLM/SEM空间回归模型'
---

# 空间回归分析 | Spatial Regression

空间回归全流程：OLS 拟合 → 残差空间自相关诊断（Moran's I + Lagrange Multiplier）→ 空间滞后模型（SLM）与空间误差模型（SEM）的极大似然估计（网格搜索 + 精确 log 行列式）。

## 核心算法

- OLS 与残差 Moran's I 诊断
- LM-lag / LM-error 检验
- SLM/SEM 极大似然估计

## 依赖

```bash
pip install numpy rasterio scipy geopandas shapely scikit-learn
```

## 使用方法

### 示例 1（合成数据，离线）

```bash
python geoskill-spatial-regression.py --bbox 116.0 39.0 117.0 40.0 --synthetic --output-dir ./out
```

### 示例 2（指定区域 + 静默）

```bash
python geoskill-spatial-regression.py --bbox 121.0 31.0 122.0 32.0 --synthetic --output-dir ./out2 --quiet
```

### 示例 3（真实输入）

```bash
python geoskill-spatial-regression.py --input <你的数据文件> --output-dir ./out3
```

### 示例 4（极小区域边界测试）

```bash
python geoskill-spatial-regression.py --bbox 116.39 39.90 116.40 39.91 --synthetic --output-dir ./out4 --quiet
```

## 输出

| 文件 | 格式 | 说明 |
|---|---|---|
| `regression_stats.json` | GeoTIFF/GeoJSON/JSON | 主产物 |
| `output-manifest.json` | JSON | 运行清单 |

## 数据源 / Source

- 合成模式：本地生成物理一致的模拟数据，无外部数据源。
- 真实模式：读取本地输入文件，无网络请求。

## 隐私声明 / Privacy

- 默认完全离线运行，不发起任何网络请求。
- `--synthetic` 模式不读取任何外部数据。
- 所有计算在本地完成，不上传用户数据。

## License

MIT
