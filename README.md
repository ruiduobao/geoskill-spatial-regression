# Spatial Regression (geoskill-spatial-regression)

> OLS spatial diagnostics + SLM/SEM spatial regression models

---

## 1. Overview

Full spatial regression pipeline: OLS fitting → residual spatial autocorrelation diagnostics (Moran's I + Lagrange Multiplier) → maximum likelihood estimation of the spatial lag model (SLM) and the spatial error model (SEM) (grid search + exact log determinant).

## 2. Features

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

## 7. Methodology

This skill has been methodologically reviewed. See [`REVIEW.md`](./REVIEW.md) for:

- P0/P1/P2 issue counts and verdicts
- Reproduction commands
- Known limitations and edge cases

## 8. License

MIT License. See [`LICENSE`](./LICENSE) for full text.

---

# 空间回归分析（中文版）

> 本部分为完整中文文档，与上方英文部分对应。

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
