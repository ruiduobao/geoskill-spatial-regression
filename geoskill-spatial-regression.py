#!/usr/bin/env python3
"""spatial-regression — 空间回归分析

对点/面数据拟合回归模型并做空间诊断：
1. 普通最小二乘（OLS）
2. 残差空间自相关诊断（Moran's I + Lagrange Multiplier）
3. 空间滞后模型（SLM/SAR）与空间误差模型（SEM）的极大似然估计

数据源：本地 CSV/GeoJSON，或 --synthetic 生成模拟样本。

隐私声明 / Privacy：
- 默认离线运行，不发起任何网络请求。
- 所有处理在本地完成，不上传用户数据。

Usage:
    python spatial-regression.py --bbox 116 39 117 40 --synthetic --output-dir ./out

License: MIT
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

VERSION = "1.0.0"
SKILL_NAME = "spatial-regression"

try:
    from _geoskill_core.errors import (
        GeoSkillError, UsageError, ValidationError, ProcessError, to_exit_code,
    )
    from _geoskill_core.manifest import OutputManifest, OutputFile
except ImportError:  # pragma: no cover
    class GeoSkillError(Exception):
        def __init__(self, message: str, code: int = 7, kind: str = "EGeo", **kw):
            super().__init__(message)
            self.message, self.code, self.kind = message, code, kind

    class UsageError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=2, kind="EUsage", **k)

    class ValidationError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=6, kind="EValidate", **k)

    class ProcessError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=7, kind="EProcess", **k)

    def to_exit_code(exc):
        return getattr(exc, "code", 7)

    OutputManifest = None
    OutputFile = None


def _utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat().replace("+00:00", "Z")


# ---------------------------------------------------------------------------
# 校验
# ---------------------------------------------------------------------------
def validate_bbox(bbox: Optional[List[float]]) -> None:
    """校验地理 bbox 合法性（W>E 视为跨 180° 不支持）。"""
    if bbox is None:
        return
    w, s, e, n = bbox
    if w > e:
        raise ValidationError(
            "invalid bbox: minLon > maxLon; crossing the 180° meridian "
            "is not supported, split the extent and run twice")
    if s > n:
        raise ValidationError("invalid bbox: minLat > maxLat")
    for lon, lat in ((w, s), (e, n)):
        if not (-180.0 <= lon <= 180.0):
            raise ValidationError(f"invalid bbox: longitude {lon} out of range [-180, 180]")
        if not (-90.0 <= lat <= 90.0):
            raise ValidationError(f"invalid bbox: latitude {lat} out of range [-90, 90]")


def validate_params(n_points: int, knn: int, rho_grid: int) -> None:
    if n_points < 5:
        raise ValidationError(f"--n-points must be >= 5, got {n_points}")
    if knn < 1:
        raise ValidationError(f"--knn must be >= 1, got {knn}")
    if rho_grid < 10:
        raise ValidationError(f"--rho-grid must be >= 10, got {rho_grid}")


# ---------------------------------------------------------------------------
# 空间权重
# ---------------------------------------------------------------------------
def knn_weights(coords: np.ndarray, k: int = 6) -> np.ndarray:
    """k 最近邻行标准化权重矩阵。"""
    from scipy.spatial import cKDTree
    n = coords.shape[0]
    k = min(k, n - 1)
    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=k + 1)
    W = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        neigh = idx[i][1:]
        W[i, neigh] = 1.0 / k
    return W


# ---------------------------------------------------------------------------
# 核心算法：OLS
# ---------------------------------------------------------------------------
def ols_fit(X: np.ndarray, y: np.ndarray) -> Dict[str, Any]:
    """普通最小二乘。X 应已含截距列。返回系数、R²、残差等。"""
    X = np.asarray(X, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64).ravel()
    n, k = X.shape
    if n != y.shape[0]:
        raise ValidationError("X and y row count mismatch")
    if n <= k:
        raise ValidationError("not enough observations for OLS")
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    resid = y - yhat
    sse = float(resid @ resid)
    sst = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - sse / sst if sst > 0 else 0.0
    sigma2 = sse / (n - k)
    # 系数方差
    try:
        xtx_inv = np.linalg.inv(X.T @ X)
        se = np.sqrt(np.clip(np.diag(xtx_inv) * sigma2, 0, None))
    except np.linalg.LinAlgError:
        se = np.full(k, np.nan)
    tvals = beta / np.where(se > 0, se, np.nan)
    return {"beta": beta, "resid": resid, "yhat": yhat, "r2": float(r2),
            "sigma2": float(sigma2), "std_err": se, "t_values": tvals,
            "n": n, "k": k}


def morans_i_residuals(resid: np.ndarray, W: np.ndarray) -> Dict[str, float]:
    """残差 Global Moran's I（标准化残差，Cliff & Ord 1981 方差）。"""
    z = resid - resid.mean()
    n = z.shape[0]
    S0 = float(W.sum())
    den = float((z ** 2).sum())
    if S0 <= 1e-12 or den <= 1e-300:
        return {"I": 0.0, "z_score": 0.0, "p_value": 1.0}
    I = (n / S0) * (float(z @ W @ z) / den)
    EI = -1.0 / (n - 1)
    # Cliff & Ord (1981): S1 = 1/2 * sum((w_ij + w_ji)^2), S2 = sum((w_i. + w_.i)^2)
    S1 = 0.5 * float(((W + W.T) ** 2).sum())
    row = W.sum(axis=1)
    col = W.sum(axis=0)
    S2 = float(((row + col) ** 2).sum())
    VarI = (n * n * S1 - n * S2 + 3 * S0 * S0) / (S0 * S0 * (n * n - 1)) - EI * EI
    VarI = max(VarI, 1e-12)
    zscore = (I - EI) / np.sqrt(VarI)
    from scipy.stats import norm
    p = 2.0 * (1.0 - norm.cdf(abs(zscore)))
    return {"I": float(I), "z_score": float(zscore), "p_value": float(p)}


def lagrange_multipliers(resid: np.ndarray, W: np.ndarray, X: np.ndarray,
                         y: Optional[np.ndarray] = None) -> Dict[str, float]:
    """LM-lag / LM-error 及其稳健形式（Anselin 1988, ch.6）。

    LM_err  = (e'We/s2)^2 / T,            T = tr((W+W')W)
    LM_lag  = (e'Wy/s2)^2 / nJ,           nJ = T + (WXb)' M (WXb)/s2, M = I - X(X'X)^-1 X'
    RLM_lag = (e'Wy/s2 - e'We/s2)^2 / (nJ - T)
    RLM_err = (e'We/s2 - e'Wy/s2)^2 / (T - T^2/nJ)
    """
    n = X.shape[0]
    e = resid - resid.mean()
    s2 = float((e ** 2).sum() / n)
    eWe = float(e @ (W @ e)) / s2
    T = float(np.trace((W + W.T) @ W))
    LM_error = (eWe ** 2) / T if T > 1e-12 else 0.0
    if y is not None:
        eWy = float(e @ (W @ y.ravel())) / s2
        beta, *_ = np.linalg.lstsq(X, y.ravel(), rcond=None)
        WXb = W @ (X @ beta)
        M = np.eye(n) - X @ np.linalg.pinv(X)
        nJ = T + float(WXb @ M @ WXb) / s2 if s2 > 1e-300 else T
    else:  # API 兼容：无 y 时退化为 e'Wy ≈ e'We（CLI 路径总是传入 y）
        eWy = eWe
        nJ = T
    if nJ <= 1e-12:
        LM_lag, RLM_lag, RLM_error = 0.0, 0.0, 0.0
    else:
        diff = eWy - eWe
        LM_lag = (eWy ** 2) / nJ
        RLM_lag = (diff ** 2) / max(nJ - T, 1e-12)
        RLM_error = (diff ** 2) / max(T - T * T / nJ, 1e-12)
    return {"LM_lag": float(LM_lag), "LM_error": float(LM_error),
            "RLM_lag": float(RLM_lag), "RLM_error": float(RLM_error)}


# ---------------------------------------------------------------------------
# 核心算法：SLM / SEM 极大似然
# ---------------------------------------------------------------------------
def _logdet_spd(A: np.ndarray) -> float:
    sign, logdet = np.linalg.slogdet(A)
    return float(logdet) if sign > 0 else -np.inf


def slm_mle(y: np.ndarray, X: np.ndarray, W: np.ndarray,
            n_grid: int = 200) -> Dict[str, Any]:
    """空间滞后模型 y = rho*Wy + X*beta + eps 的极大似然（网格搜索 rho）。"""
    y = y.ravel().astype(np.float64)
    X = X.astype(np.float64)
    n = y.shape[0]
    Wy = W @ y
    # 特征值用于精确 log|I - rho W|
    evals = np.linalg.eigvals(W)
    evals = evals[np.abs(evals.imag) < 1e-9].real
    rho_min = 1.0 / evals.min() + 1e-3 if evals.min() < 0 else -0.999
    rho_max = 1.0 / evals.max() - 1e-3 if evals.max() > 0 else 0.999
    rho_min, rho_max = max(rho_min, -0.999), min(rho_max, 0.999)
    grid = np.linspace(rho_min, rho_max, n_grid)
    best = {"ll": -np.inf, "rho": 0.0}
    for rho in grid:
        A = y - rho * Wy
        beta, *_ = np.linalg.lstsq(X, A, rcond=None)
        e = A - X @ beta
        s2 = (e @ e) / n
        logdet = np.log(np.abs(np.real(np.linalg.det(np.eye(n) - rho * W))) + 1e-300)
        ll = -n / 2 * (np.log(2 * np.pi) + np.log(s2) + 1) + logdet
        if ll > best["ll"]:
            best = {"ll": float(ll), "rho": float(rho), "beta": beta, "sigma2": float(s2)}
    return best


def sem_mle(y: np.ndarray, X: np.ndarray, W: np.ndarray,
            n_grid: int = 200) -> Dict[str, Any]:
    """空间误差模型 y = X*beta + u, u = lambda*Wu + eps 的极大似然（网格搜索 lambda）。"""
    y = y.ravel().astype(np.float64)
    X = X.astype(np.float64)
    n = y.shape[0]
    evals = np.linalg.eigvals(W)
    evals = evals[np.abs(evals.imag) < 1e-9].real
    lam_min = 1.0 / evals.min() + 1e-3 if evals.min() < 0 else -0.999
    lam_max = 1.0 / evals.max() - 1e-3 if evals.max() > 0 else 0.999
    lam_min, lam_max = max(lam_min, -0.999), min(lam_max, 0.999)
    grid = np.linspace(lam_min, lam_max, n_grid)
    beta_ols, *_ = np.linalg.lstsq(X, y, rcond=None)
    e = y - X @ beta_ols
    We = W @ e
    best = {"ll": -np.inf, "lambda": 0.0}
    for lam in grid:
        e_star = e - lam * We
        # 广义最小二乘：用 (I - lam W) 变换
        A = np.eye(n) - lam * W
        yt = A @ y
        Xt = A @ X
        beta, *_ = np.linalg.lstsq(Xt, yt, rcond=None)
        resid = yt - Xt @ beta
        s2 = (resid @ resid) / n
        logdet = np.log(np.abs(np.real(np.linalg.det(A))) + 1e-300)
        ll = -n / 2 * (np.log(2 * np.pi) + np.log(s2) + 1) + logdet
        if ll > best["ll"]:
            best = {"ll": float(ll), "lambda": float(lam), "beta": beta, "sigma2": float(s2)}
    return best


# ---------------------------------------------------------------------------
# 合成数据
# ---------------------------------------------------------------------------
def generate_synthetic(bbox: List[float], n_points: int = 120, model: str = "slm",
                       seed: int = 42) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
    """生成带空间结构的模拟样本。

    Returns (coords, X (含截距), y, info)。
    """
    rng = np.random.default_rng(seed)
    w, s, e, n = bbox
    coords = rng.uniform([w, s], [e, n], (n_points, 2))
    x1 = rng.normal(0, 1, n_points)
    x2 = rng.normal(0, 1, n_points)
    X = np.column_stack([np.ones(n_points), x1, x2])
    beta_true = np.array([2.0, 1.5, -1.0])
    W = knn_weights(coords, k=6)
    if model == "slm":
        # y = rho Wy + X beta + eps → 解 (I - rho W) y = X beta + eps
        eps = rng.normal(0, 0.3, n_points)
        rho = 0.6
        y = np.linalg.solve(np.eye(n_points) - rho * W, X @ beta_true + eps)
    elif model == "sem":
        eps = rng.normal(0, 0.3, n_points)
        lam = 0.6
        u = np.linalg.solve(np.eye(n_points) - lam * W, eps)
        y = X @ beta_true + u
    else:  # ols (无空间结构)
        y = X @ beta_true + rng.normal(0, 0.5, n_points)
    info = {"n_points": n_points, "model": model, "beta_true": beta_true.tolist()}
    return coords, X, y, info


# ---------------------------------------------------------------------------
# GeoTIFF I/O（输出系数栅格用）/ Manifest
# ---------------------------------------------------------------------------
def write_geotiff(path: str, array: np.ndarray, bbox: List[float], nodata: float = -9999.0) -> None:
    import rasterio
    from rasterio.transform import from_bounds

    if array.ndim == 2:
        array = array[np.newaxis, ...]
    if array.ndim != 3:
        raise ValidationError("write_geotiff expects a 2D or 3D array")
    nb, hh, ww = array.shape
    transform = from_bounds(bbox[0], bbox[1], bbox[2], bbox[3], ww, hh)
    profile = {
        "driver": "GTiff", "height": hh, "width": ww, "count": nb,
        "dtype": "float32", "crs": "EPSG:4326", "transform": transform,
        "nodata": nodata, "compress": "deflate",
    }
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with rasterio.open(path, "w", **profile) as dst:
        for b in range(nb):
            dst.write(array[b].astype("float32"), b + 1)


def write_manifest(
    output_dir: str, args: argparse.Namespace, outputs: List[Dict[str, Any]],
    qa: Dict[str, Any], started_at: str, exit_code: int, bbox: List[float],
) -> Optional[str]:
    if OutputManifest is None:
        return None
    cmd = " ".join([SKILL_NAME] + sys.argv[1:])
    man = OutputManifest(
        skill=SKILL_NAME, skill_version=VERSION, command=cmd,
        started_at=started_at, finished_at=_utc_now(), exit_code=exit_code,
        inputs={"input": getattr(args, "input", None), "synthetic": bool(getattr(args, "synthetic", False))},
        outputs=[OutputFile(**o) for o in outputs],
        qa=qa,
        software={"python": sys.version.split()[0], "skill": SKILL_NAME},
    )
    path = os.path.join(output_dir, "output-manifest.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(man.to_dict(), f, ensure_ascii=False, indent=2, default=str)
    return path


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def process(args: argparse.Namespace) -> int:
    started_at = _utc_now()
    output_dir = args.output_dir
    bbox = list(args.bbox) if args.bbox else None
    validate_params(args.n_points, args.knn, args.rho_grid)
    validate_bbox(bbox)

    if args.input and not args.synthetic:
        if not os.path.exists(args.input):
            raise UsageError(f"input not found: {args.input}", path=args.input)
        try:
            coords, X, y = _load_csv(args.input)
        except (ValueError, ValidationError) as exc:
            if isinstance(exc, ValidationError):
                raise
            raise ValidationError(
                f"failed to read input '{args.input}' (must be numeric CSV "
                f"with x,y,dep[,ind*]): {exc}") from exc
        if bbox is None:
            bbox = [float(coords[:, 0].min()), float(coords[:, 1].min()),
                    float(coords[:, 0].max()), float(coords[:, 1].max())]
        source_note = args.input
    else:
        if bbox is None:
            raise UsageError("provide --bbox (synthetic mode) or --input <csv>")
        coords, X, y, _ = generate_synthetic(bbox, n_points=args.n_points, model=args.model)
        source_note = "synthetic"

    if coords.shape[0] < X.shape[1] + 2:
        raise ValidationError("too few observations")

    W = knn_weights(coords, k=args.knn)
    ols = ols_fit(X, y)
    diag = morans_i_residuals(ols["resid"], W)
    lm = lagrange_multipliers(ols["resid"], W, X, y)
    slm = slm_mle(y, X, W, n_grid=args.rho_grid)
    sem = sem_mle(y, X, W, n_grid=args.rho_grid)

    os.makedirs(output_dir, exist_ok=True)
    stats_path = os.path.join(output_dir, "regression_stats.json")
    summary = {
        "ols": {"beta": ols["beta"].tolist(), "r2": ols["r2"],
                "std_err": ols["std_err"].tolist(), "t_values": ols["t_values"].tolist()},
        "diagnostics": {"morans_i": diag, "lagrange_multipliers": lm},
        "slm": {"rho": slm["rho"], "beta": np.asarray(slm["beta"]).tolist(),
                "sigma2": slm["sigma2"], "log_likelihood": slm["ll"]},
        "sem": {"lambda": sem["lambda"], "beta": np.asarray(sem["beta"]).tolist(),
                "sigma2": sem["sigma2"], "log_likelihood": sem["ll"]},
    }
    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    qa = {
        "source": source_note,
        "n_obs": int(coords.shape[0]),
        "ols_r2": ols["r2"],
        "morans_i_resid": diag["I"],
        "slm_rho": slm["rho"],
        "sem_lambda": sem["lambda"],
    }
    outputs = [{"path": stats_path, "kind": "json"}]
    man_path = write_manifest(output_dir, args, outputs, qa, started_at, 0, bbox)

    if not args.quiet:
        print(f"[{SKILL_NAME}] source: {source_note}  n={coords.shape[0]}")
        print(f"[{SKILL_NAME}] OLS R²={ols['r2']:.4f}  beta={np.round(ols['beta'], 3).tolist()}")
        print(f"[{SKILL_NAME}] Moran's I(resid)={diag['I']:.4f}  p={diag['p_value']:.4f}")
        print(f"[{SKILL_NAME}] SLM rho={slm['rho']:.3f}  SEM lambda={sem['lambda']:.3f}")
        print(f"[{SKILL_NAME}] output: {stats_path}")
        if man_path:
            print(f"[{SKILL_NAME}] manifest: {man_path}")
    return 0


def _load_csv(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """读取 CSV：列 x, y, dep, ind1, ind2, ... → (coords, X含截距, y)。"""
    import csv
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    if len(rows) < 2:
        raise ValidationError("CSV has no data rows")
    header = [h.strip().lower() for h in rows[0]]
    data = np.array([[float(v) for v in r] for r in rows[1:] if r], dtype=np.float64)
    col = {h: i for i, h in enumerate(header)}
    if "x" not in col or "y" not in col or "dep" not in col:
        raise ValidationError("CSV must contain columns: x, y, dep")
    coords = data[:, [col["x"], col["y"]]]
    y = data[:, col["dep"]]
    ind_cols = [i for h, i in col.items() if h.startswith("ind")]
    inds = data[:, ind_cols] if ind_cols else np.zeros((data.shape[0], 0))
    X = np.column_stack([np.ones(data.shape[0]), inds])
    return coords, X, y


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=SKILL_NAME,
        description="Spatial regression: OLS diagnostics + SLM/SEM maximum-likelihood estimation.",
    )
    p.add_argument("--bbox", nargs=4, type=float, metavar=("W", "S", "E", "N"),
                   help="geographic extent [minLon minLat maxLon maxLat]")
    p.add_argument("--input", help="input CSV (x,y,dep,ind...)")
    p.add_argument("--model", default="slm", choices=["slm", "sem", "ols"],
                   help="synthetic data-generating model (default: slm)")
    p.add_argument("--n-points", type=int, default=120, help="synthetic sample size (default: 120)")
    p.add_argument("--knn", type=int, default=6, help="k-nearest neighbors for W (default: 6)")
    p.add_argument("--rho-grid", type=int, default=100, help="MLE grid resolution (default: 100)")
    p.add_argument("--synthetic", action="store_true", help="use synthetic data")
    p.add_argument("--output-dir", default="./output", help="output directory")
    p.add_argument("--quiet", action="store_true", help="suppress output")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return process(args)
    except GeoSkillError as exc:
        print(f"[{SKILL_NAME}] ERROR [{exc.kind}] {exc.message}", file=sys.stderr)
        return exc.code
    except KeyboardInterrupt:
        return 130
    except Exception as exc:  # noqa: BLE001
        print(f"[{SKILL_NAME}] ERROR {exc}", file=sys.stderr)
        return to_exit_code(exc)


if __name__ == "__main__":
    sys.exit(main())
