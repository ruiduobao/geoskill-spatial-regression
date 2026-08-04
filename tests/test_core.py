"""Core algorithm tests for spatial-regression."""
import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from conftest import mod


class TestOLS:
    def test_recovers_exact_coefficients(self):
        """无噪声线性数据 OLS 应精确复现系数。"""
        rng = np.random.default_rng(0)
        n = 100
        X = np.column_stack([np.ones(n), rng.normal(0, 1, n), rng.normal(0, 1, n)])
        beta_true = np.array([3.0, 2.0, -1.5])
        y = X @ beta_true
        res = mod.ols_fit(X, y)
        np.testing.assert_allclose(res["beta"], beta_true, atol=1e-8)
        assert res["r2"] > 0.9999

    def test_r2_range(self):
        rng = np.random.default_rng(1)
        n = 80
        X = np.column_stack([np.ones(n), rng.normal(0, 1, n)])
        y = X @ np.array([1.0, 2.0]) + rng.normal(0, 0.5, n)
        res = mod.ols_fit(X, y)
        assert 0 < res["r2"] < 1

    def test_row_mismatch_raises(self):
        with pytest.raises(mod.ValidationError):
            mod.ols_fit(np.ones((5, 2)), np.ones(4))


class TestMoranResiduals:
    def test_random_residuals_near_zero(self):
        rng = np.random.default_rng(2)
        coords = rng.uniform(0, 1, (60, 2))
        W = mod.knn_weights(coords, k=6)
        resid = rng.normal(0, 1, 60)
        res = mod.morans_i_residuals(resid, W)
        assert abs(res["I"]) < 0.2

    def test_autocorrelated_residuals_positive(self):
        rng = np.random.default_rng(3)
        coords = rng.uniform(0, 1, (80, 2))
        W = mod.knn_weights(coords, k=6)
        # 用空间平滑值当残差 → 正自相关
        resid = np.sin(coords[:, 0] * 6) + np.cos(coords[:, 1] * 6)
        res = mod.morans_i_residuals(resid, W)
        assert res["I"] > 0.1


class TestLagrange:
    def test_returns_all_keys(self):
        rng = np.random.default_rng(4)
        coords = rng.uniform(0, 1, (50, 2))
        W = mod.knn_weights(coords, k=6)
        X = np.column_stack([np.ones(50), rng.normal(0, 1, 50)])
        resid = rng.normal(0, 1, 50)
        lm = mod.lagrange_multipliers(resid, W, X)
        for key in ["LM_lag", "LM_error", "RLM_lag", "RLM_error"]:
            assert key in lm
            assert lm[key] >= 0


class TestSLM:
    def test_recovers_rho(self):
        """从 SLM 模拟数据中估计的 rho 应接近真值 0.6。"""
        rng = np.random.default_rng(5)
        coords = rng.uniform(0, 1, (100, 2))
        W = mod.knn_weights(coords, k=6)
        n = 100
        X = np.column_stack([np.ones(n), rng.normal(0, 1, n)])
        beta = np.array([1.0, 2.0])
        rho_true = 0.6
        eps = rng.normal(0, 0.2, n)
        y = np.linalg.solve(np.eye(n) - rho_true * W, X @ beta + eps)
        res = mod.slm_mle(y, X, W, n_grid=150)
        assert abs(res["rho"] - rho_true) < 0.15


class TestSEM:
    def test_recovers_lambda(self):
        """从 SEM 模拟数据中估计的 lambda 应接近真值 0.6。"""
        rng = np.random.default_rng(6)
        coords = rng.uniform(0, 1, (100, 2))
        W = mod.knn_weights(coords, k=6)
        n = 100
        X = np.column_stack([np.ones(n), rng.normal(0, 1, n)])
        beta = np.array([1.0, 2.0])
        lam_true = 0.6
        eps = rng.normal(0, 0.2, n)
        u = np.linalg.solve(np.eye(n) - lam_true * W, eps)
        y = X @ beta + u
        res = mod.sem_mle(y, X, W, n_grid=150)
        assert abs(res["lambda"] - lam_true) < 0.2


class TestSynthetic:
    def test_shapes(self):
        coords, X, y, info = mod.generate_synthetic([116, 39, 117, 40], n_points=50, model="slm")
        assert coords.shape == (50, 2)
        assert X.shape == (50, 3)
        assert y.shape == (50,)
