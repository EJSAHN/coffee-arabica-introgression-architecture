#!/usr/bin/env python3
"""Fail-fast test for the NumPy BLAS/LAPACK stack used by the workflow."""
from __future__ import annotations

import os
import sys

print("Python:", sys.executable, flush=True)
print("CONDA_PREFIX:", os.environ.get("CONDA_PREFIX", ""), flush=True)
print("PATH first entries:", os.environ.get("PATH", "").split(os.pathsep)[:8], flush=True)

import numpy as np

print("NumPy:", np.__version__, flush=True)
rng = np.random.default_rng(20260722)
x = rng.normal(size=(512, 44)).astype(np.float64)

print("Testing matrix multiplication...", flush=True)
gram = x.T @ x
print("  matmul OK; checksum =", float(gram[0, 0]), flush=True)

print("Testing LAPACK eigh...", flush=True)
vals, _ = np.linalg.eigh(gram + np.eye(44, dtype=np.float64) * 1e-8)
print("  eigh OK; leading eigenvalue =", float(vals[-1]), flush=True)

print("Testing LAPACK SVD...", flush=True)
_, singular_values, _ = np.linalg.svd(x[:44, :], full_matrices=False)
print("  SVD OK; leading singular value =", float(singular_values[0]), flush=True)
print("NATIVE_STACK_PREFLIGHT_PASS", flush=True)
