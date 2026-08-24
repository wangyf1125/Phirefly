"""MindQuantum QAIA solver wrapper for compressed Phirefly Ising graphs."""

from __future__ import annotations

import inspect
import time
from dataclasses import dataclass

import numpy as np
from scipy import sparse

from .core.spins import sign_keep_zero


@dataclass(frozen=True)
class QAIAResult:
    phaselet_spins: np.ndarray
    hyperread_spins: np.ndarray
    scores: np.ndarray
    coupling_nnz: int
    runtime_s: float


def algorithm_dt(algorithm: str) -> str:
    if algorithm != "CFC":
        raise SystemExit("Phirefly public release supports only CFC QAIA for the final method")
    return "0.1"


def build_ising_coupling(matrix: sparse.csr_matrix, phaselet_coupling: sparse.csr_matrix | None) -> sparse.csr_matrix:
    """Convert hyperread-phaselet evidence into one symmetric Ising coupling matrix."""

    n_hyper, n_phaselets = matrix.shape
    n_nodes = n_hyper + n_phaselets
    coo = matrix.tocoo()
    rows = np.concatenate([coo.row, coo.col + n_hyper])
    cols = np.concatenate([coo.col + n_hyper, coo.row])
    data = np.concatenate([coo.data, coo.data]).astype(np.float64, copy=False)
    if phaselet_coupling is not None and phaselet_coupling.nnz:
        phaselet_coupling = phaselet_coupling.tocoo()
        rows = np.concatenate([rows, phaselet_coupling.row + n_hyper])
        cols = np.concatenate([cols, phaselet_coupling.col + n_hyper])
        data = np.concatenate([data, phaselet_coupling.data.astype(np.float64, copy=False)])
    coupling = sparse.coo_matrix((data, (rows, cols)), shape=(n_nodes, n_nodes)).tocsr()
    coupling.sum_duplicates()
    return coupling


def msf_seed_state(
    matrix: sparse.csr_matrix,
    batch_size: int,
    seed: int,
    seed_scale: float,
    seed_noise: float,
) -> np.ndarray:
    """Build QAIA continuous initial state from the all-positive phaselet backbone."""

    n_hyper, n_phaselets = matrix.shape
    z_seed = np.ones((n_phaselets, batch_size), dtype=np.int8)
    eta_seed = sign_keep_zero(matrix @ z_seed)
    seed_spins = np.vstack([eta_seed, z_seed]).astype(np.float64, copy=False)
    rng = np.random.default_rng(seed)
    x0 = float(seed_scale) * seed_spins
    if seed_noise:
        x0 = x0 + float(seed_noise) * rng.standard_normal(size=x0.shape)
    return x0


def qaia_kwargs(solver_cls, coupling: sparse.csr_matrix, x0: np.ndarray, n_iter: int, batch_size: int, dt: float, backend: str) -> dict:
    solver_params = inspect.signature(solver_cls).parameters
    kwargs = {"J": coupling, "h": None, "x": x0, "n_iter": int(n_iter), "batch_size": int(batch_size)}
    if "dt" in solver_params:
        kwargs["dt"] = float(dt)
    if "xi" in solver_params:
        abs_row_sums = np.asarray(abs(coupling).sum(axis=1)).ravel()
        max_abs_row_sum = float(abs_row_sums.max()) if abs_row_sums.size else 0.0
        kwargs["xi"] = 1.0 / max_abs_row_sum if max_abs_row_sum > 0 else 1.0
    if backend:
        kwargs["backend"] = backend
    return kwargs


def score_compressed_samples(
    matrix: sparse.csr_matrix,
    phaselet_coupling: sparse.csr_matrix | None,
    z: np.ndarray,
    eta: np.ndarray,
) -> np.ndarray:
    scores = np.asarray((matrix @ z) * eta, dtype=np.float64).sum(axis=0)
    if phaselet_coupling is not None and phaselet_coupling.nnz:
        scores = scores + 0.5 * np.asarray((phaselet_coupling @ z) * z, dtype=np.float64).sum(axis=0)
    return scores


def solve_qaia(
    matrix: sparse.csr_matrix,
    phaselet_coupling: sparse.csr_matrix | None,
    algorithm: str,
    backend: str,
    batch_size: int,
    n_iter: int,
    dt: float,
    seed: int,
    seed_scale: float,
    seed_noise: float,
) -> QAIAResult:
    try:
        import mindquantum.algorithm.qaia as qaia
    except Exception as exc:
        raise SystemExit("MindQuantum QAIA is unavailable in this environment") from exc
    solver_cls = getattr(qaia, algorithm, None)
    if solver_cls is None:
        raise SystemExit(f"MindQuantum QAIA algorithm {algorithm!r} is unavailable")

    n_hyper, n_phaselets = matrix.shape
    coupling = build_ising_coupling(matrix, phaselet_coupling)
    if coupling.nnz == 0:
        z = np.ones((n_phaselets, int(batch_size)), dtype=np.int8)
        eta = np.ones((n_hyper, int(batch_size)), dtype=np.int8)
        return QAIAResult(
            phaselet_spins=z,
            hyperread_spins=eta,
            scores=np.zeros(int(batch_size), dtype=np.float64),
            coupling_nnz=0,
            runtime_s=0.0,
        )
    x0 = msf_seed_state(matrix, batch_size, seed, seed_scale, seed_noise)
    kwargs = qaia_kwargs(solver_cls, coupling, x0, n_iter, batch_size, dt, backend)
    started = time.monotonic()
    try:
        solver = solver_cls(**kwargs)
    except TypeError:
        kwargs.pop("backend", None)
        solver = solver_cls(**kwargs)
    solver.update()
    runtime_s = time.monotonic() - started
    x = solver.x
    if hasattr(x, "detach"):
        x = x.detach().cpu().numpy()
    spins = np.sign(np.asarray(x))
    spins[spins == 0] = 1
    if spins.ndim == 1:
        spins = spins.reshape(-1, 1)
    eta = spins[:n_hyper, :].astype(np.int8, copy=False)
    z = spins[n_hyper:, :].astype(np.int8, copy=False)
    return QAIAResult(
        phaselet_spins=z,
        hyperread_spins=eta,
        scores=score_compressed_samples(matrix, phaselet_coupling, z, eta),
        coupling_nnz=int(coupling.nnz),
        runtime_s=float(runtime_s),
    )


__all__ = [
    "QAIAResult",
    "algorithm_dt",
    "build_ising_coupling",
    "msf_seed_state",
    "qaia_kwargs",
    "score_compressed_samples",
    "solve_qaia",
]
