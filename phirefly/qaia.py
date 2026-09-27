"""MindQuantum QAIA solver wrapper for compressed Phirefly Ising graphs."""

from __future__ import annotations

import inspect
import csv
import time
from dataclasses import dataclass
from pathlib import Path

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
    trace_path: str | None = None


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
    coupling.eliminate_zeros()
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


def solver_spins(solver, n_hyper: int) -> tuple[np.ndarray, np.ndarray]:
    """Return signed hyperread and phaselet spins from a MindQuantum solver."""

    x = solver.x
    if hasattr(x, "detach"):
        x = x.detach().cpu().numpy()
    spins = np.sign(np.asarray(x))
    spins[spins == 0] = 1
    if spins.ndim == 1:
        spins = spins.reshape(-1, 1)
    eta = spins[:n_hyper, :].astype(np.int8, copy=False)
    z = spins[n_hyper:, :].astype(np.int8, copy=False)
    return eta, z


def set_solver_iteration_window(solver, algorithm: str, start: int, stop: int, schedules: dict[str, np.ndarray]) -> None:
    """Limit a MindQuantum solver update to one contiguous schedule window."""

    solver.n_iter = int(stop - start)
    if algorithm in {"BSB", "ASB", "DSB"} and "p" in schedules:
        solver.p = schedules["p"][start:stop]
    elif algorithm == "CFC" and "p" in schedules:
        solver.p = schedules["p"][start:stop]
    elif algorithm == "SimCIM" and "p_list" in schedules:
        solver.p_list = schedules["p_list"][start:stop]


def write_trace_row(
    writer: csv.DictWriter,
    matrix: sparse.csr_matrix,
    phaselet_coupling: sparse.csr_matrix | None,
    eta: np.ndarray,
    z: np.ndarray,
    metadata: dict[str, object],
    iteration: int,
    elapsed_s: float,
    best_objective_seen: float,
) -> float:
    """Write one QAIA trace checkpoint and return updated best objective."""

    scores = score_compressed_samples(matrix, phaselet_coupling, z, eta)
    best_objective = float(np.max(scores)) if scores.size else 0.0
    best_objective_seen = max(best_objective_seen, best_objective)
    writer.writerow(
        {
            **metadata,
            "iteration": int(iteration),
            "elapsed_time_s": f"{elapsed_s:.6f}",
            "energy": f"{-best_objective:.8f}",
            "objective": f"{best_objective:.8f}",
            "best_energy": f"{-best_objective_seen:.8f}",
            "best_objective": f"{best_objective_seen:.8f}",
            "objective_mean": f"{float(np.mean(scores)):.8f}" if scores.size else "0.00000000",
            "objective_sd": f"{float(np.std(scores)):.8f}" if scores.size else "0.00000000",
        }
    )
    return best_objective_seen


def update_solver_with_trace(
    solver,
    matrix: sparse.csr_matrix,
    phaselet_coupling: sparse.csr_matrix | None,
    algorithm: str,
    n_iter: int,
    trace_every: int,
    trace_output: Path,
    trace_metadata: dict[str, object],
    started: float,
) -> None:
    """Run a QAIA solver while writing true checkpoint traces."""

    if algorithm == "LQA":
        raise SystemExit("LQA trace is not supported by MindQuantum chunked update because Adam moments are local to update()")
    trace_output.parent.mkdir(parents=True, exist_ok=True)
    schedules: dict[str, np.ndarray] = {}
    if hasattr(solver, "p"):
        schedules["p"] = np.asarray(solver.p).copy()
    if hasattr(solver, "p_list"):
        schedules["p_list"] = np.asarray(solver.p_list).copy()

    n_hyper = matrix.shape[0]
    fields = [
        "region",
        "solver",
        "seed",
        "batch_size",
        "nodes",
        "couplings",
        "iteration",
        "elapsed_time_s",
        "energy",
        "objective",
        "best_energy",
        "best_objective",
        "objective_mean",
        "objective_sd",
    ]
    checkpoints = list(range(0, int(n_iter), int(trace_every)))
    if not checkpoints or checkpoints[-1] != int(n_iter):
        checkpoints.append(int(n_iter))
    best_seen = -np.inf
    previous = 0
    with trace_output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        eta, z = solver_spins(solver, n_hyper)
        best_seen = write_trace_row(
            writer,
            matrix,
            phaselet_coupling,
            eta,
            z,
            trace_metadata,
            iteration=0,
            elapsed_s=0.0,
            best_objective_seen=best_seen,
        )
        for checkpoint in checkpoints[1:]:
            set_solver_iteration_window(solver, algorithm, previous, checkpoint, schedules)
            solver.update()
            previous = checkpoint
            eta, z = solver_spins(solver, n_hyper)
            best_seen = write_trace_row(
                writer,
                matrix,
                phaselet_coupling,
                eta,
                z,
                trace_metadata,
                iteration=checkpoint,
                elapsed_s=time.monotonic() - started,
                best_objective_seen=best_seen,
            )


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
    trace_every: int = 0,
    trace_output: Path | None = None,
    trace_metadata: dict[str, object] | None = None,
) -> QAIAResult:
    try:
        import mindquantum.algorithm.qaia as qaia
    except ImportError as exc:
        raise SystemExit(
            "MindQuantum QAIA could not be imported. From the source directory run "
            "python -m pip install '.[qaia]' (or '.[gpu]' for GPU support). "
            f"Original error: {exc}"
        ) from exc
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
            trace_path=str(trace_output) if trace_output is not None else None,
        )
    x0 = msf_seed_state(matrix, batch_size, seed, seed_scale, seed_noise)
    kwargs = qaia_kwargs(solver_cls, coupling, x0, n_iter, batch_size, dt, backend)
    started = time.monotonic()
    solver = solver_cls(**kwargs)
    if getattr(solver, "backend", backend) != backend:
        raise RuntimeError(f"Requested {backend}, but solver selected {solver.backend}")
    if trace_every and trace_output is not None:
        metadata = {
            "region": "",
            "solver": algorithm,
            "seed": seed,
            "batch_size": batch_size,
            "nodes": coupling.shape[0],
            "couplings": coupling.nnz // 2,
        }
        if trace_metadata:
            metadata.update(trace_metadata)
        update_solver_with_trace(
            solver,
            matrix,
            phaselet_coupling,
            algorithm=algorithm,
            n_iter=n_iter,
            trace_every=trace_every,
            trace_output=trace_output,
            trace_metadata=metadata,
            started=started,
        )
    else:
        solver.update()
    eta, z = solver_spins(solver, n_hyper)
    runtime_s = time.monotonic() - started
    return QAIAResult(
        phaselet_spins=z,
        hyperread_spins=eta,
        scores=score_compressed_samples(matrix, phaselet_coupling, z, eta),
        coupling_nnz=int(coupling.nnz),
        runtime_s=float(runtime_s),
        trace_path=str(trace_output) if trace_output is not None else None,
    )


__all__ = [
    "QAIAResult",
    "algorithm_dt",
    "build_ising_coupling",
    "msf_seed_state",
    "qaia_kwargs",
    "score_compressed_samples",
    "solve_qaia",
    "solver_spins",
]
