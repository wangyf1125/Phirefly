"""Risk scoring and hard/soft/drop classification for SNP parity edges."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PhaseEdge:
    snp_i: int
    snp_j: int
    score: float
    support: float
    count: int
    confidence: float
    conflict_ratio: float
    domination_score: float
    distance_bp: int
    risk: float
    expected_sign: int
    observed_sign: int


def phase_edge_risk(
    confidence: float,
    count: int,
    conflict_ratio: float,
    distance_bp: int,
    domination_score: float,
    risk_alpha_conf: float,
    risk_beta_count: float,
    risk_gamma_conflict: float,
    risk_delta_distance: float,
    risk_epsilon_domination: float,
    risk_distance_norm_bp: int,
) -> float:
    count_penalty = 1.0 / math.log1p(count) if count > 0 else 1.0
    distance_norm = min(float(distance_bp) / max(1.0, float(risk_distance_norm_bp)), 1.0)
    return (
        float(risk_alpha_conf) * (1.0 - confidence)
        + float(risk_beta_count) * count_penalty
        + float(risk_gamma_conflict) * conflict_ratio
        + float(risk_delta_distance) * distance_norm
        + float(risk_epsilon_domination) * domination_score
    )


def make_phase_edge(
    i: int,
    j: int,
    metrics: tuple[float, float, int, float],
    positions: np.ndarray,
    tau_vec: np.ndarray,
    risk_alpha_conf: float,
    risk_beta_count: float,
    risk_gamma_conflict: float,
    risk_delta_distance: float,
    risk_epsilon_domination: float,
    risk_distance_norm_bp: int,
) -> PhaseEdge:
    score, support, count, max_abs_contribution = metrics
    distance_bp = abs(int(positions[j]) - int(positions[i]))
    confidence = abs(score) / support if support > 0 else 0.0
    conflict_ratio = 1.0 - confidence
    domination_score = max_abs_contribution / support if support > 0 else 1.0
    risk = phase_edge_risk(
        confidence,
        count,
        conflict_ratio,
        distance_bp,
        domination_score,
        risk_alpha_conf,
        risk_beta_count,
        risk_gamma_conflict,
        risk_delta_distance,
        risk_epsilon_domination,
        risk_distance_norm_bp,
    )
    return PhaseEdge(
        snp_i=int(i),
        snp_j=int(j),
        score=float(score),
        support=float(support),
        count=int(count),
        confidence=float(confidence),
        conflict_ratio=float(conflict_ratio),
        domination_score=float(domination_score),
        distance_bp=int(distance_bp),
        risk=float(risk),
        expected_sign=int(tau_vec[i]) * int(tau_vec[j]),
        observed_sign=1 if score > 0 else -1 if score < 0 else 0,
    )


def hard_edge_passes(
    edge: PhaseEdge,
    min_abs_support: float,
    min_count: int,
    min_confidence: float,
    bridge_mode: str,
    risk_threshold: float,
) -> tuple[bool, bool, bool]:
    support_ok = edge.support >= min_abs_support and edge.count >= min_count
    filter_ok = edge.confidence >= min_confidence if bridge_mode == "confidence" else edge.risk <= risk_threshold
    sign_ok = edge.score * edge.expected_sign > 0
    return support_ok, filter_ok, sign_ok


def soft_edge_passes(
    edge: PhaseEdge,
    soft_min_abs_support: float,
    soft_min_count: int,
    soft_min_confidence: float,
    soft_risk_threshold: float,
) -> bool:
    return (
        edge.support >= soft_min_abs_support
        and edge.count >= soft_min_count
        and edge.confidence >= soft_min_confidence
        and edge.risk <= soft_risk_threshold
    )


def classify_phase_edge(
    edge: PhaseEdge,
    min_abs_support: float,
    min_count: int,
    min_confidence: float,
    bridge_mode: str,
    risk_threshold: float,
    soft_bridge_mode: str,
    soft_min_abs_support: float,
    soft_min_count: int,
    soft_min_confidence: float,
    soft_risk_threshold: float,
) -> str:
    support_ok, filter_ok, sign_ok = hard_edge_passes(edge, min_abs_support, min_count, min_confidence, bridge_mode, risk_threshold)
    if support_ok and filter_ok and sign_ok:
        return "hard"
    if soft_bridge_mode != "none" and soft_edge_passes(edge, soft_min_abs_support, soft_min_count, soft_min_confidence, soft_risk_threshold):
        return "soft" if sign_ok else "soft_conflicting_tau"
    return "drop"


def drop_status(
    edge: PhaseEdge,
    min_abs_support: float,
    min_count: int,
    min_confidence: float,
    bridge_mode: str,
    risk_threshold: float,
) -> str:
    support_ok, _filter_ok, sign_ok = hard_edge_passes(edge, min_abs_support, min_count, min_confidence, bridge_mode, risk_threshold)
    if not support_ok:
        return "drop_weak_hard_support"
    if bridge_mode == "confidence" and edge.confidence < min_confidence:
        return "drop_weak_hard_confidence"
    if bridge_mode == "risk" and edge.risk > risk_threshold:
        return "drop_high_hard_risk"
    return "drop_conflicting_tau" if not sign_ok else "drop_noisy_soft"


def edge_row(edge: PhaseEdge, snp_ids: list[str], positions: np.ndarray, status: str) -> dict[str, object]:
    return {
        "snp_i": edge.snp_i,
        "snp_j": edge.snp_j,
        "chrom": snp_ids[edge.snp_i].split(":", 1)[0],
        "pos_i": int(positions[edge.snp_i]),
        "pos_j": int(positions[edge.snp_j]),
        "distance_bp": edge.distance_bp,
        "score": f"{edge.score:.8f}",
        "support": f"{edge.support:.8f}",
        "count": edge.count,
        "confidence": f"{edge.confidence:.8f}",
        "conflict_ratio": f"{edge.conflict_ratio:.8f}",
        "domination_score": f"{edge.domination_score:.8f}",
        "risk": f"{edge.risk:.8f}",
        "expected_sign": edge.expected_sign,
        "observed_sign": edge.observed_sign,
        "status": status,
    }


def soft_edge_rows(soft_edges: list[PhaseEdge]) -> list[dict[str, object]]:
    return [
        {
            "snp_i": edge.snp_i,
            "snp_j": edge.snp_j,
            "score": edge.score,
            "support": edge.support,
            "count": edge.count,
            "confidence": edge.confidence,
            "risk": edge.risk,
            "expected_sign": edge.expected_sign,
            "observed_sign": edge.observed_sign,
        }
        for edge in soft_edges
    ]


def classify_parity_edges(
    parity_edges: dict[tuple[int, int], tuple[float, float, int, float]],
    snp_ids: list[str],
    tau_vec: np.ndarray,
    positions: np.ndarray,
    max_distance_bp: int,
    min_abs_support: float,
    min_count: int,
    min_confidence: float,
    bridge_mode: str,
    risk_threshold: float,
    risk_alpha_conf: float,
    risk_beta_count: float,
    risk_gamma_conflict: float,
    risk_delta_distance: float,
    risk_epsilon_domination: float,
    risk_distance_norm_bp: int,
    soft_bridge_mode: str,
    soft_min_abs_support: float,
    soft_min_count: int,
    soft_min_confidence: float,
    soft_risk_threshold: float,
) -> tuple[list[PhaseEdge], list[PhaseEdge], list[dict[str, object]], dict[str, object]]:
    hard_edges: list[PhaseEdge] = []
    soft_edges: list[PhaseEdge] = []
    edge_rows: list[dict[str, object]] = []
    conflicting = 0
    weak = 0

    for (i, j), metrics in parity_edges.items():
        edge = make_phase_edge(
            i,
            j,
            metrics,
            positions,
            tau_vec,
            risk_alpha_conf,
            risk_beta_count,
            risk_gamma_conflict,
            risk_delta_distance,
            risk_epsilon_domination,
            risk_distance_norm_bp,
        )
        if edge.distance_bp > max_distance_bp:
            continue
        support_ok, filter_ok, sign_ok = hard_edge_passes(edge, min_abs_support, min_count, min_confidence, bridge_mode, risk_threshold)
        edge_class = classify_phase_edge(
            edge,
            min_abs_support,
            min_count,
            min_confidence,
            bridge_mode,
            risk_threshold,
            soft_bridge_mode,
            soft_min_abs_support,
            soft_min_count,
            soft_min_confidence,
            soft_risk_threshold,
        )
        if edge_class == "hard":
            hard_edges.append(edge)
            continue
        if support_ok and filter_ok and not sign_ok:
            conflicting += 1
        if edge_class.startswith("soft"):
            soft_edges.append(edge)
            status = "soft_bridge" if edge_class == "soft" else "soft_bridge_conflicting_tau"
            edge_rows.append(edge_row(edge, snp_ids, positions, status))
            continue
        status = drop_status(edge, min_abs_support, min_count, min_confidence, bridge_mode, risk_threshold)
        if status in {"drop_weak_hard_support", "drop_weak_hard_confidence", "drop_high_hard_risk"}:
            weak += 1
        edge_rows.append(edge_row(edge, snp_ids, positions, status))

    return hard_edges, soft_edges, edge_rows, {"weak_phaselet_edges": weak, "conflicting_phaselet_edges": conflicting}


__all__ = [
    "PhaseEdge",
    "classify_parity_edges",
    "classify_phase_edge",
    "drop_status",
    "edge_row",
    "hard_edge_passes",
    "make_phase_edge",
    "phase_edge_risk",
    "soft_edge_passes",
    "soft_edge_rows",
]
