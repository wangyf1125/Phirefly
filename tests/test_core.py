from pathlib import Path
import math

import numpy as np
import pysam

from phirefly.hyperread import build_hyperread_matrix
from phirefly.core.matrix import load_npz_graph
from phirefly.extract import SnpSite, build_adjacent_parity_edges, write_parity_edges
from phirefly.hic.orient import orient_components_msf, retained_edges
from phirefly.metrics import compute_metrics_rows
from phirefly.phaselet import build_phaselets, load_parity_edges
from phirefly.vcf import export_phased_vcf


DATA = Path(__file__).resolve().parents[1] / "examples" / "mhc_smoke" / "data"


def test_npz_graph_loads_packaged_mhc_smoke():
    read_ids, snp_ids, matrix, read_entries, observation_count = load_npz_graph(
        str(DATA / "phirefly_reads.observations.npz"),
        weighted=False,
    )
    assert len(read_ids) == 12259
    assert len(snp_ids) == 10200
    assert observation_count == 498889
    assert matrix.shape == (12259, 10200)
    assert matrix.nnz > 0
    assert len(read_entries) == len(read_ids)


def test_phaselet_build_toy_hard_and_soft_links():
    snp_ids = ["chr1:10:A>C", "chr1:20:A>C", "chr1:30:A>C", "chr1:40:A>C"]
    tau = np.ones(4, dtype=np.int8)
    read_entries = [
        [(0, 1, 1.0), (1, 1, 1.0)],
        [(0, 1, 1.0), (1, 1, 1.0)],
        [(1, 1, 1.0), (2, 1, 1.0)],
        [(2, 1, 1.0), (3, 1, 1.0)],
    ]
    phaselet_ids, phaselets, stats, edge_rows, soft_edges = build_phaselets(
        read_entries,
        snp_ids=snp_ids,
        tau_vec=tau,
        min_abs_support=2.0,
        min_count=2,
        min_confidence=0.9,
        max_distance_bp=100,
        soft_bridge_mode="phaselet",
        soft_min_abs_support=1.0,
        soft_min_count=1,
        soft_min_confidence=0.5,
        soft_risk_threshold=1.0,
    )
    assert int(stats["retained_phaselet_edges"]) == 1
    assert int(stats["soft_bridge_edges"]) >= 1
    assert phaselet_ids[0] == phaselet_ids[1]
    assert phaselet_ids[2] != phaselet_ids[3]
    assert len(phaselets) == len(set(map(int, phaselet_ids)))
    assert edge_rows
    assert soft_edges


def test_hyperread_profiles_single_phaselet_reads():
    read_entries = [
        [(0, 1, 2.0), (1, 1, 1.0)],
        [(2, 1, 1.0), (3, 1, 1.0)],
        [(1, 1, 1.0), (2, 1, 1.0)],
        [(1, -1, 1.0), (2, -1, 1.0)],
    ]
    tau = np.ones(4, dtype=np.int8)
    phaselet_ids = np.asarray([0, 0, 1, 1], dtype=np.int64)
    matrix, hyperread_rows, hyperread_edge_rows, stats = build_hyperread_matrix(
        read_entries,
        tau_vec=tau,
        phaselet_ids=phaselet_ids,
        n_phaselets=2,
        hyperread_norm="sum",
    )
    assert stats["single_phaselet_reads"] == 2
    assert stats["multi_phaselet_reads"] == 2
    assert stats["hyperreads"] == 1
    assert matrix.shape == (1, 2)
    assert matrix.nnz == 2
    assert len(hyperread_rows) == 1
    assert len(hyperread_edge_rows) == 2


def test_extraction_time_parity_edge_roundtrip(tmp_path):
    sites = [
        SnpSite("chr1", 9, "chr1:10:A>C", "A", "C"),
        SnpSite("chr1", 19, "chr1:20:A>C", "A", "C"),
    ]
    rows = [
        (0, 0, 1, 30, 60, 0.001, 2.0),
        (0, 1, 1, 30, 60, 0.001, 1.0),
        (1, 0, -1, 30, 60, 0.001, 1.0),
        (1, 1, -1, 30, 60, 0.001, 1.0),
    ]
    edge_rows = build_adjacent_parity_edges(rows, sites, max_distance_bp=100)
    path = tmp_path / "edges.tsv.gz"
    write_parity_edges(path, edge_rows)
    loaded = load_parity_edges(path)
    assert len(edge_rows) == 1
    score, support, count, max_abs = loaded[(0, 1)]
    assert math.isclose(score, 2.0)
    assert math.isclose(support, 2.0)
    assert count == 2
    assert math.isclose(max_abs, 1.0)


def write_vcf(path: Path, sample_columns: list[str], rows: list[str]) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=1000>",
                "##INFO=<ID=ZZ,Number=1,Type=String,Description=\"test info\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_columns),
                *rows,
            ]
        )
        + "\n"
    )


def test_vcf_export_preserves_header_sample_and_adds_ps(tmp_path):
    input_vcf = tmp_path / "input.vcf"
    phases = tmp_path / "phases.tsv"
    output_vcf = tmp_path / "phased.vcf"
    write_vcf(
        input_vcf,
        ["S1", "S2"],
        [
            "chr1\t10\t.\tA\tC\t.\tPASS\tZZ=x\tGT\t0/1\t0/1",
            "chr1\t20\t.\tG\tT\t.\tPASS\tZZ=y\tGT\t0/1\t0/1",
        ],
    )
    phases.write_text("snp_id\tdelta\nchr1:10:A>C\t1\nchr1:20:G>T\t-1\n")

    written, skipped, _ = export_phased_vcf(
        input_vcf=input_vcf,
        pred_snp_phases=phases,
        region="chr1:1-30",
        output_vcf=output_vcf,
        sample="S2",
        ps_start=1,
    )

    assert written == 2
    assert skipped == 0
    with pysam.VariantFile(output_vcf) as vcf:
        assert "chr1" in vcf.header.contigs
        assert "ZZ" in vcf.header.info
        assert "PS" in vcf.header.formats
        assert list(vcf.header.samples) == ["S2"]
        calls = list(vcf)
    assert calls[0].samples["S2"]["GT"] == (0, 1)
    assert calls[0].samples["S2"].phased
    assert calls[0].samples["S2"]["PS"] == 1
    assert calls[1].samples["S2"]["GT"] == (1, 0)


def test_metrics_uses_requested_pred_sample(tmp_path):
    input_vcf = tmp_path / "input.vcf"
    truth_vcf = tmp_path / "truth.vcf"
    pred_vcf = tmp_path / "pred.vcf"
    write_vcf(
        input_vcf,
        ["S1"],
        [
            "chr1\t10\t.\tA\tC\t.\tPASS\tZZ=x\tGT\t0/1",
            "chr1\t20\t.\tG\tT\t.\tPASS\tZZ=y\tGT\t0/1",
        ],
    )
    truth_vcf.write_text(
        input_vcf.read_text().replace("GT\t0/1", "GT\t0|1").replace("GT\t0/1", "GT\t0|1")
    )
    write_vcf(
        pred_vcf,
        ["wrong", "right"],
        [
            "chr1\t10\t.\tA\tC\t.\tPASS\tZZ=x\tGT\t1|0\t0|1",
            "chr1\t20\t.\tG\tT\t.\tPASS\tZZ=y\tGT\t1|0\t0|1",
        ],
    )
    rows = compute_metrics_rows(
        input_vcf=input_vcf,
        truth_vcf=truth_vcf,
        region_text="chr1:1-30",
        input_sample="S1",
        truth_sample="S1",
        pred_sample="right",
        methods=[("pred", pred_vcf, [])],
    )
    assert rows[0]["HE_percent"] == 0.0
    assert rows[0]["phased_block_N50_kb"] == rows[0]["haplotype_N50_kb"]


def test_hic_component_orientation_is_msf_only():
    rows = [
        {
            "component_i": "10",
            "component_j": "20",
            "score": "5.0",
            "abs_score": "5.0",
            "contacts": "20",
            "confidence": "1.0",
        },
        {
            "component_i": "20",
            "component_j": "30",
            "score": "-4.0",
            "abs_score": "4.0",
            "contacts": "20",
            "confidence": "1.0",
        },
    ]
    edges = retained_edges(rows, min_contacts=10, min_confidence=0.8, min_abs_score=0.0)
    orientation, phase_sets, selected = orient_components_msf(edges)
    assert len(selected) == 2
    assert orientation["20"] == orientation["10"]
    assert orientation["30"] == -orientation["20"]
    assert phase_sets["10"] == "10"
    assert phase_sets["20"] == "10"
    assert phase_sets["30"] == "10"
