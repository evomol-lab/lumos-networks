import unittest

import numpy as np
import pandas as pd

from lumos_core import (
    align_supplementary_columns,
    align_named_contrasts,
    annotate_cross_dataset_patterns,
    assemble_geo_tar_tables,
    benchmark_concordance,
    build_design_metadata,
    classify_differential_expression,
    collapse_symbol_results,
    conserved_signature,
    make_unique_sample_labels,
    network_sensitivity_analysis,
    orient_gene_sample_matrix,
    prepare_ora_inputs,
    quantile_normalize,
    random_walk_with_restart,
    resolve_sample_columns,
    validate_count_matrix,
)


class LumosCoreTests(unittest.TestCase):
    def test_quantile_normalize_preserves_ties_and_equal_distributions(self):
        matrix = np.array([[1, 4], [1, 3], [4, 2], [5, 1]], dtype=float)
        normalized = quantile_normalize(matrix)
        self.assertEqual(normalized[0, 0], normalized[1, 0])
        # The two tied observations receive the mean of the first two target
        # quantiles: mean([1.0, 1.5]) = 1.25.
        self.assertAlmostEqual(normalized[0, 0], 1.25)
        np.testing.assert_allclose(np.sort(normalized[:, 1]), [1.0, 1.5, 3.5, 4.5])

    def test_validate_count_matrix_rejects_fractional_values(self):
        with self.assertRaisesRegex(ValueError, "raw integer counts"):
            validate_count_matrix(pd.DataFrame({"S1": [1, 2.5], "S2": [3, 4]}))

    def test_validate_count_matrix_accepts_nonnegative_integers(self):
        result = validate_count_matrix(pd.DataFrame({"S1": [0, 2], "S2": [3, 4]}))
        self.assertTrue(result.dtypes.eq(np.dtype("int64")).all())

    def test_de_classification_uses_fdr_and_effect_thresholds(self):
        frame = pd.DataFrame({"effect": [2, -2, 0.5, 3], "fdr": [0.01, 0.01, 0.01, 0.2]})
        status = classify_differential_expression(frame, "effect", "fdr", 0.05, 1.0)
        self.assertEqual(
            status.tolist(),
            ["Up-regulated", "Down-regulated", "Not Significant", "Not Significant"],
        )

    def test_conserved_signature_excludes_discordant_genes(self):
        frame = pd.DataFrame({
            "Symbol": ["A", "B", "C"],
            "Log2FC_GSE1": [2.0, 2.0, -2.0],
            "FDR_GSE1": [0.01, 0.01, 0.01],
            "Log2FC_GSE2": [1.5, -1.5, -1.2],
            "FDR_GSE2": [0.03, 0.01, 0.20],
        })
        conserved, genes = conserved_signature(frame, ["GSE1", "GSE2"], 0.05, 1.0)
        self.assertEqual(conserved["Symbol"].tolist(), ["A"])
        self.assertAlmostEqual(conserved.iloc[0]["FDR"], 0.03)
        self.assertEqual(genes, {"Up-regulated": ["A"], "Down-regulated": []})

    def test_rwr_conserves_probability_with_dangling_nodes(self):
        adjacency = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
        score, iterations, converged = random_walk_with_restart(adjacency, [1, 0, 0])
        self.assertTrue(converged)
        self.assertGreater(iterations, 0)
        self.assertAlmostEqual(score.sum(), 1.0)
        self.assertTrue((score >= 0).all())

    def test_paired_design_builds_subject_block_and_condition_contrast(self):
        metadata = pd.DataFrame({
            "Accession": ["S1", "S2", "S3", "S4"],
            "participant id": ["P1", "P1", "P2", "P2"],
            "age": [40, 40, 55, 55],
        })
        design, formula, contrast, names = build_design_metadata(
            ["S1", "S2", "S3", "S4"], ["control", "case", "control", "case"],
            metadata, pair_column="participant id",
        )
        self.assertEqual(contrast, ["condition", "case", "control"])
        self.assertIn(names["participant id"], formula)
        self.assertEqual(design.loc["S2", "condition"], "case")

    def test_design_rejects_batch_confounded_with_condition(self):
        metadata = pd.DataFrame({
            "Accession": ["S1", "S2", "S3", "S4"],
            "batch": ["B1", "B1", "B2", "B2"],
        })
        with self.assertRaisesRegex(ValueError, "not full rank"):
            build_design_metadata(
                ["S1", "S2", "S3", "S4"], ["control", "control", "case", "case"],
                metadata, batch_column="batch",
            )

    def test_ora_background_excludes_query_genes_outside_universe(self):
        query, background, excluded = prepare_ora_inputs(
            ["tp53", "IL6", "NOT_MEASURED"], ["TP53", "IL6", "STAT3"]
        )
        self.assertEqual(query, ["IL6", "TP53"])
        self.assertEqual(excluded, ["NOT_MEASURED"])
        self.assertEqual(len(background), 3)

    def test_benchmark_concordance_reports_expected_identity(self):
        frame = pd.DataFrame({
            "Symbol": ["A", "B", "C"], "Log2FC": [2.0, -1.0, 0.2],
            "FDR": [0.01, 0.03, 0.8],
        })
        metrics = benchmark_concordance(frame, frame.copy(), top_k=2)
        self.assertAlmostEqual(metrics["pearson_log2fc"], 1.0)
        self.assertAlmostEqual(metrics["significant_jaccard"], 1.0)
        self.assertAlmostEqual(metrics["top_k_jaccard"], 1.0)

    def test_network_sensitivity_returns_rank_stability(self):
        edges = pd.DataFrame({"TF": ["TF1", "TF1", "TF2"], "Target": ["A", "B", "B"]})
        evidence = pd.DataFrame({
            "Symbol": ["A", "B"], "Effect": [2.0, -1.2], "FDR": [0.01, 0.04]
        })
        rankings, stability = network_sensitivity_analysis(
            edges, evidence, fdr_thresholds=[0.05], effect_thresholds=[1.0],
            restart_probabilities=[0.5, 0.7], propagation_decays=[0.5], top_k=2,
        )
        self.assertEqual(rankings["Configuration"].nunique(), 2)
        self.assertFalse(stability.empty)
        self.assertTrue(stability["Configuration_Coverage"].between(0, 1).all())

    def test_geo_sample_alignment_is_unique_and_accession_aware(self):
        columns = ["counts_GSM5630028", "counts_GSM5630036", "unrelated"]
        self.assertEqual(
            resolve_sample_columns(columns, ["GSM5630036", "GSM5630028"]),
            ["counts_GSM5630036", "counts_GSM5630028"],
        )
        with self.assertRaisesRegex(ValueError, "duplicate sample columns"):
            resolve_sample_columns(["GSM1", "GSM1"], ["GSM1"])

    def test_geo_sample_alignment_rejects_ambiguous_accession(self):
        with self.assertRaisesRegex(ValueError, "maps to 2"):
            resolve_sample_columns(
                ["rep1_GSM5630036", "rep2_GSM5630036"], ["GSM5630036"]
            )

    def test_vst_orientation_prevents_gene_sample_shape_mismatch(self):
        samples_by_genes = np.zeros((3, 5))
        oriented = orient_gene_sample_matrix(samples_by_genes, 5, ["S1", "S2", "S3"])
        self.assertEqual(oriented.shape, (5, 3))
        with self.assertRaisesRegex(ValueError, "Unexpected normalized matrix shape"):
            orient_gene_sample_matrix(np.zeros((7, 5)), 5, ["S1", "S2", "S3"])

    def test_repeated_geo_group_labels_do_not_collapse_to_last_gsm(self):
        metadata = pd.DataFrame({
            "Accession": ["GSM5630028", "GSM5630036", "GSM5630036"],
            "source_name_ch1": ["normal skin", "lesion", "lesion"],
        })
        labelled = make_unique_sample_labels(metadata, ["source_name_ch1"])
        self.assertEqual(labelled["Accession"].tolist(), ["GSM5630028", "GSM5630036"])
        self.assertEqual(
            labelled["display_label"].tolist(),
            ["GSM5630028 | normal skin", "GSM5630036 | lesion"],
        )

    def test_geo_count_title_columns_align_positionally_only_with_equal_cardinality(self):
        counts = pd.DataFrame({"donor-A": [1, 2], "donor-B": [3, 4]})
        aligned, method = align_supplementary_columns(counts, ["GSM1", "GSM2"])
        self.assertEqual(aligned.columns.tolist(), ["GSM1", "GSM2"])
        self.assertIn("cardinality", method)
        with self.assertRaisesRegex(ValueError, "Cannot align"):
            align_supplementary_columns(counts, ["GSM1", "GSM2", "GSM3"])

    def test_tar_one_column_dataframes_are_assembled_in_gsm_order(self):
        parsed_members = [
            ("GSM3_counts.txt.gz", pd.DataFrame({"raw_counts": [30, 40]}, index=["A", "B"])),
            ("GSM1_counts.txt.gz", pd.DataFrame({"raw_counts": [10, 20]}, index=["A", "B"])),
            ("GSM2_counts.txt.gz", pd.DataFrame({"raw_counts": [50, 60]}, index=["A", "B"])),
        ]
        aligned, method = assemble_geo_tar_tables(
            parsed_members, ["GSM1", "GSM2", "GSM3"]
        )
        self.assertEqual(aligned.columns.tolist(), ["GSM1", "GSM2", "GSM3"])
        self.assertEqual(aligned.loc["A"].tolist(), [10, 50, 30])
        self.assertIn("per-sample", method)

    def test_partial_direct_gsm_mapping_is_rejected(self):
        counts = pd.DataFrame({"GSM1": [1], "GSM2": [2], "other": [3]})
        with self.assertRaisesRegex(ValueError, "Cannot align"):
            align_supplementary_columns(counts, ["GSM1", "GSM2", "GSM3"])

    def test_symbol_collapse_ignores_an_equally_named_index(self):
        results = pd.DataFrame({
            "Symbol": ["A", "A", "B"],
            "Log2FC": [1.0, 3.0, -2.0],
            "PValue": [0.01, 0.03, 0.2],
            "FDR": [0.02, 0.04, 0.3],
        })
        results.index = pd.Index(["x", "y", "z"], name="Symbol")
        collapsed = collapse_symbol_results(results)
        self.assertEqual(collapsed["Symbol"].tolist(), ["A", "B"])
        self.assertAlmostEqual(collapsed.loc[0, "Log2FC"], 2.0)

    def test_cross_dataset_patterns_support_more_than_two_tissues(self):
        frame = pd.DataFrame({
            "Symbol": ["A", "B", "C", "D"],
            "Status_GSE1": ["Up-regulated", "Up-regulated", "Up-regulated", "Not Significant"],
            "Status_GSE2": ["Up-regulated", "Down-regulated", "Not Significant", "Not Significant"],
            "Status_GSE3": ["Up-regulated", "Not Significant", "Not Significant", "Down-regulated"],
        })
        annotated = annotate_cross_dataset_patterns(frame, ["GSE1", "GSE2", "GSE3"])
        self.assertEqual(
            annotated["Cross_Tissue_Pattern"].tolist(),
            ["Conserved Up", "Discordant", "Specific GSE1 (Up)", "Specific GSE3 (Down)"],
        )
        self.assertEqual(annotated["N_Significant_Datasets"].tolist(), [3, 2, 1, 1])

    def test_multiple_contrasts_from_one_dataset_align_without_meta_analysis(self):
        first = pd.DataFrame({
            "Symbol": ["A", "B"], "Log2FC": [2.0, -2.0],
            "PValue": [0.001, 0.002], "FDR": [0.01, 0.02],
        })
        second = pd.DataFrame({
            "Symbol": ["A", "C"], "Log2FC": [1.5, 2.5],
            "PValue": [0.003, 0.004], "FDR": [0.03, 0.04],
        })
        aligned, names = align_named_contrasts(
            {"Lesion vs Control": first, "Non-lesion vs Control": second}, 0.05, 1.0
        )
        self.assertEqual(len(names), 2)
        self.assertEqual(aligned.loc[aligned["Symbol"] == "A", "N_Significant_Contrasts"].iloc[0], 2)
        self.assertEqual(
            aligned.loc[aligned["Symbol"] == "A", "Cross_Contrast_Pattern"].iloc[0],
            "Consistent Up (all contrasts)",
        )
        self.assertTrue((aligned["Analysis_Scope"] == "single_dataset_multi_contrast").all())


if __name__ == "__main__":
    unittest.main()
