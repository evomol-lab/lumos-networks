"""Pure, testable analytical helpers used by Lumos Networks.

The Streamlit pages should remain presentation/orchestration layers.  Functions
in this module deliberately avoid Streamlit state so that the numerical
behaviour can be unit tested and cited precisely.
"""

from __future__ import annotations

from importlib import metadata
from itertools import product
import re
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


__version__ = "1.3.0-beta"


def sanitize_design_name(name: str) -> str:
    """Return a deterministic Formulaic/R-safe column name."""

    value = re.sub(r"[^A-Za-z0-9_]", "_", str(name).strip())
    value = re.sub(r"_+", "_", value).strip("_") or "variable"
    if value[0].isdigit():
        value = f"v_{value}"
    if value in {"condition", "Intercept"}:
        value = f"cov_{value}"
    return value


def collapse_symbol_results(
    results: pd.DataFrame,
    value_columns: Sequence[str] = ("Log2FC", "PValue", "FDR"),
) -> pd.DataFrame:
    """Collapse duplicate symbols without index/column-name ambiguity."""

    required = ["Symbol", *value_columns]
    missing = [column for column in required if column not in results.columns]
    if missing:
        raise ValueError("Result columns are missing: " + ", ".join(missing))
    frame = results.loc[:, required].copy().reset_index(drop=True)
    frame["Symbol"] = frame["Symbol"].astype("string").str.strip()
    frame = frame.dropna(subset=["Symbol"])
    frame = frame[frame["Symbol"] != ""]
    for column in value_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return (
        frame.groupby("Symbol", as_index=False, sort=False)[list(value_columns)]
        .mean()
    )


def annotate_cross_dataset_patterns(
    results: pd.DataFrame, datasets: Sequence[str]
) -> pd.DataFrame:
    """Annotate conserved, shared, discordant and dataset-specific DE patterns."""

    datasets = list(dict.fromkeys(str(dataset) for dataset in datasets))
    if len(datasets) < 2:
        raise ValueError("Cross-dataset annotation requires at least two datasets.")
    required = [f"Status_{dataset}" for dataset in datasets]
    missing = [column for column in required if column not in results.columns]
    if missing:
        raise ValueError("Cross-dataset status columns are missing: " + ", ".join(missing))

    annotated = results.copy()

    def summarize(row: pd.Series) -> pd.Series:
        significant: list[str] = []
        up: list[str] = []
        down: list[str] = []
        for dataset in datasets:
            status = row.get(f"Status_{dataset}")
            if status == "Up-regulated":
                significant.append(dataset)
                up.append(dataset)
            elif status == "Down-regulated":
                significant.append(dataset)
                down.append(dataset)

        n_significant = len(significant)
        if n_significant == 0:
            pattern = "Not significant"
        elif up and down:
            pattern = "Discordant"
        elif n_significant == len(datasets) and len(up) == len(datasets):
            pattern = "Conserved Up"
        elif n_significant == len(datasets) and len(down) == len(datasets):
            pattern = "Conserved Down"
        elif n_significant == 1:
            direction = "Up" if up else "Down"
            pattern = f"Specific {significant[0]} ({direction})"
        elif up:
            pattern = "Shared Up"
        else:
            pattern = "Shared Down"

        return pd.Series({
            "Cross_Tissue_Pattern": pattern,
            "N_Significant_Datasets": n_significant,
            "N_Up": len(up),
            "N_Down": len(down),
            "Significant_Datasets": ";".join(significant),
        })

    annotations = annotated.apply(summarize, axis=1)
    return pd.concat([annotated.reset_index(drop=True), annotations.reset_index(drop=True)], axis=1)


def align_named_contrasts(
    contrast_results: Mapping[str, pd.DataFrame],
    fdr_threshold: float,
    effect_threshold: float,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Align multiple contrasts from one dataset without treating them as studies."""

    if len(contrast_results) < 2:
        raise ValueError("At least two completed contrasts are required for alignment.")
    merged = pd.DataFrame()
    safe_to_label: dict[str, str] = {}
    for label, results in contrast_results.items():
        safe = sanitize_design_name(label)
        if safe in safe_to_label and safe_to_label[safe] != label:
            raise ValueError(
                f"Contrast labels '{safe_to_label[safe]}' and '{label}' resolve to the same name."
            )
        safe_to_label[safe] = str(label)
        collapsed = collapse_symbol_results(results)
        collapsed["Status"] = classify_differential_expression(
            collapsed, "Log2FC", "FDR", fdr_threshold, effect_threshold
        )
        collapsed = collapsed.rename(columns={
            "Log2FC": f"Log2FC_{safe}", "PValue": f"PValue_{safe}",
            "FDR": f"FDR_{safe}", "Status": f"Status_{safe}",
        })
        merged = (
            collapsed if merged.empty
            else pd.merge(merged, collapsed, on="Symbol", how="outer")
        )

    safe_names = list(safe_to_label)
    merged = annotate_cross_dataset_patterns(merged, safe_names).rename(columns={
        "Cross_Tissue_Pattern": "Cross_Contrast_Pattern",
        "N_Significant_Datasets": "N_Significant_Contrasts",
        "Significant_Datasets": "Significant_Contrasts",
    })
    merged["Cross_Contrast_Pattern"] = merged["Cross_Contrast_Pattern"].replace({
        "Conserved Up": "Consistent Up (all contrasts)",
        "Conserved Down": "Consistent Down (all contrasts)",
    })
    for safe, label in safe_to_label.items():
        merged["Cross_Contrast_Pattern"] = merged["Cross_Contrast_Pattern"].replace({
            f"Specific {safe} (Up)": f"Specific {label} (Up)",
            f"Specific {safe} (Down)": f"Specific {label} (Down)",
        })
    merged["Analysis_Scope"] = "single_dataset_multi_contrast"
    return merged, safe_to_label


def resolve_sample_columns(matrix_columns: Sequence[str], accessions: Sequence[str]) -> list[str]:
    """Map each requested GSM accession to one and only one matrix column."""

    columns = [str(column) for column in matrix_columns]
    if len(columns) != len(set(columns)):
        duplicates = sorted({column for column in columns if columns.count(column) > 1})
        raise ValueError("Expression matrix has duplicate sample columns: " + ", ".join(duplicates[:8]))
    resolved: list[str] = []
    for accession in dict.fromkeys(str(x) for x in accessions):
        exact = [column for column in columns if column == accession]
        matches = exact or [
            column for column in columns
            if re.search(rf"(?<!\d){re.escape(accession)}(?!\d)", column, re.IGNORECASE)
        ]
        if len(matches) != 1:
            raise ValueError(
                f"{accession} maps to {len(matches)} expression columns; expected exactly one."
            )
        resolved.append(matches[0])
    return resolved


def orient_gene_sample_matrix(
    values: np.ndarray, n_genes: int, sample_ids: Sequence[str]
) -> np.ndarray:
    """Return backend output in genes-by-samples orientation or fail clearly."""

    matrix = np.asarray(values)
    expected = (int(n_genes), len(sample_ids))
    if matrix.shape == expected:
        return matrix
    if matrix.shape == expected[::-1]:
        return matrix.T
    raise ValueError(
        f"Unexpected normalized matrix shape {matrix.shape}; expected {expected} "
        f"(genes×samples) or {expected[::-1]} (samples×genes)."
    )


def make_unique_sample_labels(
    metadata_frame: pd.DataFrame, label_columns: Sequence[str]
) -> pd.DataFrame:
    """Deduplicate GEO metadata by GSM and create accession-prefixed UI labels."""

    if "Accession" not in metadata_frame.columns:
        raise ValueError("GEO metadata do not contain an Accession column.")
    metadata = metadata_frame.copy()
    metadata["Accession"] = metadata["Accession"].astype(str).str.extract(
        r"(GSM\d+)", expand=False
    )
    metadata = metadata.dropna(subset=["Accession"]).drop_duplicates(
        "Accession", keep="first"
    ).reset_index(drop=True)
    labels: list[str] = []
    for _, row in metadata.iterrows():
        details = [
            str(row[column]).strip() for column in label_columns
            if column in row.index and column != "Accession" and str(row[column]).strip()
        ]
        labels.append(" | ".join([row["Accession"]] + details))
    metadata["display_label"] = labels
    if metadata["display_label"].duplicated().any():
        raise ValueError("Sample labels remain duplicated after adding GSM accessions.")
    return metadata


def align_supplementary_columns(
    supplementary: pd.DataFrame, gsm_order: Sequence[str]
) -> tuple[pd.DataFrame, str]:
    """Align heterogeneous GEO supplementary columns to unique GSM accessions."""

    frame = supplementary.copy()
    frame.columns = [str(column).strip().replace('"', '') for column in frame.columns]
    unique_gsms = list(dict.fromkeys(str(x) for x in gsm_order))
    gsm_set = set(unique_gsms)
    mapped = []
    for column in frame.columns:
        match = re.search(r"GSM\d+", column.upper())
        mapped.append(match.group(0) if match else column)
    positions_by_gsm: dict[str, int] = {}
    for position, value in enumerate(mapped):
        if value not in gsm_set:
            continue
        if value in positions_by_gsm:
            raise ValueError("Multiple supplementary columns map to the same GSM accession.")
        positions_by_gsm[value] = position
    if positions_by_gsm and set(positions_by_gsm) == gsm_set:
        synced = frame.iloc[:, [positions_by_gsm[gsm] for gsm in unique_gsms]].copy()
        synced.columns = unique_gsms
        return synced, "GSM identified in column name"
    if (
        not positions_by_gsm
        and len(frame.columns) == len(unique_gsms)
        and len(unique_gsms) == len(gsm_order)
    ):
        frame.columns = unique_gsms
        return frame, "positional mapping validated by cardinality"
    raise ValueError(
        f"Cannot align {len(frame.columns)} supplementary columns with "
        f"{len(unique_gsms)} unique GEO samples without an unambiguous GSM identifier."
    )


def assemble_geo_tar_tables(
    parsed_members: Iterable[tuple[str, pd.DataFrame | pd.Series]],
    gsm_order: Sequence[str],
) -> tuple[pd.DataFrame, str]:
    """Assemble a GEO RAW TAR containing either a matrix or one file per GSM.

    A one-column DataFrame is a single-sample table, not a consolidated matrix.
    The previous loader stopped at the first such table and then attempted to
    align one column against every sample in the GEO Series.
    """

    unique_gsms = list(dict.fromkeys(str(value) for value in gsm_order))
    gsm_set = set(unique_gsms)
    single_samples: dict[str, pd.Series] = {}
    matrix_errors: list[str] = []

    for filename, parsed in parsed_members:
        if parsed is None:
            continue

        if isinstance(parsed, pd.DataFrame) and parsed.shape[1] >= 2:
            try:
                aligned, method = align_supplementary_columns(parsed, unique_gsms)
                return aligned, f"consolidated matrix; {method}"
            except ValueError as exc:
                matrix_errors.append(str(exc))
                continue

        if isinstance(parsed, pd.DataFrame):
            if parsed.shape[1] != 1:
                continue
            column_label = str(parsed.columns[0])
            sample_values = parsed.iloc[:, 0].copy()
        elif isinstance(parsed, pd.Series):
            column_label = str(parsed.name or "")
            sample_values = parsed.copy()
        else:
            continue

        accession = None
        for candidate in (filename, column_label):
            match = re.search(r"GSM\d+", str(candidate).upper())
            if match and match.group(0) in gsm_set:
                accession = match.group(0)
                break
        if accession is None:
            continue
        if accession in single_samples:
            raise ValueError(
                f"Multiple one-sample TAR members map to {accession}; "
                "the count source is ambiguous."
            )
        single_samples[accession] = sample_values.rename(accession)

    if single_samples:
        ordered = [single_samples[gsm] for gsm in unique_gsms if gsm in single_samples]
        combined = pd.concat(ordered, axis=1, join="outer")
        aligned, method = align_supplementary_columns(combined, unique_gsms)
        return aligned, f"per-sample TAR tables; {method}"

    detail = f" Last matrix error: {matrix_errors[-1]}" if matrix_errors else ""
    raise ValueError(f"No usable GEO count tables were found in the TAR.{detail}")


def build_design_metadata(
    sample_ids: Sequence[str],
    conditions: Sequence[str],
    source_metadata: pd.DataFrame | None = None,
    sample_key: str = "Accession",
    categorical_covariates: Sequence[str] = (),
    continuous_covariates: Sequence[str] = (),
    batch_column: str | None = None,
    pair_column: str | None = None,
) -> tuple[pd.DataFrame, str, list[str], dict[str, str]]:
    """Build and validate a full-rank differential-expression design.

    The returned frame is indexed in the exact expression-matrix order. Batch,
    categorical covariates and subject/pair identifiers are encoded as strings;
    continuous covariates are numeric. Paired designs require each retained
    subject to contribute samples to both conditions.
    """

    sample_ids = [str(x) for x in sample_ids]
    conditions = [str(x) for x in conditions]
    if len(sample_ids) != len(conditions) or not sample_ids:
        raise ValueError("Sample identifiers and conditions must have equal non-zero length.")
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Sample identifiers must be unique in the design.")
    levels = list(dict.fromkeys(conditions))
    if len(levels) != 2:
        raise ValueError("Differential expression currently requires exactly two conditions.")

    requested = list(dict.fromkeys(
        [x for x in [batch_column, pair_column] if x]
        + list(categorical_covariates)
        + list(continuous_covariates)
    ))
    design = pd.DataFrame(index=pd.Index(sample_ids, name="sample"))
    design["condition"] = conditions

    if requested:
        if source_metadata is None or source_metadata.empty:
            raise ValueError("Covariates were requested but sample metadata are unavailable.")
        if sample_key not in source_metadata.columns:
            raise ValueError(f"Metadata sample key '{sample_key}' is missing.")
        indexed = source_metadata.copy()
        indexed[sample_key] = indexed[sample_key].astype(str)
        if indexed[sample_key].duplicated().any():
            raise ValueError(f"Metadata sample key '{sample_key}' contains duplicates.")
        indexed = indexed.set_index(sample_key)
        missing_samples = sorted(set(sample_ids).difference(indexed.index))
        if missing_samples:
            raise ValueError("Metadata are missing selected samples: " + ", ".join(missing_samples[:8]))
        missing_columns = [x for x in requested if x not in indexed.columns]
        if missing_columns:
            raise ValueError("Metadata columns are missing: " + ", ".join(missing_columns))
    else:
        indexed = pd.DataFrame(index=sample_ids)

    name_map: dict[str, str] = {}
    used_names = {"condition"}
    terms: list[str] = []

    def unique_name(original: str) -> str:
        candidate = sanitize_design_name(original)
        base, suffix = candidate, 2
        while candidate in used_names:
            candidate = f"{base}_{suffix}"
            suffix += 1
        used_names.add(candidate)
        name_map[original] = candidate
        return candidate

    categorical = list(dict.fromkeys(
        [x for x in [batch_column, pair_column] if x] + list(categorical_covariates)
    ))
    for column in categorical:
        safe = unique_name(column)
        values = indexed.loc[sample_ids, column]
        if values.isna().any() or values.astype(str).str.strip().isin(["", "nan", "None"]).any():
            raise ValueError(f"Categorical variable '{column}' contains missing values.")
        design[safe] = values.astype(str).to_numpy()
        if design[safe].nunique() < 2:
            raise ValueError(f"Categorical variable '{column}' has fewer than two levels.")
        terms.append(safe)

    for column in dict.fromkeys(continuous_covariates):
        if column in categorical:
            raise ValueError(f"Variable '{column}' cannot be both categorical and continuous.")
        safe = unique_name(column)
        values = pd.to_numeric(indexed.loc[sample_ids, column], errors="coerce")
        if values.isna().any() or not np.isfinite(values.to_numpy()).all():
            raise ValueError(f"Continuous variable '{column}' contains missing or non-numeric values.")
        if values.nunique() < 2:
            raise ValueError(f"Continuous variable '{column}' has no variation.")
        design[safe] = values.to_numpy(dtype=float)
        terms.append(safe)

    if pair_column:
        pair_name = name_map[pair_column]
        pairing = design.groupby(pair_name)["condition"].nunique()
        invalid = pairing[pairing < 2].index.astype(str).tolist()
        if invalid:
            raise ValueError(
                "Paired design requires both conditions for every subject; incomplete subjects: "
                + ", ".join(invalid[:8])
            )

    # A simple treatment-coded model-rank diagnostic catches confounded batch,
    # subject and condition variables before a backend produces opaque errors.
    encoded = pd.get_dummies(design, drop_first=True, dtype=float)
    matrix = np.column_stack([np.ones(len(encoded)), encoded.to_numpy(dtype=float)])
    if np.linalg.matrix_rank(matrix) < matrix.shape[1]:
        raise ValueError(
            "The requested design is not full rank. A covariate is redundant or confounded "
            "with condition; revise batch/covariate/pair selections."
        )

    formula = "~" + " + ".join(terms + ["condition"])
    contrast = ["condition", levels[1], levels[0]]
    return design, formula, contrast, name_map


def prepare_ora_inputs(
    query_genes: Iterable[str], background_genes: Iterable[str]
) -> tuple[list[str], list[str], list[str]]:
    """Normalize an ORA query/universe and exclude query genes outside it."""

    def normalize(values: Iterable[str]) -> list[str]:
        return sorted({str(x).strip().upper() for x in values if str(x).strip()})

    query = normalize(query_genes)
    background = normalize(background_genes)
    if not query:
        raise ValueError("ORA query gene list is empty.")
    if not background:
        raise ValueError("ORA background/universe is empty.")
    background_set = set(background)
    excluded = sorted(set(query).difference(background_set))
    retained = sorted(set(query).intersection(background_set))
    if not retained:
        raise ValueError("No query genes occur in the selected ORA background.")
    return retained, background, excluded


def benchmark_concordance(
    left: pd.DataFrame,
    right: pd.DataFrame,
    gene_column: str = "Symbol",
    effect_column: str = "Log2FC",
    fdr_column: str = "FDR",
    alpha: float = 0.05,
    top_k: int = 100,
) -> dict[str, float | int]:
    """Compute transparent concordance metrics for two DE result tables."""

    required = {gene_column, effect_column, fdr_column}
    for label, frame in (("left", left), ("right", right)):
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{label} benchmark table is missing: {', '.join(sorted(missing))}")
    merged = left[list(required)].merge(
        right[list(required)], on=gene_column, suffixes=("_left", "_right")
    ).dropna()
    if merged.empty:
        raise ValueError("Benchmark tables have no complete genes in common.")
    lfc_l = pd.to_numeric(merged[f"{effect_column}_left"], errors="coerce")
    lfc_r = pd.to_numeric(merged[f"{effect_column}_right"], errors="coerce")
    fdr_l = pd.to_numeric(merged[f"{fdr_column}_left"], errors="coerce")
    fdr_r = pd.to_numeric(merged[f"{fdr_column}_right"], errors="coerce")
    valid = lfc_l.notna() & lfc_r.notna() & fdr_l.notna() & fdr_r.notna()
    merged, lfc_l, lfc_r, fdr_l, fdr_r = (
        merged.loc[valid], lfc_l.loc[valid], lfc_r.loc[valid], fdr_l.loc[valid], fdr_r.loc[valid]
    )
    sig_l, sig_r = set(merged.loc[fdr_l < alpha, gene_column]), set(merged.loc[fdr_r < alpha, gene_column])
    union_sig = sig_l | sig_r
    top_l = set(merged.assign(_f=fdr_l).nsmallest(min(top_k, len(merged)), "_f")[gene_column])
    top_r = set(merged.assign(_f=fdr_r).nsmallest(min(top_k, len(merged)), "_f")[gene_column])
    direction = np.sign(lfc_l) == np.sign(lfc_r)
    return {
        "n_common": int(len(merged)),
        "pearson_log2fc": float(lfc_l.corr(lfc_r, method="pearson")),
        "spearman_log2fc": float(lfc_l.corr(lfc_r, method="spearman")),
        "rmse_log2fc": float(np.sqrt(np.mean(np.square(lfc_l - lfc_r)))),
        "direction_agreement": float(direction.mean()),
        "significant_jaccard": float(len(sig_l & sig_r) / len(union_sig)) if union_sig else 1.0,
        "top_k_jaccard": float(len(top_l & top_r) / len(top_l | top_r)),
    }


def quantile_normalize(values: np.ndarray) -> np.ndarray:
    """Quantile-normalize a genes-by-samples matrix, averaging tied ranks.

    Parameters
    ----------
    values
        Two-dimensional numeric matrix without missing or infinite values.

    Returns
    -------
    numpy.ndarray
        Float64 matrix with the same shape. Equal input values within a sample
        receive equal normalized values.
    """

    matrix = np.asarray(values, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("Quantile normalization requires a two-dimensional matrix.")
    if matrix.size == 0:
        return matrix.copy()
    if not np.isfinite(matrix).all():
        raise ValueError("Quantile normalization does not accept missing or infinite values.")

    rank_means = np.sort(matrix, axis=0).mean(axis=1)
    normalized = np.empty_like(matrix, dtype=float)

    for column_index in range(matrix.shape[1]):
        column = matrix[:, column_index]
        order = np.argsort(column, kind="mergesort")
        sorted_column = column[order]
        start = 0
        while start < len(order):
            stop = start + 1
            while stop < len(order) and sorted_column[stop] == sorted_column[start]:
                stop += 1
            normalized[order[start:stop], column_index] = rank_means[start:stop].mean()
            start = stop

    return normalized


def validate_count_matrix(counts: pd.DataFrame) -> pd.DataFrame:
    """Validate and return an integer RNA-seq count matrix.

    PyDESeq2 requires finite, non-negative integer counts.  Silent rounding or
    truncation is intentionally rejected because it can turn normalized
    expression values into invalid count data.
    """

    numeric = counts.apply(pd.to_numeric, errors="coerce")
    if numeric.empty:
        raise ValueError("The RNA-seq count matrix is empty.")
    if numeric.isna().any().any() or not np.isfinite(numeric.to_numpy()).all():
        raise ValueError("RNA-seq counts contain missing, non-numeric, or infinite values.")
    if (numeric.to_numpy() < 0).any():
        raise ValueError("RNA-seq counts must be non-negative.")
    if not np.allclose(numeric.to_numpy(), np.rint(numeric.to_numpy()), rtol=0, atol=1e-8):
        raise ValueError(
            "PyDESeq2 requires raw integer counts; transformed or normalized values were detected."
        )
    max_count = numeric.to_numpy().max(initial=0)
    if max_count > np.iinfo(np.int64).max:
        raise ValueError("RNA-seq counts exceed the supported 64-bit integer range.")
    return numeric.astype(np.int64)


def classify_differential_expression(
    frame: pd.DataFrame,
    effect_column: str,
    fdr_column: str,
    fdr_threshold: float,
    absolute_effect_threshold: float,
) -> pd.Series:
    """Classify rows as Up-regulated, Down-regulated, or Not Significant."""

    effect = pd.to_numeric(frame[effect_column], errors="coerce")
    fdr = pd.to_numeric(frame[fdr_column], errors="coerce")
    status = pd.Series("Not Significant", index=frame.index, dtype="object")
    significant = fdr.lt(fdr_threshold)
    status.loc[significant & effect.ge(absolute_effect_threshold)] = "Up-regulated"
    status.loc[significant & effect.le(-absolute_effect_threshold)] = "Down-regulated"
    return status


def conserved_signature(
    frame: pd.DataFrame,
    dataset_ids: Sequence[str],
    fdr_threshold: float,
    absolute_effect_threshold: float,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Select genes significant in every dataset with a concordant direction.

    This is intentionally stricter than averaging effects or taking the minimum
    FDR, either of which can mislabel discordant evidence as conserved.
    """

    if len(dataset_ids) < 2:
        raise ValueError("A conserved signature requires at least two datasets.")
    required = ["Symbol"]
    for dataset_id in dataset_ids:
        required.extend([f"Log2FC_{dataset_id}", f"FDR_{dataset_id}"])
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError("Missing multi-dataset columns: " + ", ".join(missing))

    result = frame.copy()
    effects = result[[f"Log2FC_{x}" for x in dataset_ids]].apply(pd.to_numeric, errors="coerce")
    fdrs = result[[f"FDR_{x}" for x in dataset_ids]].apply(pd.to_numeric, errors="coerce")
    significant_everywhere = fdrs.lt(fdr_threshold).all(axis=1)
    large_everywhere = effects.abs().ge(absolute_effect_threshold).all(axis=1)
    all_up = effects.gt(0).all(axis=1)
    all_down = effects.lt(0).all(axis=1)

    result["Effect"] = effects.mean(axis=1)
    result["FDR"] = fdrs.max(axis=1)
    result["Status"] = "Not conserved"
    result.loc[significant_everywhere & large_everywhere & all_up, "Status"] = "Conserved up"
    result.loc[significant_everywhere & large_everywhere & all_down, "Status"] = "Conserved down"
    conserved = result[result["Status"].ne("Not conserved")].copy()

    genes = {
        "Up-regulated": sorted(
            conserved.loc[conserved["Status"].eq("Conserved up"), "Symbol"]
            .dropna().astype(str).str.strip().str.upper().unique().tolist()
        ),
        "Down-regulated": sorted(
            conserved.loc[conserved["Status"].eq("Conserved down"), "Symbol"]
            .dropna().astype(str).str.strip().str.upper().unique().tolist()
        ),
    }
    return conserved, genes


def random_walk_with_restart(
    adjacency: np.ndarray,
    seeds: Iterable[float],
    restart_probability: float = 0.7,
    tolerance: float = 1e-10,
    max_iterations: int = 1000,
) -> tuple[np.ndarray, int, bool]:
    """Run a personalized random walk on a column-oriented adjacency matrix.

    Dangling columns redistribute probability according to the seed vector, so
    probability mass is conserved at every iteration.
    """

    matrix = np.asarray(adjacency, dtype=float)
    seed = np.asarray(list(seeds), dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Adjacency matrix must be square.")
    if seed.shape != (matrix.shape[0],):
        raise ValueError("Seed vector length must match the adjacency matrix.")
    if not 0 < restart_probability < 1:
        raise ValueError("Restart probability must be between 0 and 1.")
    if (matrix < 0).any() or not np.isfinite(matrix).all():
        raise ValueError("Adjacency weights must be finite and non-negative.")
    if (seed < 0).any() or not np.isfinite(seed).all():
        raise ValueError("Seed weights must be finite and non-negative.")

    if seed.sum() == 0:
        seed = np.full(matrix.shape[0], 1.0 / matrix.shape[0])
    else:
        seed = seed / seed.sum()

    transition = np.zeros_like(matrix, dtype=float)
    column_sums = matrix.sum(axis=0)
    non_dangling = column_sums > 0
    transition[:, non_dangling] = matrix[:, non_dangling] / column_sums[non_dangling]
    transition[:, ~non_dangling] = seed[:, None]

    probability = seed.copy()
    for iteration in range(1, max_iterations + 1):
        updated = (
            (1.0 - restart_probability) * transition.dot(probability)
            + restart_probability * seed
        )
        if np.linalg.norm(updated - probability, ord=1) < tolerance:
            return updated / updated.sum(), iteration, True
        probability = updated
    return probability / probability.sum(), max_iterations, False


def network_sensitivity_analysis(
    edges: pd.DataFrame,
    evidence: pd.DataFrame,
    fdr_thresholds: Sequence[float] = (0.01, 0.05, 0.10),
    effect_thresholds: Sequence[float] = (0.5, 1.0, 1.5),
    restart_probabilities: Sequence[float] = (0.5, 0.7, 0.9),
    propagation_decays: Sequence[float] = (0.25, 0.5, 0.75),
    top_k: int = 20,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Perturb network-selection and propagation parameters.

    Returns a long per-configuration ranking table and a per-node stability
    summary. This is a parameter sensitivity analysis, not a bootstrap or an
    estimate of biological uncertainty.
    """

    edge_required = {"TF", "Target"}
    evidence_required = {"Symbol", "Effect", "FDR"}
    if missing := edge_required.difference(edges.columns):
        raise ValueError("Network edges are missing: " + ", ".join(sorted(missing)))
    if missing := evidence_required.difference(evidence.columns):
        raise ValueError("Network evidence is missing: " + ", ".join(sorted(missing)))
    if top_k < 1:
        raise ValueError("top_k must be positive.")

    clean_edges = edges[["TF", "Target"]].dropna().astype(str).drop_duplicates()
    clean_edges = clean_edges[(clean_edges["TF"] != "") & (clean_edges["Target"] != "")]
    ev = evidence[["Symbol", "Effect", "FDR"]].copy()
    ev["Symbol"] = ev["Symbol"].astype(str)
    ev["Effect"] = pd.to_numeric(ev["Effect"], errors="coerce")
    ev["FDR"] = pd.to_numeric(ev["FDR"], errors="coerce")
    ev = ev.dropna().drop_duplicates("Symbol")
    rows: list[dict[str, float | int | str | bool]] = []

    config_id = 0
    for fdr, effect, restart, decay in product(
        fdr_thresholds, effect_thresholds, restart_probabilities, propagation_decays
    ):
        if not 0 < float(fdr) <= 1 or float(effect) < 0:
            raise ValueError("FDR must be in (0, 1] and effect thresholds non-negative.")
        if not 0 < float(restart) < 1 or not 0 <= float(decay) < 1:
            raise ValueError("Restart must be in (0, 1) and propagation decay in [0, 1).")
        selected = ev[(ev["FDR"] <= fdr) & (ev["Effect"].abs() >= effect)]
        selected_genes = set(selected["Symbol"])
        relevant = clean_edges[clean_edges["Target"].isin(selected_genes)]
        if relevant.empty:
            continue
        nodes = sorted(set(relevant["TF"]) | set(relevant["Target"]))
        index = {node: idx for idx, node in enumerate(nodes)}
        adjacency = np.zeros((len(nodes), len(nodes)), dtype=float)
        for row in relevant.itertuples(index=False):
            adjacency[index[row.Target], index[row.TF]] = 1.0
        effect_map = selected.set_index("Symbol")["Effect"].to_dict()
        signed_seed = np.array([effect_map.get(node, 0.0) for node in nodes], dtype=float)
        rwr, iterations, converged = random_walk_with_restart(
            adjacency, np.abs(signed_seed), restart_probability=float(restart)
        )
        column_sums = adjacency.sum(axis=0)
        transition = np.zeros_like(adjacency)
        non_dangling = column_sums > 0
        transition[:, non_dangling] = adjacency[:, non_dangling] / column_sums[non_dangling]
        try:
            propagated = np.linalg.solve(
                np.eye(len(nodes)) - float(decay) * transition, signed_seed
            )
        except np.linalg.LinAlgError:
            propagated = np.full(len(nodes), np.nan)
        config_id += 1
        ordering = np.argsort(-rwr, kind="stable")
        ranks = np.empty(len(nodes), dtype=int)
        ranks[ordering] = np.arange(1, len(nodes) + 1)
        for idx, node in enumerate(nodes):
            rows.append({
                "Configuration": config_id,
                "FDR_Threshold": float(fdr),
                "Effect_Threshold": float(effect),
                "Restart_Probability": float(restart),
                "Propagation_Decay": float(decay),
                "Node": node,
                "RWR_Score": float(rwr[idx]),
                "RWR_Rank": int(ranks[idx]),
                "Linear_Propagation": float(propagated[idx]),
                "TopK": bool(ranks[idx] <= min(top_k, len(nodes))),
                "RWR_Converged": bool(converged),
                "RWR_Iterations": int(iterations),
            })

    rankings = pd.DataFrame(rows)
    if rankings.empty:
        return rankings, pd.DataFrame()
    total_configs = rankings["Configuration"].nunique()
    stability = rankings.groupby("Node", as_index=False).agg(
        Configurations_Present=("Configuration", "nunique"),
        TopK_Frequency=("TopK", "mean"),
        Median_RWR_Rank=("RWR_Rank", "median"),
        IQR_RWR_Rank=("RWR_Rank", lambda x: float(x.quantile(0.75) - x.quantile(0.25))),
        Mean_RWR_Score=("RWR_Score", "mean"),
        SD_RWR_Score=("RWR_Score", "std"),
        Median_Linear_Propagation=("Linear_Propagation", "median"),
    )
    stability["Configuration_Coverage"] = stability["Configurations_Present"] / total_configs
    stability = stability.sort_values(
        ["TopK_Frequency", "Configuration_Coverage", "Median_RWR_Rank"],
        ascending=[False, False, True],
    )
    return rankings, stability


def package_versions(packages: Iterable[str]) -> Mapping[str, str]:
    """Return installed package versions without importing the packages."""

    versions: dict[str, str] = {}
    for package in packages:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = "not installed"
    return versions
