#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qPCR pipeline for Bio-Rad CFX files: per-well amplification efficiency + Cq merge + efficiency-corrected relative expression.

What it does:
  1) Reads Bio-Rad amplification-curve Excel export, including Bio-Rad/CFX xlsx files
     that some Python Excel readers fail to open because of non-standard OOXML names.
  2) For each well, automatically baseline-corrects fluorescence.
  3) Automatically searches for the log-linear amplification window.
  4) Fits: log10(F_corrected) = slope * cycle + intercept
  5) Calculates:
       E_amp_factor = 10^slope
       Efficiency_% = (E_amp_factor - 1) * 100
  6) Writes one clean per-well result CSV by default.
  7) Internally self-checks the selected window; optional --debug writes diagnostic curve/window files.

This is an independent LinRegPCR/RDML-Tools-inspired implementation for Bio-Rad Excel curves, now with iterative baseline optimization by log-window slope balance.
It is intended for local analysis when RDML-Tools cannot read your RDML file.

Usage examples:
  python qpcr_linreg_efficiency_from_biorad_excel.py --amp "Quantification Amplification Results.xlsx"
  python qpcr_linreg_efficiency_from_biorad_excel.py --amp "file.xlsx" --sheet SYBR --out my_results

Output:
  my_results_per_well_efficiency.csv  (now annotated with Target/Sample/Condition from Cq file)
  my_results_cq_merged_with_efficiency.csv
  my_results_relative_expression_efficiency_corrected.csv  (if --reference is given)
  my_results_summary_by_condition.csv  (if --reference is given)
  my_results_relative_expression_plot.png  (if --reference is given; disable with --no-plot)
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import statistics
import sys
import zipfile
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional, Tuple

NS_MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
NS_REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def col_letters_to_index(col: str) -> int:
    n = 0
    for ch in col.upper():
        if 'A' <= ch <= 'Z':
            n = n * 26 + (ord(ch) - ord('A') + 1)
    return n


def cell_ref_to_rc(ref: str) -> Tuple[int, int]:
    m = re.match(r"([A-Za-z]+)([0-9]+)", ref)
    if not m:
        raise ValueError(f"Bad cell reference: {ref}")
    return int(m.group(2)), col_letters_to_index(m.group(1))


def find_zip_name(z: zipfile.ZipFile, wanted: str) -> str:
    """Find OOXML member ignoring case and slash/backslash differences."""
    norm_wanted = wanted.replace('\\', '/').lower()
    for name in z.namelist():
        if name.replace('\\', '/').lower() == norm_wanted:
            return name
    raise KeyError(f"Cannot find {wanted} in xlsx archive")


def read_shared_strings(z: zipfile.ZipFile) -> List[str]:
    try:
        name = find_zip_name(z, "xl/sharedStrings.xml")
    except KeyError:
        return []
    root = ET.fromstring(z.read(name))
    out: List[str] = []
    for si in root.findall(NS_MAIN + "si"):
        # Handles normal and rich text shared strings
        text = "".join((t.text or "") for t in si.iter(NS_MAIN + "t"))
        out.append(text)
    return out


def workbook_sheet_map(z: zipfile.ZipFile) -> Dict[str, str]:
    """Return sheet_name -> worksheet xml path."""
    wb_name = find_zip_name(z, "xl/workbook.xml")
    rels_name = find_zip_name(z, "xl/_rels/workbook.xml.rels")
    wb_root = ET.fromstring(z.read(wb_name))
    rels_root = ET.fromstring(z.read(rels_name))

    # relationship id -> target
    rel_map: Dict[str, str] = {}
    for rel in rels_root:
        rid = rel.attrib.get("Id")
        target = rel.attrib.get("Target")
        if rid and target:
            target = target.replace('\\', '/')
            if target.startswith("/"):
                target = target.lstrip("/")
            elif not target.startswith("xl/"):
                target = "xl/" + target.lstrip("/")
            rel_map[rid] = target

    sheet_map: Dict[str, str] = {}
    sheets = wb_root.find(NS_MAIN + "sheets")
    if sheets is None:
        return sheet_map
    for sh in sheets.findall(NS_MAIN + "sheet"):
        name = sh.attrib.get("name")
        rid = sh.attrib.get(NS_REL + "id") or sh.attrib.get("id")
        if name and rid and rid in rel_map:
            sheet_map[name] = rel_map[rid]
    return sheet_map


def parse_cell_value(cell: ET.Element, shared: List[str]) -> Any:
    cell_type = cell.attrib.get("t")

    if cell_type == "inlineStr":
        is_el = cell.find(NS_MAIN + "is")
        if is_el is None:
            return ""
        return "".join((t.text or "") for t in is_el.iter(NS_MAIN + "t"))

    v = cell.find(NS_MAIN + "v")
    if v is None or v.text is None:
        return ""

    raw = v.text
    if cell_type == "s":
        try:
            return shared[int(raw)]
        except Exception:
            return raw
    if cell_type == "b":
        return raw == "1"

    try:
        f = float(raw)
        if f.is_integer():
            return int(f)
        return f
    except Exception:
        return raw


def read_xlsx_sheet(path: str, sheet_name: Optional[str] = None) -> Tuple[str, List[List[Any]]]:
    """Read worksheet into list-of-lists using direct OOXML parsing."""
    with zipfile.ZipFile(path) as z:
        shared = read_shared_strings(z)
        sheet_map = workbook_sheet_map(z)
        if not sheet_map:
            raise RuntimeError("No sheets found in workbook")

        if sheet_name is None:
            # Prefer SYBR if present, otherwise first sheet
            chosen = "SYBR" if "SYBR" in sheet_map else next(iter(sheet_map))
        else:
            # Case-insensitive sheet match
            matches = [s for s in sheet_map if s.lower() == sheet_name.lower()]
            if not matches:
                raise RuntimeError(f"Sheet '{sheet_name}' not found. Available: {list(sheet_map)}")
            chosen = matches[0]

        ws_name = find_zip_name(z, sheet_map[chosen])
        root = ET.fromstring(z.read(ws_name))
        sheet_data = root.find(NS_MAIN + "sheetData")
        if sheet_data is None:
            raise RuntimeError("No sheetData in worksheet")

        cells: Dict[Tuple[int, int], Any] = {}
        max_r = 0
        max_c = 0
        for row in sheet_data.findall(NS_MAIN + "row"):
            for c in row.findall(NS_MAIN + "c"):
                ref = c.attrib.get("r")
                if not ref:
                    continue
                r, col = cell_ref_to_rc(ref)
                val = parse_cell_value(c, shared)
                cells[(r, col)] = val
                max_r = max(max_r, r)
                max_c = max(max_c, col)

        table = []
        for r in range(1, max_r + 1):
            table.append([cells.get((r, c), "") for c in range(1, max_c + 1)])
        return chosen, table


def to_float(x: Any) -> Optional[float]:
    if x is None or x == "":
        return None
    if isinstance(x, (int, float)):
        if math.isfinite(float(x)):
            return float(x)
        return None
    s = str(x).strip().replace(",", ".")
    if s in {"", "NaN", "nan", "N/A", "n/a", "None", "No Ct", "Undetermined", "no ct", "undetermined"}:
        return None
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except Exception:
        return None


def extract_amp_curves(table: List[List[Any]]) -> Dict[str, List[Tuple[float, float]]]:
    """Extract Cycle + well fluorescence columns from Bio-Rad table."""
    # Find header row containing Cycle
    header_i = None
    for i, row in enumerate(table[:20]):
        vals = [str(x).strip().lower() for x in row]
        if "cycle" in vals:
            header_i = i
            break
    if header_i is None:
        raise RuntimeError("Cannot find a header row containing 'Cycle'")

    header = [str(x).strip() for x in table[header_i]]
    try:
        cycle_col = [h.lower() for h in header].index("cycle")
    except ValueError:
        cycle_col = 0

    well_cols = []
    well_re = re.compile(r"^[A-Ha-h](?:[1-9]|1[0-2])$")
    for j, h in enumerate(header):
        if j == cycle_col:
            continue
        if well_re.match(h):
            well_cols.append((j, h.upper()))

    if not well_cols:
        # fallback: all numeric-looking columns except cycle
        for j, h in enumerate(header):
            if j != cycle_col and h:
                well_cols.append((j, h))

    curves: Dict[str, List[Tuple[float, float]]] = {well: [] for _, well in well_cols}
    for row in table[header_i + 1:]:
        if cycle_col >= len(row):
            continue
        cyc = to_float(row[cycle_col])
        if cyc is None:
            continue
        for j, well in well_cols:
            if j < len(row):
                fl = to_float(row[j])
                if fl is not None:
                    curves[well].append((cyc, fl))
    return curves


def linreg(x: List[float], y: List[float]) -> Tuple[float, float, float, float]:
    """Return slope, intercept, r2, residual_sd."""
    n = len(x)
    mx = sum(x) / n
    my = sum(y) / n
    sxx = sum((xi - mx) ** 2 for xi in x)
    syy = sum((yi - my) ** 2 for yi in y)
    if sxx == 0:
        return float('nan'), float('nan'), float('nan'), float('nan')
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    slope = sxy / sxx
    intercept = my - slope * mx
    yhat = [slope * xi + intercept for xi in x]
    ss_res = sum((yi - yh) ** 2 for yi, yh in zip(y, yhat))
    r2 = 1.0 - ss_res / syy if syy > 0 else float('nan')
    resid_sd = math.sqrt(ss_res / max(n - 2, 1))
    return slope, intercept, r2, resid_sd


def median(vals: List[float]) -> float:
    return statistics.median(vals) if vals else float('nan')


def mad(vals: List[float], center: Optional[float] = None) -> float:
    if not vals:
        return float('nan')
    c = median(vals) if center is None else center
    return median([abs(v - c) for v in vals])


def choose_baseline(points: List[Tuple[float, float]], max_baseline_cycle: int = 15) -> Dict[str, Any]:
    """Automatic early-cycle baseline estimate. Uses longest stable early region before signal rises."""
    points = sorted(points)
    early = [(c, f) for c, f in points if c <= max_baseline_cycle]
    if len(early) < 5:
        early = points[:min(10, len(points))]

    # Default: median of cycles 3-10 if present; else first 10 cycles
    preferred = [f for c, f in points if 3 <= c <= 10]
    if len(preferred) < 5:
        preferred = [f for c, f in points[:min(10, len(points))]]

    base = median(preferred)
    noise = 1.4826 * mad(preferred, base)
    if not math.isfinite(noise) or noise == 0:
        # robust fallback from first differences
        diffs = [abs(points[i+1][1] - points[i][1]) for i in range(min(10, len(points)-1))]
        noise = median(diffs) if diffs else 1.0
    if not math.isfinite(noise) or noise == 0:
        noise = max(abs(base) * 0.001, 1e-9)

    return {
        "baseline": base,
        "baseline_cycles": "3-10" if len([1 for c, _ in points if 3 <= c <= 10]) >= 5 else f"first_{len(preferred)}",
        "noise": noise,
    }




def _smooth_values(vals: List[float], radius: int = 2) -> List[float]:
    """Small dependency-free moving-average smoother used only for QC."""
    out: List[float] = []
    n = len(vals)
    for i in range(n):
        lo = max(0, i - radius)
        hi = min(n, i + radius + 1)
        out.append(sum(vals[lo:hi]) / max(hi - lo, 1))
    return out


def _step_slopes(cycles: List[float], logs: List[float]) -> List[float]:
    out: List[float] = []
    for i in range(len(cycles) - 1):
        dx = cycles[i + 1] - cycles[i]
        if dx != 0:
            out.append((logs[i + 1] - logs[i]) / dx)
    return out


def _sd(vals: List[float]) -> float:
    vals = [float(v) for v in vals if v is not None and math.isfinite(float(v))]
    if len(vals) < 2:
        return 0.0
    m = sum(vals) / len(vals)
    return math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1))


def _linreg_slopes_by_halves(cycles: List[float], logs: List[float]) -> Tuple[float, float, float]:
    """Return lower-half slope, upper-half slope, and abs difference.

    This implements the LinRegPCR/Ruijter-like baseline criterion: after baseline
    subtraction, the lower and upper halves of the log-linear segment should have
    nearly the same slope. With 4 data points, the two halves are the first and
    last two points; for longer windows, overlapping halves are avoided.
    """
    n = len(cycles)
    if n < 4:
        return float('nan'), float('nan'), float('nan')
    mid = n // 2
    x_low = cycles[:mid]
    y_low = logs[:mid]
    x_high = cycles[n - mid:]
    y_high = logs[n - mid:]
    sl_low, _, _, _ = linreg(x_low, y_low)
    sl_high, _, _, _ = linreg(x_high, y_high)
    return sl_low, sl_high, abs(sl_high - sl_low)


def _quantile(vals: List[float], q: float) -> Optional[float]:
    vals = sorted([v for v in vals if v is not None and math.isfinite(v)])
    if not vals:
        return None
    if len(vals) == 1:
        return vals[0]
    k = (len(vals) - 1) * q
    lo = math.floor(k)
    hi = math.ceil(k)
    if lo == hi:
        return vals[int(k)]
    return vals[lo] * (hi - k) + vals[hi] * (k - lo)


def _baseline_balance_metrics(
    cycles: List[float],
    raw_vals: List[float],
    baseline: float,
) -> Optional[Dict[str, float]]:
    """Metrics for one candidate baseline on one putative exponential window.

    This is the core LinRegPCR-like idea used below: after baseline subtraction,
    the lower and upper halves of the log-linear window should have nearly the
    same slope. If the baseline is too high or too low, the log-transformed curve
    becomes visibly bent and the two half-window slopes diverge.
    """
    fcs = [float(f) - float(baseline) for f in raw_vals]
    if any((not math.isfinite(fc)) or fc <= 0 for fc in fcs):
        return None
    logs = [math.log10(fc) for fc in fcs]
    slope, intercept, r2, resid_sd = linreg(cycles, logs)
    if not (math.isfinite(slope) and slope > 0 and math.isfinite(r2)):
        return None

    step_slopes = _step_slopes(cycles, logs)
    step_sd = _sd(step_slopes)
    rel_step_sd = (step_sd / slope) if slope > 0 else 999.0
    positive_steps = sum(1 for st in step_slopes if math.isfinite(st) and st > 0)
    frac_positive_steps = positive_steps / max(len(step_slopes), 1)

    sl_low, sl_high, balance = _linreg_slopes_by_halves(cycles, logs)
    rel_balance = (balance / slope) if slope > 0 and math.isfinite(balance) else 999.0
    e_amp = 10.0 ** slope

    return {
        "slope": slope,
        "intercept": intercept,
        "r2": r2,
        "residual_sd": resid_sd,
        "E_amp_factor": e_amp,
        "relative_slope_balance": rel_balance,
        "relative_step_sd": rel_step_sd,
        "frac_positive_steps": frac_positive_steps,
    }


def _optimize_baseline_for_window(
    cycles: List[float],
    raw_vals: List[float],
    early_base: float,
    noise: float,
    raw_range: float,
    min_e: float,
    max_e: float,
) -> Optional[Tuple[float, Dict[str, float]]]:
    """Iteratively find the baseline that best straightens one window.

    This is not copied LinRegPCR source code. It implements the same published
    principle in a transparent way: search for a baseline that minimizes the
    difference between the lower-half and upper-half slopes of the selected
    log-linear window, while penalizing unstable or implausible fits.
    """
    if len(cycles) < 4 or len(cycles) != len(raw_vals):
        return None

    min_raw_in_window = min(raw_vals)
    eps = max(abs(min_raw_in_window) * 1e-12, noise * 1e-6, raw_range * 1e-12, 1e-12)

    # Upper bound must keep every corrected fluorescence in the window positive.
    upper = min(min_raw_in_window - eps, early_base + 0.10 * raw_range)
    # Allow baseline to move substantially below the early-cycle background. This
    # mimics LinRegPCR's stepwise lowering of baseline when the log curve is bent.
    lower = min(upper - eps, early_base - 0.35 * raw_range)
    # But do not allow a wildly negative baseline on already baseline-corrected
    # exports; that can manufacture unrealistically straight log windows.
    raw_min = min(raw_vals)
    lower = max(lower, raw_min - 0.50 * raw_range)

    if not (math.isfinite(lower) and math.isfinite(upper)) or upper <= lower:
        return None

    best_baseline: Optional[float] = None
    best_metrics: Optional[Dict[str, float]] = None
    best_objective = float("inf")

    lo, hi = lower, upper
    # Coarse pass + repeated local refinements around the best value.
    for level, n_grid in enumerate((121, 81, 81, 61)):
        if hi <= lo:
            break
        for i in range(n_grid):
            b = lo + (hi - lo) * i / max(n_grid - 1, 1)
            m = _baseline_balance_metrics(cycles, raw_vals, b)
            if m is None:
                continue

            e = m["E_amp_factor"]
            e_penalty = 0.0
            # The strict E gate is applied later. Here use a soft penalty so the
            # optimizer does not pick a mathematically straight but impossible fit.
            if e < min_e:
                e_penalty = (min_e - e) / max(min_e, 1e-9)
            elif e > max_e:
                e_penalty = (e - max_e) / max(max_e, 1e-9)

            objective = (
                m["relative_slope_balance"]
                + 0.35 * min(m["relative_step_sd"], 5.0)
                + 2.50 * max(0.0, 0.995 - m["r2"])
                + 0.50 * max(0.0, 0.80 - m["frac_positive_steps"])
                + 0.75 * e_penalty
            )
            if objective < best_objective:
                best_objective = objective
                best_baseline = b
                best_metrics = dict(m)
                best_metrics["baseline_objective"] = objective

        if best_baseline is None:
            break
        span = (hi - lo) / (7.0 if level == 0 else 5.0)
        lo = max(lower, best_baseline - span / 2.0)
        hi = min(upper, best_baseline + span / 2.0)

    if best_baseline is None or best_metrics is None:
        return None
    return best_baseline, best_metrics


def _baseline_candidate_values(
    points: List[Tuple[float, float]],
    early_base: float,
    noise: float,
    min_window: int,
    max_window: int,
    min_e: float,
    max_e: float,
) -> List[float]:
    """Generate baseline candidates using iterative window straightening.

    Older versions used a small grid around the early-cycle median. This version
    first finds rough exponential-window candidates and then, for each such
    window, iteratively adjusts baseline until the lower and upper halves of the
    log-linear window have balanced slopes. The final strict window selection is
    still performed in calc_efficiency_for_well().
    """
    points = sorted(points)
    raws_all = [f for _, f in points]
    raw_min = min(raws_all)
    raw_max = max(raws_all)
    raw_range = max(raw_max - raw_min, abs(raw_max) * 1e-9, 1e-9)

    base_candidates: List[float] = [early_base]
    # Conservative local candidates keep the method usable when no rough window
    # can be optimized, e.g. almost flat or very noisy curves.
    for mult in (-3, -2, -1, -0.5, 0.5, 1, 2, 3):
        base_candidates.append(early_base + mult * noise)
    for frac in (-0.010, -0.005, 0.005, 0.010):
        base_candidates.append(early_base + frac * raw_range)

    # Rough pass with the early baseline only: identify windows that plausibly
    # represent rising signal and are worth baseline optimization.
    corrected = [(c, f, f - early_base) for c, f in points]
    positive = [fc for _, _, fc in corrected if fc > 0]
    if len(positive) < min_window:
        return sorted({round(b, 12) for b in base_candidates if math.isfinite(b)})

    rough_max_delta = max(positive)
    rough_low = max(6.0 * noise, 0.004 * rough_max_delta, abs(early_base) * 0.0002, 1e-12)
    rough_high = 0.92 * rough_max_delta
    log_points = [(c, f, fc, math.log10(fc)) for c, f, fc in corrected if fc > 0]

    optimized: List[Tuple[float, float, Dict[str, float]]] = []
    seen_windows = set()
    for window_n in range(max(4, min_window), max_window + 1):
        if len(log_points) < window_n:
            continue
        for start_idx in range(0, len(log_points) - window_n + 1):
            window = log_points[start_idx:start_idx + window_n]
            cycles = [c for c, _, _, _ in window]
            if not all(abs((cycles[i + 1] - cycles[i]) - 1.0) <= 1e-6 for i in range(len(cycles) - 1)):
                continue
            key = (cycles[0], cycles[-1], window_n)
            if key in seen_windows:
                continue
            seen_windows.add(key)

            fcs = [fc for _, _, fc, _ in window]
            # Relaxed rough gate: do not over-filter before baseline optimization.
            if median(fcs) < rough_low or median(fcs) > rough_high:
                continue
            logs = [lg for _, _, _, lg in window]
            slope, _, r2, _ = linreg(cycles, logs)
            if not (math.isfinite(slope) and slope > 0 and math.isfinite(r2)):
                continue
            e_rough = 10.0 ** slope
            step_slopes = _step_slopes(cycles, logs)
            frac_pos = sum(1 for st in step_slopes if math.isfinite(st) and st > 0) / max(len(step_slopes), 1)
            if r2 < 0.965 or e_rough < 1.03 or e_rough > 3.20 or frac_pos < 0.60:
                continue

            raw_vals = [f for _, f, _, _ in window]
            opt = _optimize_baseline_for_window(
                cycles, raw_vals, early_base, noise, raw_range, min_e, max_e
            )
            if opt is None:
                continue
            b_opt, metrics = opt
            # Prefer baselines that make a stable, balanced, high-r2 window.
            rank = (
                metrics.get("baseline_objective", 999.0)
                + 0.15 * abs(metrics.get("E_amp_factor", 2.0) - 1.80)
                + 0.005 * min(cycles)
            )
            optimized.append((rank, b_opt, metrics))

    optimized.sort(key=lambda x: x[0])
    local_step = max(0.25 * noise, 0.0005 * raw_range, 1e-12)
    for _, b, _ in optimized[:45]:
        base_candidates.extend([b, b - local_step, b + local_step])

    # Keep the list finite and deterministic. Existing strict filters decide
    # whether a candidate baseline/window combination is acceptable.
    unique = sorted({round(b, 12) for b in base_candidates if math.isfinite(float(b))})
    if len(unique) > 160:
        # Preserve candidates closest to optimized baselines and early baseline.
        anchors = [early_base] + [b for _, b, _ in optimized[:20]]
        unique = sorted(unique, key=lambda b: min(abs(b - a) for a in anchors))[:160]
    return unique


def _rough_noise_and_base(points: List[Tuple[float, float]]) -> Tuple[float, float, str]:
    preferred = [f for c, f in points if 3 <= c <= 10]
    label = "3-10"
    if len(preferred) < 5:
        preferred = [f for _, f in points[:min(10, len(points))]]
        label = f"first_{len(preferred)}"
    base = median(preferred)
    noise = 1.4826 * mad(preferred, base)
    if not math.isfinite(noise) or noise == 0:
        diffs = [abs(points[i + 1][1] - points[i][1]) for i in range(min(10, len(points) - 1))]
        noise = median(diffs) if diffs else 1.0
    if not math.isfinite(noise) or noise == 0:
        noise = max(abs(base) * 0.001, 1e-9)
    return base, noise, label


def calc_efficiency_for_well(
    points: List[Tuple[float, float]],
    min_window: int = 5,
    max_window: int = 7,
    min_e: float = 1.20,
    max_e: float = 2.40,
    good_r2: float = 0.990,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Strict-gated LinReg/RDML-like efficiency estimate for one qPCR well.

    Important principle:
      A window is never allowed to be OK only because it has a nice E/R2.
      It must first be inside the real amplification signal region:
        corrected fluorescence > noise / low-amplitude gate,
        corrected fluorescence < plateau gate.

    Selection order:
      1) robust early baseline and noise;
      2) hard dynamic-range gates: above noise, below plateau;
      3) consecutive candidate windows of min_window..max_window cycles;
      4) strict OK criteria: E range, R2, slope-balance, stable step slopes;
      5) if no strict window exists, return a diagnostic BAD_FIT window, but not OK.

    Expression rows are not deleted because of BAD_FIT; however, BAD_FIT E values
    should not be used for target mean efficiency in rdml_like mode.
    """
    points = sorted([(float(c), float(f)) for c, f in points if c is not None and f is not None])
    well_curve: List[Dict[str, Any]] = []
    candidate_window_rows: List[Dict[str, Any]] = []

    min_window = max(4, int(min_window))
    max_window = max(min_window, int(max_window))

    def empty_result(status: str, details: str, base: Optional[float] = None, noise_val: Optional[float] = None,
                     max_delta_val: Optional[float] = None, low_val: Optional[float] = None,
                     high_val: Optional[float] = None) -> Dict[str, Any]:
        return {
            "status": status,
            "E_amp_factor": None,
            "Efficiency_percent": None,
            "slope": None,
            "intercept": None,
            "r2": None,
            "residual_sd": None,
            "cycle_start": None,
            "cycle_end": None,
            "window_cycles": None,
            "window_n": 0,
            "baseline": round(base, 6) if base is not None and math.isfinite(float(base)) else None,
            "baseline_cycles": None,
            "noise": round(noise_val, 6) if noise_val is not None and math.isfinite(float(noise_val)) else None,
            "max_delta": round(max_delta_val, 6) if max_delta_val is not None and math.isfinite(float(max_delta_val)) else None,
            "threshold_low": round(low_val, 6) if low_val is not None and math.isfinite(float(low_val)) else None,
            "threshold_high": round(high_val, 6) if high_val is not None and math.isfinite(float(high_val)) else None,
            "derivative_peak_cycle": None,
            "derivative_peak_slope": None,
            "slope_quality_vs_peak": None,
            "self_check_status": "FAIL",
            "self_check_details": details,
            "baseline_slope_low": None,
            "baseline_slope_high": None,
            "baseline_slope_balance": None,
        }

    if len(points) < min_window + 3:
        return empty_result("TOO_FEW_POINTS", "too_few_points"), well_curve, candidate_window_rows

    early_base, noise, baseline_label = _rough_noise_and_base(points)
    raw_vals = [f for _, f in points]
    raw_min, raw_max = min(raw_vals), max(raw_vals)
    raw_range = max(raw_max - raw_min, 1e-9)

    # LinRegPCR-like baseline optimization. Instead of using only a small grid
    # around the early-cycle baseline, generate candidate baselines by iteratively
    # straightening plausible log-linear windows: the best baseline is the one for
    # which lower-half and upper-half slopes of the window are most balanced.
    baseline_candidates = _baseline_candidate_values(
        points,
        early_base=early_base,
        noise=noise,
        min_window=min_window,
        max_window=max_window,
        min_e=min_e,
        max_e=max_e,
    )
    baseline_candidates = sorted(
        {round(b, 12) for b in baseline_candidates if math.isfinite(float(b))},
        key=lambda b: abs(float(b) - early_base),
    )

    strict_best: Optional[Dict[str, Any]] = None
    diagnostic_best: Optional[Dict[str, Any]] = None

    for baseline in baseline_candidates:
        corrected = [(c, f, f - baseline) for c, f in points]
        positive = [fc for _, _, fc in corrected if fc > 0]
        if len(positive) < min_window:
            continue
        max_delta = max(positive)
        if not math.isfinite(max_delta) or max_delta <= 0:
            continue

        # HARD gates. A candidate point below low is treated as noise; above high as plateau.
        # Do not allow such windows to become OK even when R2 looks good.
        low = max(10.0 * noise, 0.01 * max_delta, abs(baseline) * 0.0005, 1e-12)
        high = 0.80 * max_delta
        no_amp_threshold = max(10.0 * noise, abs(baseline) * 0.005, 1e-9)
        if max_delta <= no_amp_threshold:
            continue

        log_points = [(c, fc, math.log10(fc)) for c, _, fc in corrected if fc > 0]

        # Diagnostic curve rows for the first/closest baseline only.
        if not well_curve and abs(float(baseline) - early_base) == min(abs(float(b) - early_base) for b in baseline_candidates):
            for c, raw, fc in corrected:
                well_curve.append({
                    "Cycle": c,
                    "Raw_Fluorescence": raw,
                    "Baseline": baseline,
                    "Corrected_Fluorescence": fc,
                    "Log10_Corrected": math.log10(fc) if fc > 0 else None,
                    "Selected_For_Window": False,
                })

        for window_n in range(min_window, max_window + 1):
            if len(log_points) < window_n:
                continue
            for start_idx in range(0, len(log_points) - window_n + 1):
                window = log_points[start_idx:start_idx + window_n]
                cycles = [c for c, _, _ in window]
                fcs = [fc for _, fc, _ in window]
                logs = [lg for _, _, lg in window]
                window_cycles = ",".join(str(int(c)) if float(c).is_integer() else str(c) for c in cycles)

                # Must be truly consecutive cycles. No spanning over missing data.
                consecutive = all(abs((cycles[i + 1] - cycles[i]) - 1.0) <= 1e-6 for i in range(len(cycles) - 1))
                if not consecutive:
                    candidate_window_rows.append({
                        "baseline": baseline, "cycle_start": min(cycles), "cycle_end": max(cycles),
                        "window_cycles": window_cycles, "window_n": window_n, "status": "REJECTED_GAP",
                        "E_amp_factor": None, "Efficiency_percent": None, "r2": None,
                        "slope": None, "score": None,
                    })
                    continue

                in_dynamic_range = all((fc > low and fc < high) for fc in fcs)
                slope, intercept, r2, resid_sd = linreg(cycles, logs)
                if not (math.isfinite(slope) and math.isfinite(r2) and slope > 0):
                    continue
                E = 10 ** slope
                eff = (E - 1.0) * 100.0

                step_slopes = _step_slopes(cycles, logs)
                positive_steps = sum(1 for s in step_slopes if math.isfinite(float(s)) and s > 0)
                frac_positive_steps = positive_steps / max(len(step_slopes), 1)
                step_sd = _sd(step_slopes)
                rel_step_sd = (step_sd / slope) if slope > 0 else 999.0
                sl_low, sl_high, balance = _linreg_slopes_by_halves(cycles, logs)
                rel_balance = (balance / slope) if slope > 0 and math.isfinite(balance) else 999.0
                center_frac = median(fcs) / max_delta if max_delta > 0 else 0.0

                # Strict biological/technical validity of the selected window.
                passes_noise_plateau = in_dynamic_range
                passes_e_r2 = (E >= min_e and E <= max_e and r2 >= good_r2)
                passes_balance = rel_balance <= 0.30
                passes_step_stability = rel_step_sd <= 0.35 and frac_positive_steps >= 0.80
                strict_ok = passes_noise_plateau and passes_e_r2 and passes_balance and passes_step_stability

                # Score among candidates. Dynamic range and strict filters are not merely soft;
                # strict_ok is required for OK. Score only ranks windows within a class.
                score = (
                    1000.0 * (r2 if math.isfinite(r2) else -1.0)
                    + 4.0 * window_n
                    - 70.0 * min(rel_balance, 2.0)
                    - 45.0 * min(rel_step_sd, 2.0)
                    - 25.0 * abs(center_frac - 0.32)
                    - 0.20 * min(cycles)
                )
                if not passes_noise_plateau:
                    score -= 10000.0
                if not (E >= min_e and E <= max_e):
                    score -= 3000.0
                if r2 < good_r2:
                    score -= 2000.0 * min(good_r2 - r2, 0.2)

                fail_parts = []
                if not passes_noise_plateau:
                    fail_parts.append("below_noise_or_on_plateau")
                if E < min_e:
                    fail_parts.append("low_E")
                if E > max_e:
                    fail_parts.append("high_E")
                if r2 < good_r2:
                    fail_parts.append("low_R2")
                if not passes_balance:
                    fail_parts.append("unbalanced_window_halves")
                if not passes_step_stability:
                    fail_parts.append("unstable_or_nonmonotonic_step_slopes")

                cand = {
                    "baseline": baseline,
                    "baseline_cycles": baseline_label,
                    "noise": noise,
                    "max_delta": max_delta,
                    "threshold_low": low,
                    "threshold_high": high,
                    "window_n": window_n,
                    "cycles": cycles,
                    "cycle_start": min(cycles),
                    "cycle_end": max(cycles),
                    "window_cycles": window_cycles,
                    "slope": slope,
                    "intercept": intercept,
                    "r2": r2,
                    "residual_sd": resid_sd,
                    "E_amp_factor": E,
                    "Efficiency_percent": eff,
                    "baseline_slope_low": sl_low,
                    "baseline_slope_high": sl_high,
                    "baseline_slope_balance": balance,
                    "relative_slope_balance": rel_balance,
                    "relative_step_sd": rel_step_sd,
                    "frac_positive_steps": frac_positive_steps,
                    "center_frac": center_frac,
                    "in_dynamic_range": in_dynamic_range,
                    "strict_ok": strict_ok,
                    "score": score,
                    "fail_details": ";".join(fail_parts) if fail_parts else "strict_window_ok",
                }
                candidate_window_rows.append({
                    "baseline": baseline,
                    "cycle_start": min(cycles),
                    "cycle_end": max(cycles),
                    "window_cycles": window_cycles,
                    "window_n": window_n,
                    "slope": slope,
                    "E_amp_factor": E,
                    "Efficiency_percent": eff,
                    "r2": r2,
                    "residual_sd": resid_sd,
                    "relative_slope_balance": rel_balance,
                    "relative_step_sd": rel_step_sd,
                    "frac_positive_steps": frac_positive_steps,
                    "center_frac": center_frac,
                    "in_dynamic_range": in_dynamic_range,
                    "strict_ok": strict_ok,
                    "status": "ACCEPTABLE" if strict_ok else "REJECTED_CRITERIA",
                    "fail_details": cand["fail_details"],
                    "score": score,
                })

                if strict_ok:
                    if strict_best is None or cand["score"] > strict_best["score"]:
                        strict_best = cand
                # Diagnostic window: prefer dynamic-range windows first, then highest score.
                diagnostic_rank = (1 if in_dynamic_range else 0, score)
                if diagnostic_best is None or diagnostic_rank > diagnostic_best.get("diagnostic_rank", (-1, -float("inf"))):
                    cand["diagnostic_rank"] = diagnostic_rank
                    diagnostic_best = cand

    selected = strict_best or diagnostic_best
    if selected is None:
        # Use early baseline diagnostics if nothing usable exists.
        max_delta0 = max([f - early_base for _, f in points], default=float("nan"))
        return empty_result(
            "NO_LOG_WINDOW",
            "no_candidate_window_after_noise_plateau_filter",
            base=early_base,
            noise_val=noise,
            max_delta_val=max_delta0,
            low_val=None,
            high_val=None,
        ), well_curve, candidate_window_rows

    # Rebuild well_curve for the selected baseline so debug output matches the chosen result.
    baseline = float(selected["baseline"])
    selected_set = set(float(c) for c in selected["cycles"])
    well_curve = []
    log_points_for_peak: List[Tuple[float, float]] = []
    for c, f in points:
        fc = f - baseline
        log_fc = math.log10(fc) if fc > 0 else None
        if log_fc is not None and math.isfinite(log_fc):
            log_points_for_peak.append((c, log_fc))
        well_curve.append({
            "Cycle": c,
            "Raw_Fluorescence": f,
            "Baseline": baseline,
            "Corrected_Fluorescence": fc,
            "Log10_Corrected": log_fc,
            "Selected_For_Window": float(c) in selected_set,
        })

    derivative_peak_cycle = None
    derivative_peak_slope = None
    if len(log_points_for_peak) >= 2:
        slopes = []
        for i in range(len(log_points_for_peak) - 1):
            c0, l0 = log_points_for_peak[i]
            c1, l1 = log_points_for_peak[i + 1]
            if c1 != c0:
                s = (l1 - l0) / (c1 - c0)
                if math.isfinite(s) and s > 0:
                    slopes.append(((c0 + c1) / 2.0, s))
        if slopes:
            derivative_peak_cycle, derivative_peak_slope = max(slopes, key=lambda x: x[1])

    E = float(selected["E_amp_factor"])
    r2 = float(selected["r2"])
    rel_balance = float(selected["relative_slope_balance"])
    rel_step_sd = float(selected["relative_step_sd"])
    slope_quality_vs_peak = None
    if derivative_peak_slope is not None and derivative_peak_slope > 0:
        slope_quality_vs_peak = float(selected["slope"]) / derivative_peak_slope

    if selected.get("strict_ok"):
        status = "OK"
        self_check_status = "PASS"
        self_check_details = "strict_window_ok;above_noise_below_plateau;E_R2_balance_step_checks_passed"
    else:
        status = "BAD_FIT"
        self_check_status = "FAIL"
        self_check_details = selected.get("fail_details", "no_strict_window_found")

    result = {
        "status": status,
        "E_amp_factor": round(E, 6),
        "Efficiency_percent": round(float(selected["Efficiency_percent"]), 3),
        "slope": round(float(selected["slope"]), 6),
        "intercept": round(float(selected["intercept"]), 6),
        "r2": round(r2, 6),
        "residual_sd": round(float(selected["residual_sd"]), 8),
        "cycle_start": selected["cycle_start"],
        "cycle_end": selected["cycle_end"],
        "window_cycles": selected["window_cycles"],
        "window_n": selected["window_n"],
        "baseline": round(baseline, 6),
        "baseline_cycles": selected["baseline_cycles"],
        "noise": round(float(selected["noise"]), 6),
        "max_delta": round(float(selected["max_delta"]), 6),
        "threshold_low": round(float(selected["threshold_low"]), 6),
        "threshold_high": round(float(selected["threshold_high"]), 6),
        "derivative_peak_cycle": derivative_peak_cycle,
        "derivative_peak_slope": round(derivative_peak_slope, 6) if derivative_peak_slope is not None else None,
        "slope_quality_vs_peak": round(slope_quality_vs_peak, 6) if slope_quality_vs_peak is not None else None,
        "self_check_status": self_check_status,
        "self_check_details": self_check_details,
        "baseline_slope_low": round(float(selected["baseline_slope_low"]), 6) if math.isfinite(float(selected["baseline_slope_low"])) else None,
        "baseline_slope_high": round(float(selected["baseline_slope_high"]), 6) if math.isfinite(float(selected["baseline_slope_high"])) else None,
        "baseline_slope_balance": round(float(selected["baseline_slope_balance"]), 6) if math.isfinite(float(selected["baseline_slope_balance"])) else None,
    }
    return result, well_curve, candidate_window_rows

def safe_out_prefix(path: str, prefix: Optional[str]) -> str:
    if prefix:
        return prefix
    base = os.path.splitext(os.path.basename(path))[0]
    base = re.sub(r"[^A-Za-z0-9А-Яа-я_.-]+", "_", base).strip("_")
    return base + "_linreg"


def write_csv(path: str, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)




def normalize_well(well: Any) -> str:
    s = str(well).strip().upper()
    m = re.match(r"^([A-P])0*([0-9]{1,2})$", s)
    if m:
        return f"{m.group(1)}{int(m.group(2)):02d}"
    return s


def find_header_row_for_columns(table: List[List[Any]], required: List[str], max_scan: int = 40) -> int:
    required_l = [x.strip().lower() for x in required]
    best_i = -1
    best_score = -1
    for i, row in enumerate(table[:max_scan]):
        vals = [str(x).strip().lower() for x in row]
        score = sum(1 for r in required_l if r in vals)
        # extra signal for Bio-Rad Cq table
        score += sum(1 for v in vals if v in {"well", "cq", "target", "sample"})
        if score > best_score:
            best_score = score
            best_i = i
    if best_score <= 0:
        raise RuntimeError(f"Cannot find header row containing any of: {required}")
    return best_i


def table_to_records(table: List[List[Any]], required_cols: List[str]) -> Tuple[List[str], List[Dict[str, Any]]]:
    hi = find_header_row_for_columns(table, required_cols)
    raw_header = [str(x).strip() for x in table[hi]]
    header: List[str] = []
    seen: Dict[str, int] = {}
    for idx, h in enumerate(raw_header):
        if not h:
            h = f"_blank_{idx+1}"
        if h in seen:
            seen[h] += 1
            h = f"{h}_{seen[h]}"
        else:
            seen[h] = 1
        header.append(h)
    rows: List[Dict[str, Any]] = []
    for row in table[hi+1:]:
        rec = {header[j]: (row[j] if j < len(row) else "") for j in range(len(header))}
        if any(str(v).strip() != "" for v in rec.values()):
            rows.append(rec)
    return header, rows


def get_col_ci(row: Dict[str, Any], col: Optional[str], default: Any = "") -> Any:
    if not col:
        return default
    if col in row:
        return row[col]
    col_l = col.strip().lower()
    for k, v in row.items():
        if k.strip().lower() == col_l:
            return v
    return default


def has_col_ci(header: List[str], col: str) -> bool:
    col_l = col.strip().lower()
    return any(h.strip().lower() == col_l for h in header)


def mean(vals: List[float]) -> Optional[float]:
    vals = [v for v in vals if v is not None and math.isfinite(v)]
    return sum(vals) / len(vals) if vals else None


def stdev(vals: List[float]) -> Optional[float]:
    vals = [v for v in vals if v is not None and math.isfinite(v)]
    if len(vals) < 2:
        return None
    m = sum(vals) / len(vals)
    return math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1))


def write_csv_dynamic(path: str, rows: List[Dict[str, Any]], preferred: Optional[List[str]] = None) -> None:
    fields: List[str] = []
    if preferred:
        for f in preferred:
            if f not in fields:
                fields.append(f)
    for row in rows:
        for k in row.keys():
            if k not in fields:
                fields.append(k)
    if not fields:
        fields = ["empty"]
    write_csv(path, rows, fields)



def plot_relative_expression_summary(
    summary_rows: List[Dict[str, Any]],
    out_path: str,
    metric: str = "relative",
    error: str = "sem",
    title: Optional[str] = None,
) -> bool:
    """
    Save a grouped bar plot from final_summary_by_condition rows.

    metric:
      relative -> mean_relative_expression
      log2     -> mean_log2_relative_expression
    error:
      sem, sd, none

    The function imports matplotlib lazily, so the numerical pipeline still works
    on computers where matplotlib is not installed.
    """
    rows = [r for r in summary_rows if str(r.get("Target", "")).strip() and str(r.get("Condition", "")).strip()]
    if not rows:
        print("WARNING: plot skipped: summary table is empty.", file=sys.stderr)
        return False

    if metric == "log2":
        mean_col = "mean_log2_relative_expression"
        sd_col = "sd_log2_relative_expression"
        sem_col = "sem_log2_relative_expression"
        y_label = "log2(relative expression)"
        default_title = "Efficiency-corrected qPCR expression, log2 scale"
    else:
        mean_col = "mean_relative_expression"
        sd_col = "sd_relative_expression"
        sem_col = "sem_relative_expression"
        y_label = "Relative expression"
        default_title = "Efficiency-corrected qPCR relative expression"

    err_col = None
    if error == "sem":
        err_col = sem_col
    elif error == "sd":
        err_col = sd_col

    targets = sorted({str(r.get("Target", "")).strip() for r in rows})
    conditions = sorted({str(r.get("Condition", "")).strip() for r in rows})

    # Try to keep the control condition first when it is present.
    controls = [str(r.get("control_condition", "")).strip() for r in rows if str(r.get("control_condition", "")).strip()]
    control = controls[0] if controls else None
    if control in conditions:
        conditions = [control] + [c for c in conditions if c != control]

    data: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for r in rows:
        t = str(r.get("Target", "")).strip()
        c = str(r.get("Condition", "")).strip()
        val = to_float(r.get(mean_col))
        err_val = to_float(r.get(err_col)) if err_col else None
        if val is not None and math.isfinite(val):
            data[(t, c)] = {"mean": val, "error": err_val if err_val is not None and math.isfinite(err_val) else 0.0}

    if not data:
        print(f"WARNING: plot skipped: no finite values in column {mean_col}.", file=sys.stderr)
        return False

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as e:
        print("WARNING: plot skipped: matplotlib/numpy is not available. Install with: pip install matplotlib numpy", file=sys.stderr)
        print(f"Matplotlib import error: {e}", file=sys.stderr)
        return False

    n_targets = len(targets)
    n_conditions = len(conditions)
    x = np.arange(n_targets, dtype=float)
    total_width = 0.82
    bar_width = total_width / max(n_conditions, 1)

    fig_width = max(7.0, 1.15 * n_targets * max(n_conditions, 2))
    fig_height = 5.2
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=300)

    for j, condition in enumerate(conditions):
        offsets = x - total_width / 2 + bar_width / 2 + j * bar_width
        means = [data.get((target, condition), {}).get("mean", float("nan")) for target in targets]
        errors = [data.get((target, condition), {}).get("error", 0.0) for target in targets]
        yerr = errors if err_col else None
        ax.bar(offsets, means, width=bar_width, label=condition, yerr=yerr, capsize=4 if err_col else 0)

    ax.set_xticks(x)
    ax.set_xticklabels(targets)
    ax.set_ylabel(y_label)
    ax.set_title(title or default_title)
    ax.legend(title="Condition", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.axhline(0 if metric == "log2" else 1, linewidth=1, linestyle="--")

    if metric == "relative":
        ax.set_ylim(bottom=0)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return True


def amplification_sheet_names(path: str) -> List[str]:
    """Return sheet names from workbook."""
    with zipfile.ZipFile(path) as z:
        return list(workbook_sheet_map(z).keys())


def compute_per_well_efficiency(args: argparse.Namespace) -> List[Dict[str, Any]]:
    """
    Compute per-well efficiency.

    In older Bio-Rad exports all amplification curves may be on SYBR.
    In newer/other exports curves may be split across target sheets: tfeb, PPIA, tbp, cyc1, etc.
    Therefore --sheet auto/all reads all non-Run Information sheets and merges wells.
    """
    requested = str(args.sheet or "auto").strip()
    if requested.lower() in {"auto", "all", "*"}:
        all_sheets = amplification_sheet_names(args.amp)
        sheet_list = [s for s in all_sheets if "run information" not in s.lower()]
    else:
        sheet_list = [requested]

    results_by_well: Dict[str, Dict[str, Any]] = {}
    seen_sources: Dict[str, str] = {}

    for sheet_name in sheet_list:
        sheet, table = read_xlsx_sheet(args.amp, sheet_name)
        try:
            curves = extract_amp_curves(table)
        except Exception as e:
            print(f"WARNING: amplification sheet '{sheet}' skipped: {e}", file=sys.stderr)
            continue
        if not curves:
            print(f"WARNING: amplification sheet '{sheet}' has no well curves; skipped.", file=sys.stderr)
            continue

        print(f"Read amplification sheet: {sheet}; wells found: {len(curves)}")

        sorted_wells = sorted(curves.keys(), key=lambda w: (normalize_well(w)[0], int(re.sub(r"[^0-9]", "", normalize_well(w)) or 0), normalize_well(w)))
        for idx, well in enumerate(sorted_wells, start=1):
            norm_well = normalize_well(well)
            if idx == 1 or idx % 12 == 0 or idx == len(sorted_wells):
                print(f"  efficiency window search: {idx}/{len(sorted_wells)} wells processed", flush=True)

            res, _corr, _cand = calc_efficiency_for_well(
                curves[well],
                min_window=args.min_window,
                max_window=args.max_window,
                min_e=args.min_e,
                max_e=args.max_e,
                good_r2=args.good_r2,
            )
            res = {"Well": norm_well, "Amp_Sheet": sheet, **res}

            # If duplicate well appears, prefer a target-specific sheet over generic SYBR.
            if norm_well in results_by_well:
                old_sheet = seen_sources.get(norm_well, "")
                if old_sheet.lower() == "sybr" and sheet.lower() != "sybr":
                    results_by_well[norm_well] = res
                    seen_sources[norm_well] = sheet
                else:
                    print(
                        f"WARNING: duplicate well {norm_well} in amplification sheets "
                        f"'{old_sheet}' and '{sheet}'. Keeping '{old_sheet}'.",
                        file=sys.stderr,
                    )
            else:
                results_by_well[norm_well] = res
                seen_sources[norm_well] = sheet

    if not results_by_well:
        raise RuntimeError("No well curves found in amplification file")

    results = list(results_by_well.values())
    results.sort(key=lambda r: (str(r.get("Well", ""))[0], int(re.sub(r"[^0-9]", "", str(r.get("Well", ""))) or 0), str(r.get("Well", ""))))
    print(f"Amplification sheets used: {', '.join(sheet_list)}")
    return results

def read_cq_records(path: str, sheet: Optional[str]) -> Tuple[str, List[str], List[Dict[str, Any]]]:
    chosen, table = read_xlsx_sheet(path, sheet)
    header, rows = table_to_records(table, ["Well", "Cq", "Target", "Sample"])
    return chosen, header, rows


def build_merged_cq_efficiency(args: argparse.Namespace, eff_rows: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], List[str]]:
    cq_sheet, header, cq_rows = read_cq_records(args.cq, args.cq_sheet)
    print(f"Read Cq sheet: {cq_sheet}")
    print("Cq file columns: " + ", ".join(header))

    for required in [args.well_col, args.cq_col]:
        if required and not has_col_ci(header, required):
            raise RuntimeError(f"Column '{required}' was not found in Cq file. Available columns: {header}")
    for optional in [args.target_col, args.sample_col, args.condition_col, args.biological_rep_col, getattr(args, "analysis_col", None)]:
        if optional and not has_col_ci(header, optional):
            print(f"WARNING: optional column '{optional}' was not found in Cq file; corresponding output values will be blank.", file=sys.stderr)

    eff_by_well: Dict[str, Dict[str, Any]] = {normalize_well(r.get("Well", "")): r for r in eff_rows}
    merged: List[Dict[str, Any]] = []

    for row in cq_rows:
        well = normalize_well(get_col_ci(row, args.well_col))
        target = str(get_col_ci(row, args.target_col, "")).strip() if args.target_col else ""
        if (not target) and args.target_fallback_col:
            target = str(get_col_ci(row, args.target_fallback_col, "")).strip()
        sample = str(get_col_ci(row, args.sample_col, "")).strip() if args.sample_col else ""
        condition_raw = str(get_col_ci(row, args.condition_col, sample)).strip() if args.condition_col else sample
        condition = apply_condition_regex(condition_raw, getattr(args, "condition_regex", None), int(getattr(args, "condition_regex_group", 1)))
        biorep = str(get_col_ci(row, args.biological_rep_col, "")).strip() if args.biological_rep_col else ""
        analysis = str(get_col_ci(row, getattr(args, "analysis_col", None), "")).strip() if getattr(args, "analysis_col", None) else ""
        content = str(get_col_ci(row, args.content_col, "")).strip() if args.content_col else ""
        fluor = str(get_col_ci(row, args.fluor_col, "")).strip() if args.fluor_col else ""
        cq = to_float(get_col_ci(row, args.cq_col, ""))
        if cq is not None and cq <= 0:
            cq = None

        eff = eff_by_well.get(well, {})
        e_amp = to_float(eff.get("E_amp_factor"))
        e_pct = to_float(eff.get("Efficiency_percent"))
        eff_status = str(eff.get("status", ""))
        eff_self = str(eff.get("self_check_status", ""))
        eff_slope = to_float(eff.get("slope"))
        eff_intercept = to_float(eff.get("intercept"))
        eff_cycle_start = to_float(eff.get("cycle_start"))
        eff_cycle_end = to_float(eff.get("cycle_end"))
        eff_window_cycles = str(eff.get("window_cycles", ""))
        amp_sheet = str(eff.get("Amp_Sheet", ""))

        merged.append({
            "Well": well,
            "Amp_Sheet": amp_sheet,
            "Fluor": fluor,
            "Target": target,
            "Sample": sample,
            "Condition": condition,
            "BiologicalRep": biorep,
            "Analysis": analysis,
            "Content": content,
            "Cq": cq,
            "Cq_missing": cq is None,
            "E_amp_factor_per_well": e_amp,
            "Efficiency_percent_per_well": e_pct,
            "Efficiency_status": eff_status,
            "Efficiency_self_check": eff_self,
            "Eff_slope": eff_slope,
            "Eff_intercept": eff_intercept,
            "Eff_cycle_start": eff_cycle_start,
            "Eff_cycle_end": eff_cycle_end,
            "Eff_window_cycles": eff_window_cycles,
        })
    return merged, header


def _format_join_value(value: Any) -> str:
    """Stable CSV-friendly representation for per-well annotation values."""
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.6g}"
    return str(value).strip()


def _join_unique_values(rows: List[Dict[str, Any]], field: str) -> str:
    """Join unique non-empty values while preserving first-seen order."""
    out: List[str] = []
    seen = set()
    for row in rows:
        val = _format_join_value(row.get(field))
        if val and val not in seen:
            out.append(val)
            seen.add(val)
    return ";".join(out)


def annotate_efficiency_rows_from_cq(eff_rows: List[Dict[str, Any]], merged_rows: List[Dict[str, Any]]) -> None:
    """
    Add Cq-layout annotations to the per-well efficiency table.

    The raw amplification file only knows wells. The Cq/layout file knows what each
    well contains (target, sample, condition, biological replicate, etc.). This
    function copies those annotations back to *_per_well_efficiency.csv so BAD_FIT
    and low/high efficiency wells can be inspected by gene and sample instead of
    by well coordinate only.
    """
    rows_by_dataset_well: Dict[Tuple[str, str], List[Dict[str, Any]]] = {}
    rows_by_well: Dict[str, List[Dict[str, Any]]] = {}

    for row in merged_rows:
        well = normalize_well(row.get("Well", ""))
        if not well:
            continue
        dataset = str(row.get("Dataset", ""))
        rows_by_dataset_well.setdefault((dataset, well), []).append(row)
        rows_by_well.setdefault(well, []).append(row)

    fields_to_copy = [
        "Target", "Target_original", "Sample", "Condition", "BiologicalRep",
        "Analysis", "Content", "Fluor", "Cq", "Cq_source", "Cq_missing",
    ]

    for eff in eff_rows:
        well = normalize_well(eff.get("Well", ""))
        dataset = str(eff.get("Dataset", ""))
        matches = rows_by_dataset_well.get((dataset, well)) or rows_by_well.get(well, [])
        eff["Cq_annotation_found"] = bool(matches)
        for field in fields_to_copy:
            eff[field] = _join_unique_values(matches, field) if matches else ""



def normalize_key(s: Any) -> str:
    """Case-insensitive key for target/condition matching."""
    return str(s).strip().lower()


def apply_condition_regex(value: str, pattern: Optional[str], group: int = 1) -> str:
    """Optionally extract condition from Sample/Condition text, e.g. 'tfeb 10-1' -> '10'."""
    value = str(value).strip()
    if not pattern:
        return value
    try:
        m = re.search(pattern, value)
        if not m:
            return value
        return str(m.group(group)).strip()
    except Exception:
        return value


def calculate_missing_cq_from_amp_fit(rows: List[Dict[str, Any]], enabled: bool = True) -> None:
    """
    Fill missing Cq values using the same log-linear regression fit used for efficiency.

    Method:
      1) For each OK well, estimate the log10(corrected fluorescence) at the midpoint
         of the selected log-linear window: y_mid = slope * mid_cycle + intercept.
      2) For each target, use the median y_mid as an automatic target-specific threshold.
      3) For any row with missing Cq and a usable fit, calculate:
             Cq = (target_log10_threshold - intercept) / slope

    This is not a Bio-Rad Cq. It is a reproducible internally calculated Cq from
    the same baseline-corrected amplification curves. For relative expression,
    a target-specific constant threshold is acceptable because target-specific
    threshold constants cancel during control normalization, provided the same
    method is applied consistently across all samples for that target.
    """
    for r in rows:
        r["Cq_original"] = r.get("Cq")
        if r.get("Cq") is not None:
            r["Cq_source"] = "file"
        else:
            r["Cq_source"] = "missing"
        r["Cq_auto_log10_threshold"] = None
        r["Cq_auto_method"] = None

    if not enabled:
        return

    target_to_log_mids: Dict[str, List[float]] = {}
    target_original_name: Dict[str, str] = {}

    for r in rows:
        target = str(r.get("Target", "")).strip()
        if not target:
            continue
        slope = to_float(r.get("Eff_slope"))
        intercept = to_float(r.get("Eff_intercept"))
        c0 = to_float(r.get("Eff_cycle_start"))
        c1 = to_float(r.get("Eff_cycle_end"))
        status = str(r.get("Efficiency_status", ""))
        self_check = str(r.get("Efficiency_self_check", ""))
        if slope is None or intercept is None or c0 is None or c1 is None:
            continue
        if slope <= 0:
            continue
        if status != "OK":
            continue
        if self_check not in {"PASS", "", "NA"}:
            continue
        mid_cycle = (c0 + c1) / 2.0
        log_mid = slope * mid_cycle + intercept
        if math.isfinite(log_mid):
            key = normalize_key(target)
            target_to_log_mids.setdefault(key, []).append(log_mid)
            target_original_name.setdefault(key, target)

    target_to_threshold = {k: median(v) for k, v in target_to_log_mids.items() if v}

    for r in rows:
        if r.get("Cq") is not None:
            continue
        target = str(r.get("Target", "")).strip()
        key = normalize_key(target)
        threshold = target_to_threshold.get(key)
        slope = to_float(r.get("Eff_slope"))
        intercept = to_float(r.get("Eff_intercept"))
        if threshold is None or slope is None or intercept is None or slope <= 0:
            continue
        cq_calc = (threshold - intercept) / slope
        if math.isfinite(cq_calc) and cq_calc > 0:
            r["Cq"] = cq_calc
            r["Cq_missing"] = False
            r["Cq_source"] = "calculated_from_amp_curve"
            r["Cq_auto_log10_threshold"] = threshold
            r["Cq_auto_method"] = "target_median_selected_window_midpoint_log10"


def resolve_actual_target_names(rows: List[Dict[str, Any]], reference: Optional[str], targets: Optional[List[str]]) -> Tuple[Optional[str], Optional[List[str]]]:
    """Resolve reference/target names case-insensitively against actual Target values in the file."""
    actuals = []
    seen = set()
    for r in rows:
        t = str(r.get("Target", "")).strip()
        if t and normalize_key(t) not in seen:
            actuals.append(t)
            seen.add(normalize_key(t))
    key_to_actual = {normalize_key(t): t for t in actuals}

    ref_actual = None
    if reference:
        ref_actual = key_to_actual.get(normalize_key(reference), reference)

    targets_actual = None
    if targets:
        targets_actual = [key_to_actual.get(normalize_key(t), t) for t in targets]
    return ref_actual, targets_actual

def percentile(vals: List[float], p: float) -> Optional[float]:
    vals = sorted([float(v) for v in vals if v is not None and math.isfinite(float(v))])
    if not vals:
        return None
    if len(vals) == 1:
        return vals[0]
    k = (len(vals) - 1) * (p / 100.0)
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return vals[int(k)]
    return vals[f] * (c - k) + vals[c] * (k - f)


def compute_target_efficiency_qc(args: argparse.Namespace, rows: List[Dict[str, Any]]) -> Tuple[Dict[str, Optional[float]], List[Dict[str, Any]]]:
    """
    RDML-like target efficiency:
      individual reaction E is used for QC/outlier filtering;
      final quantities use filtered mean E per target.
    """
    per_target: Dict[str, List[Dict[str, Any]]] = {}

    for r in rows:
        target = str(r.get("Target", "")).strip()
        if not target:
            continue
        e = r.get("E_amp_factor_per_well")
        # For rdml_like target efficiency, use only wells whose log-linear window
        # passed strict noise/plateau + E/R2 + slope-stability checks. Cq values from
        # BAD_FIT wells remain in expression analysis, but their unreliable E does not
        # contaminate the target mean efficiency.
        status = str(r.get("Efficiency_status", ""))
        self_check = str(r.get("Efficiency_self_check", ""))
        ok = (
            e is not None
            and math.isfinite(float(e))
            and float(e) > 0
            and status == "OK"
            and self_check in {"PASS", "", "NA"}
        )
        if ok:
            per_target.setdefault(target, []).append(r)

    target_mean: Dict[str, Optional[float]] = {}
    qc_rows: List[Dict[str, Any]] = []

    for target in sorted({str(r.get("Target", "")).strip() for r in rows if str(r.get("Target", "")).strip()}):
        candidates = per_target.get(target, [])
        es = [float(r["E_amp_factor_per_well"]) for r in candidates]
        n_raw = len(es)

        kept_rows = candidates[:]
        outlier_low = None
        outlier_high = None
        method = getattr(args, "efficiency_outlier_mode", "iqr")

        if n_raw >= 4 and method == "iqr":
            q1 = percentile(es, 25)
            q3 = percentile(es, 75)
            iqr = (q3 - q1) if q1 is not None and q3 is not None else None
            if iqr is not None:
                k = float(getattr(args, "efficiency_iqr_k", 1.5))
                outlier_low = q1 - k * iqr
                outlier_high = q3 + k * iqr
                kept_rows = [r for r in candidates if outlier_low <= float(r["E_amp_factor_per_well"]) <= outlier_high]
        elif n_raw >= 4 and method == "mad":
            med = median(es)
            m = mad(es, med)
            k = float(getattr(args, "efficiency_mad_k", 3.5))
            if m and math.isfinite(m) and m > 0:
                # approximate robust z threshold in raw units
                outlier_low = med - k * 1.4826 * m
                outlier_high = med + k * 1.4826 * m
                kept_rows = [r for r in candidates if outlier_low <= float(r["E_amp_factor_per_well"]) <= outlier_high]
        elif method == "none":
            kept_rows = candidates[:]

        kept_es = [float(r["E_amp_factor_per_well"]) for r in kept_rows]
        e_mean = mean(kept_es)
        target_mean[target] = e_mean

        all_rows_for_target = [r for r in rows if str(r.get("Target", "")).strip() == target]
        rejected_status = len(all_rows_for_target) - n_raw
        outlier_wells = []
        kept_wells = []
        for r in candidates:
            well = str(r.get("Well", ""))
            e = float(r["E_amp_factor_per_well"])
            if r in kept_rows:
                kept_wells.append(well)
            else:
                outlier_wells.append(well)

        e_sd = stdev(kept_es)
        qc_rows.append({
            "Target": target,
            "E_target_mean": e_mean,
            "Efficiency_percent_target_mean": ((e_mean - 1.0) * 100.0) if e_mean is not None else None,
            "n_efficiency_OK_raw": n_raw,
            "n_efficiency_kept_after_outlier_filter": len(kept_es),
            "n_efficiency_rejected_bad_curve": rejected_status,
            "n_efficiency_rejected_outlier": len(outlier_wells),
            "sd_E_kept": e_sd,
            "cv_E_percent_kept": (100.0 * e_sd / e_mean) if e_mean not in (None, 0) and e_sd is not None else None,
            "outlier_method": method,
            "outlier_low": outlier_low,
            "outlier_high": outlier_high,
            "kept_wells": ";".join(kept_wells),
            "outlier_wells": ";".join(outlier_wells),
        })

    return target_mean, qc_rows



def parse_fixed_target_efficiency_map(spec: Optional[str]) -> Dict[str, float]:
    """Parse --fixed-target-efficiency values.

    Accepted formats:
      ppia=1.542,atg9b=1.739,lamp1=1.735,p62=1.698
      ppia=54.2%,atg9b=73.9%   # percent is converted to E factor = 1 + percent/100

    Target matching is case-insensitive; the original target spelling from the
    Cq file is preserved in outputs.
    """
    if spec is None or str(spec).strip() == "":
        raise RuntimeError("--fixed-target-efficiency is required when --efficiency-mode fixed_target")

    out: Dict[str, float] = {}
    items = [x.strip() for x in str(spec).replace(";", ",").split(",") if x.strip()]
    if not items:
        raise RuntimeError("--fixed-target-efficiency is empty. Expected format: ppia=1.542,atg9b=1.739")

    for item in items:
        if "=" not in item:
            raise RuntimeError(f"Bad --fixed-target-efficiency item '{item}'. Expected target=E, e.g. ppia=1.542")
        key, raw_val = item.split("=", 1)
        key = key.strip()
        raw_val = raw_val.strip().replace(",", ".")
        if not key:
            raise RuntimeError(f"Bad --fixed-target-efficiency item '{item}': empty target name")
        try:
            if raw_val.endswith("%"):
                pct = float(raw_val[:-1].strip())
                e = 1.0 + pct / 100.0
            else:
                e = float(raw_val)
        except Exception:
            raise RuntimeError(f"Bad efficiency value for target '{key}': '{raw_val}'")
        if not math.isfinite(e) or e <= 1.0:
            raise RuntimeError(f"Efficiency for target '{key}' must be an amplification factor > 1.0, got: {e}")
        out[normalize_key(key)] = e
    return out


def assign_efficiency_by_mode(rows: List[Dict[str, Any]], mode: str, args: Optional[argparse.Namespace] = None) -> Tuple[Dict[str, Optional[float]], List[Dict[str, Any]]]:
    """
    Mutates rows: add E_used, E_source, Quantity_E_minus_Cq and efficiency QC fields.

    Modes:
      fixed: use one user-defined efficiency factor for all targets, e.g. E=2.0.
      fixed_target: use user-defined efficiency factors per target/gene.
      target_mean: one mean E per target/gene.
      rdml_like: same final behavior as target_mean, but optional outlier filtering can be used for target E.
      per_well: each reaction uses its own E, with target_mean fallback.
      sample_target_mean: one mean E per biological sample x target; this is the
        sample-specific alternative for possible inhibitors/impurities.
    """
    if args is None:
        args = argparse.Namespace(efficiency_outlier_mode="none", efficiency_iqr_k=1.5, efficiency_mad_k=3.5)

    # Always compute global target means for QC and fallback.
    if mode in {"fixed", "fixed_target"}:
        targets_in_rows = sorted({str(r.get("Target", "")).strip() for r in rows if str(r.get("Target", "")).strip()})

        if mode == "fixed":
            fixed_e = float(getattr(args, "fixed_efficiency", 2.0))
            if not math.isfinite(fixed_e) or fixed_e <= 1.0:
                raise RuntimeError(f"--fixed-efficiency must be a finite amplification factor > 1.0, got: {fixed_e}")
            target_mean = {target: fixed_e for target in targets_in_rows}
            fixed_label = f"fixed_efficiency_{fixed_e:g}"
        else:
            fixed_map = parse_fixed_target_efficiency_map(getattr(args, "fixed_target_efficiency", None))
            missing = [target for target in targets_in_rows if normalize_key(target) not in fixed_map]
            if missing:
                raise RuntimeError(
                    "--efficiency-mode fixed_target requires an efficiency for every Target in the data. "
                    f"Missing: {', '.join(missing)}. Provided keys: {', '.join(sorted(fixed_map))}"
                )
            target_mean = {target: fixed_map[normalize_key(target)] for target in targets_in_rows}
            fixed_label = "fixed_target_efficiency"

        qc_rows = []
        for target in targets_in_rows:
            target_rows = [r for r in rows if str(r.get("Target", "")).strip() == target]
            ok_rows = [r for r in target_rows if str(r.get("Efficiency_status", "")) == "OK" and str(r.get("Efficiency_self_check", "")) in {"PASS", "", "NA"}]
            actual_es = [float(r.get("E_amp_factor_per_well")) for r in ok_rows if r.get("E_amp_factor_per_well") is not None and math.isfinite(float(r.get("E_amp_factor_per_well")))]
            actual_sd = stdev(actual_es)
            actual_mean = mean(actual_es)
            fixed_e_for_target = target_mean.get(target)
            qc_rows.append({
                "Target": target,
                "E_target_mean": fixed_e_for_target,
                "Efficiency_percent_target_mean": ((fixed_e_for_target - 1.0) * 100.0) if fixed_e_for_target is not None else None,
                "n_efficiency_OK_raw": len(ok_rows),
                "n_efficiency_kept_after_outlier_filter": len(ok_rows),
                "n_efficiency_rejected_bad_curve": len(target_rows) - len(ok_rows),
                "n_efficiency_rejected_outlier": 0,
                "sd_E_kept": actual_sd,
                "cv_E_percent_kept": (100.0 * actual_sd / actual_mean) if actual_mean not in (None, 0) and actual_sd is not None else None,
                "outlier_method": f"{fixed_label};actual_curve_E_kept_for_QC_only",
                "outlier_low": None,
                "outlier_high": None,
                "kept_wells": ";".join(str(r.get("Well", "")) for r in ok_rows),
                "outlier_wells": "",
            })
    elif mode == "rdml_like":
        target_mean, qc_rows = compute_target_efficiency_qc(args, rows)
    else:
        target_to_es: Dict[str, List[float]] = {}
        for r in rows:
            target = str(r.get("Target", "")).strip()
            e = r.get("E_amp_factor_per_well")
            if target and e is not None and math.isfinite(float(e)) and float(e) > 0:
                target_to_es.setdefault(target, []).append(float(e))
        target_mean = {t: mean(vs) for t, vs in target_to_es.items()}
        qc_rows = []
        for target in sorted({str(r.get("Target", "")).strip() for r in rows if str(r.get("Target", "")).strip()}):
            vals = target_to_es.get(target, [])
            em = target_mean.get(target)
            sd = stdev(vals)
            qc_rows.append({
                "Target": target,
                "E_target_mean": em,
                "Efficiency_percent_target_mean": ((em - 1.0) * 100.0) if em is not None else None,
                "n_efficiency_OK_raw": len(vals),
                "n_efficiency_kept_after_outlier_filter": len(vals),
                "n_efficiency_rejected_bad_curve": len([r for r in rows if str(r.get("Target", "")).strip() == target]) - len(vals),
                "n_efficiency_rejected_outlier": 0,
                "sd_E_kept": sd,
                "cv_E_percent_kept": (100.0 * sd / em) if em not in (None, 0) and sd is not None else None,
                "outlier_method": "none_for_" + mode,
                "outlier_low": None,
                "outlier_high": None,
                "kept_wells": ";".join(str(r.get("Well", "")) for r in rows if str(r.get("Target", "")).strip() == target and r.get("E_amp_factor_per_well") is not None),
                "outlier_wells": "",
            })

    # Optional sample-specific means: one E per Condition/BiologicalRep/Target.
    sample_target_es: Dict[Tuple[str, str, str], List[float]] = {}
    for r in rows:
        target = str(r.get("Target", "")).strip()
        condition = str(r.get("Condition", "")).strip()
        biorep = str(r.get("BiologicalRep", "")).strip()
        e = r.get("E_amp_factor_per_well")
        status = str(r.get("Efficiency_status", ""))
        self_check = str(r.get("Efficiency_self_check", ""))
        if (
            target and condition and e is not None and math.isfinite(float(e)) and float(e) > 0
            and status == "OK" and self_check in {"PASS", "", "NA"}
        ):
            sample_target_es.setdefault((condition, biorep, target), []).append(float(e))

    sample_target_stats: Dict[Tuple[str, str, str], Dict[str, Any]] = {}
    for key, vals in sample_target_es.items():
        em = mean(vals)
        esd = stdev(vals)
        target = key[2]
        global_e = target_mean.get(target)
        ratio = (em / global_e) if em is not None and global_e not in (None, 0) else None
        delta = (em - global_e) if em is not None and global_e is not None else None
        cv = (100.0 * esd / em) if em not in (None, 0) and esd is not None else None
        sample_target_stats[key] = {
            "E_sample_target_mean": em,
            "Efficiency_percent_sample_target_mean": ((em - 1.0) * 100.0) if em is not None else None,
            "E_sample_target_sd": esd,
            "E_sample_target_cv_percent": cv,
            "n_E_sample_target": len(vals),
            "E_target_global_mean": global_e,
            "E_sample_vs_target_mean_ratio": ratio,
            "E_sample_vs_target_mean_delta": delta,
        }

    for r in rows:
        target = str(r.get("Target", "")).strip()
        condition = str(r.get("Condition", "")).strip()
        biorep = str(r.get("BiologicalRep", "")).strip()
        key = (condition, biorep, target)
        st = sample_target_stats.get(key, {})

        # Store both global and sample-specific efficiency information for QC/export.
        r["E_target_global_mean"] = target_mean.get(target)
        r["Efficiency_percent_target_global_mean"] = ((target_mean[target] - 1.0) * 100.0) if target in target_mean and target_mean.get(target) is not None else None
        for fld in [
            "E_sample_target_mean", "Efficiency_percent_sample_target_mean", "E_sample_target_sd",
            "E_sample_target_cv_percent", "n_E_sample_target", "E_sample_vs_target_mean_ratio",
            "E_sample_vs_target_mean_delta",
        ]:
            r[fld] = st.get(fld)

        if mode == "fixed":
            r["E_used"] = target_mean.get(target)
            r["E_source"] = f"fixed_efficiency_{float(getattr(args, 'fixed_efficiency', 2.0)):g}"
        elif mode == "fixed_target":
            r["E_used"] = target_mean.get(target)
            r["E_source"] = f"fixed_target_efficiency_{target}={target_mean.get(target):g}" if target_mean.get(target) is not None else "fixed_target_efficiency_missing"
        elif mode == "per_well":
            e = r.get("E_amp_factor_per_well")
            if e is not None and math.isfinite(float(e)) and float(e) > 0:
                r["E_used"] = float(e)
                r["E_source"] = "per_well"
            else:
                r["E_used"] = target_mean.get(target)
                r["E_source"] = "target_mean_fallback"
        elif mode == "sample_target_mean":
            e_sample = st.get("E_sample_target_mean")
            if e_sample is not None and math.isfinite(float(e_sample)) and float(e_sample) > 0:
                r["E_used"] = float(e_sample)
                r["E_source"] = "sample_target_mean"
            else:
                r["E_used"] = target_mean.get(target)
                r["E_source"] = "target_mean_fallback_no_sample_E"
        elif mode in {"target_mean", "rdml_like"}:
            r["E_used"] = target_mean.get(target)
            r["E_source"] = "target_mean_filtered_rdml_like" if mode == "rdml_like" else "target_mean"
        else:
            raise RuntimeError(f"Unknown efficiency mode: {mode}")

        # QC warnings about sample-specific E, regardless of mode.
        e_used = r.get("E_used")
        qc_parts = []
        existing_qc = str(r.get("Efficiency_status", "")).strip()
        if existing_qc and existing_qc != "OK":
            qc_parts.append(existing_qc)
        cv = r.get("E_sample_target_cv_percent")
        ratio = r.get("E_sample_vs_target_mean_ratio")
        if cv is not None and math.isfinite(float(cv)) and float(cv) > 10.0:
            qc_parts.append("sample_target_E_CV_gt_10pct")
        if ratio is not None and math.isfinite(float(ratio)) and abs(float(ratio) - 1.0) > 0.10:
            qc_parts.append("sample_target_E_differs_from_global_gt_10pct")
        if mode == "sample_target_mean" and str(r.get("E_source")) == "target_mean_fallback_no_sample_E":
            qc_parts.append("no_sample_target_E_used_global_fallback")
        r["Efficiency_QC_warning"] = ";".join(sorted(set(qc_parts)))

        cq = r.get("Cq")
        e_used = r.get("E_used")
        if cq is not None and e_used is not None and e_used > 0:
            try:
                r["Quantity_E_minus_Cq"] = float(e_used) ** (-float(cq))
            except Exception:
                r["Quantity_E_minus_Cq"] = None
        else:
            r["Quantity_E_minus_Cq"] = None

    return target_mean, qc_rows


def safe_name_for_file(value: Any) -> str:
    """Make a short filesystem-safe name for analysis-specific output files."""
    s = str(value).strip()
    s = re.sub(r"[^A-Za-z0-9А-Яа-я_.-]+", "_", s).strip("_")
    return s or "analysis"


def sniff_delimiter(path: str) -> str:
    """Detect delimiter for small CSV/TSV mapping files."""
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        sample = f.read(4096)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        return dialect.delimiter
    except Exception:
        if "\t" in sample:
            return "\t"
        if ";" in sample:
            return ";"
        return ","


def read_mapping_table(path: str) -> List[Dict[str, str]]:
    """Read a CSV/TSV mapping table into case-preserving string records."""
    delim = sniff_delimiter(path)
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        if not reader.fieldnames:
            raise RuntimeError(f"Mapping file has no header: {path}")
        rows: List[Dict[str, str]] = []
        for rec in reader:
            clean = {str(k).strip(): str(v).strip() for k, v in rec.items() if k is not None}
            if any(v for v in clean.values()):
                rows.append(clean)
    return rows


def get_mapping_value_ci(row: Dict[str, str], names: List[str], default: str = "") -> str:
    """Case-insensitive lookup in user mapping files."""
    lowered = {k.strip().lower(): v for k, v in row.items()}
    for name in names:
        val = lowered.get(name.strip().lower())
        if val is not None:
            return str(val).strip()
    return default


def split_list_cell(value: str) -> List[str]:
    """Split cells like 'A; B, C' while preserving single ordinary condition names."""
    value = str(value).strip()
    if not value:
        return []
    parts = re.split(r"\s*[;,]\s*", value)
    return [p.strip() for p in parts if p.strip()]


def parse_control_map(value: Optional[str]) -> Dict[str, str]:
    """Parse --control-map like 'exp1=k-;exp2=mock'."""
    out: Dict[str, str] = {}
    if not value:
        return out
    for item in re.split(r"\s*;\s*", str(value).strip()):
        if not item:
            continue
        if "=" not in item:
            raise RuntimeError("Bad --control-map item. Use format: analysis1=control1;analysis2=control2")
        a, c = item.split("=", 1)
        a = a.strip()
        c = c.strip()
        if a and c:
            out[a] = c
    return out


def build_analysis_plan(args: argparse.Namespace, rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Build independent analysis groups for relative-expression normalization.

    Default behavior stays backward-compatible: one analysis containing all conditions,
    normalized to --control-condition.

    Optional modes:
      1) --analysis-map CSV/TSV with columns:
         Analysis, Condition, ControlCondition
         One row = one condition included in one analysis. The same Condition may be
         repeated in several analyses if needed.

      2) --analysis-col in the Cq/layout file plus either --control-map or
         --control-condition. This is convenient when every row already has an
         experiment/batch/group column.
    """
    available_conditions = sorted({str(r.get("Condition", "")).strip() for r in rows if str(r.get("Condition", "")).strip()})
    condition_by_key = {normalize_key(c): c for c in available_conditions}

    plans_by_name: Dict[str, Dict[str, Any]] = {}

    if getattr(args, "analysis_map", None):
        if not os.path.exists(args.analysis_map):
            raise RuntimeError(f"Analysis map file not found: {args.analysis_map}")
        map_rows = read_mapping_table(args.analysis_map)
        if not map_rows:
            raise RuntimeError(f"Analysis map is empty: {args.analysis_map}")

        for i, rec in enumerate(map_rows, start=2):
            analysis = get_mapping_value_ci(rec, ["Analysis", "Experiment", "Group", "AnalysisGroup", "RunGroup"])
            condition_cell = get_mapping_value_ci(rec, ["Condition", "Sample", "SampleCondition"])
            control = get_mapping_value_ci(rec, ["ControlCondition", "Control", "Control_Condition", "Calibrator", "ReferenceCondition"])
            include = get_mapping_value_ci(rec, ["Include", "Use"], "1").lower()
            if include in {"0", "no", "false", "exclude", "skip"}:
                continue
            if not analysis:
                analysis = "analysis_1"
            if not condition_cell:
                raise RuntimeError(f"Analysis map row {i}: missing Condition")
            if not control:
                raise RuntimeError(f"Analysis map row {i}: missing ControlCondition")

            plan = plans_by_name.setdefault(analysis, {"Analysis": analysis, "control_condition": control, "conditions": []})
            if normalize_key(plan["control_condition"]) != normalize_key(control):
                raise RuntimeError(
                    f"Analysis map row {i}: analysis '{analysis}' has conflicting controls: "
                    f"'{plan['control_condition']}' and '{control}'"
                )
            for condition in split_list_cell(condition_cell):
                actual = condition_by_key.get(normalize_key(condition), condition)
                if actual not in plan["conditions"]:
                    plan["conditions"].append(actual)

    elif getattr(args, "analysis_col", None):
        control_map = parse_control_map(getattr(args, "control_map", None))
        grouped_conditions: Dict[str, List[str]] = {}
        for r in rows:
            analysis = str(r.get("Analysis", "")).strip() or "analysis_1"
            condition = str(r.get("Condition", "")).strip()
            if condition:
                grouped_conditions.setdefault(analysis, [])
                if condition not in grouped_conditions[analysis]:
                    grouped_conditions[analysis].append(condition)
        for analysis, conditions in sorted(grouped_conditions.items()):
            control = control_map.get(analysis) or getattr(args, "control_condition", None)
            if not control:
                raise RuntimeError(
                    f"No control specified for analysis '{analysis}'. Use --control-map '{analysis}=control' "
                    f"or provide --control-condition."
                )
            plans_by_name[analysis] = {"Analysis": analysis, "control_condition": control, "conditions": conditions}

    else:
        control = getattr(args, "control_condition", None)
        if not control:
            raise RuntimeError("--control-condition is required unless --analysis-map provides ControlCondition values")
        plans_by_name["analysis_1"] = {"Analysis": "analysis_1", "control_condition": control, "conditions": available_conditions}

    plans: List[Dict[str, Any]] = []
    for analysis, plan in plans_by_name.items():
        conditions = plan.get("conditions", [])
        control = str(plan.get("control_condition", "")).strip()
        actual_control = condition_by_key.get(normalize_key(control), control)
        if actual_control not in conditions:
            conditions = [actual_control] + [c for c in conditions if normalize_key(c) != normalize_key(actual_control)]
        missing = [c for c in conditions if normalize_key(c) not in condition_by_key]
        if missing:
            print(
                f"WARNING: analysis '{analysis}' contains condition(s) not found in Cq table: {', '.join(missing)}",
                file=sys.stderr,
            )
        present_conditions = [condition_by_key.get(normalize_key(c), c) for c in conditions if normalize_key(c) in condition_by_key]
        if not present_conditions:
            print(f"WARNING: analysis '{analysis}' has no present conditions; skipped.", file=sys.stderr)
            continue
        if normalize_key(actual_control) not in {normalize_key(c) for c in present_conditions}:
            print(
                f"WARNING: control '{control}' for analysis '{analysis}' is not present in the selected conditions; "
                f"relative expression will be empty for that analysis.",
                file=sys.stderr,
            )
        plans.append({"Analysis": analysis, "control_condition": actual_control, "conditions": present_conditions})

    if not plans:
        raise RuntimeError("No valid analysis groups were produced")
    return plans

def aggregate_relative_expression(args: argparse.Namespace, rows: List[Dict[str, Any]], target_mean_map: Optional[Dict[str, Optional[float]]] = None) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    reference, resolved_targets = resolve_actual_target_names(rows, args.reference, args.targets)
    if not reference:
        return [], []

    all_targets = sorted({str(r.get("Target", "")).strip() for r in rows if str(r.get("Target", "")).strip()})
    targets = resolved_targets if resolved_targets else [t for t in all_targets if normalize_key(t) != normalize_key(reference)]
    target_mean_map = target_mean_map or {}

    analysis_plans = build_analysis_plan(args, rows)
    print("Analysis groups:")
    for plan in analysis_plans:
        print(
            f"  - {plan['Analysis']}: control={plan['control_condition']}; "
            f"conditions={', '.join(plan['conditions'])}"
        )

    # Technical replicate quantities by condition + biological replicate + target.
    # This basic aggregation is shared across all analysis groups; each group then
    # chooses its own included conditions and control mean.
    grouped: Dict[Tuple[str, str, str], List[float]] = {}
    cq_grouped: Dict[Tuple[str, str, str], List[float]] = {}
    wells_grouped: Dict[Tuple[str, str, str], List[str]] = {}
    status_grouped: Dict[Tuple[str, str, str], List[str]] = {}
    e_used_grouped: Dict[Tuple[str, str, str], List[float]] = {}
    e_sample_grouped: Dict[Tuple[str, str, str], List[float]] = {}
    e_global_grouped: Dict[Tuple[str, str, str], List[float]] = {}
    e_source_grouped: Dict[Tuple[str, str, str], List[str]] = {}
    eff_qc_grouped: Dict[Tuple[str, str, str], List[str]] = {}
    for r in rows:
        target = str(r.get("Target", "")).strip()
        if not target:
            continue
        q = r.get("Quantity_E_minus_Cq")
        cq = r.get("Cq")
        if q is None or cq is None:
            continue
        condition = str(r.get("Condition", "")).strip()
        biorep = str(r.get("BiologicalRep", "")).strip()
        key = (condition, biorep, target)
        grouped.setdefault(key, []).append(float(q))
        cq_grouped.setdefault(key, []).append(float(cq))
        wells_grouped.setdefault(key, []).append(str(r.get("Well", "")))
        status_grouped.setdefault(key, []).append(str(r.get("Efficiency_status", "")))
        eu = r.get("E_used")
        if eu is not None and math.isfinite(float(eu)):
            e_used_grouped.setdefault(key, []).append(float(eu))
        es = r.get("E_sample_target_mean")
        if es is not None and math.isfinite(float(es)):
            e_sample_grouped.setdefault(key, []).append(float(es))
        eg = r.get("E_target_global_mean")
        if eg is not None and math.isfinite(float(eg)):
            e_global_grouped.setdefault(key, []).append(float(eg))
        src = str(r.get("E_source", "")).strip()
        if src:
            e_source_grouped.setdefault(key, []).append(src)
        eqc = str(r.get("Efficiency_QC_warning", "")).strip()
        if eqc:
            eff_qc_grouped.setdefault(key, []).append(eqc)

    agg: Dict[Tuple[str, str, str], Dict[str, Any]] = {}
    for key, qs in grouped.items():
        c, b, t = key
        cq_vals = cq_grouped[key]
        e_used_vals = e_used_grouped.get(key, [])
        e_sample_vals = e_sample_grouped.get(key, [])
        e_global_vals = e_global_grouped.get(key, [])
        agg[key] = {
            "Condition": c,
            "BiologicalRep": b,
            "Unit": f"{c}|{b}" if b else c,
            "Target": t,
            "n_technical": len(qs),
            "mean_quantity_E_minus_Cq": mean(qs),
            "sd_quantity_E_minus_Cq": stdev(qs),
            "mean_Cq": mean(cq_vals),
            "sd_Cq": stdev(cq_vals),
            "wells": ";".join(wells_grouped[key]),
            "efficiency_statuses": ";".join(sorted(set(status_grouped[key]))),
            "E_used_mean": mean(e_used_vals),
            "E_used_sd": stdev(e_used_vals),
            "E_sample_target_mean": mean(e_sample_vals),
            "E_sample_target_sd": stdev(e_sample_vals),
            "E_sample_target_cv_percent": (100.0 * stdev(e_sample_vals) / mean(e_sample_vals)) if mean(e_sample_vals) not in (None, 0) and stdev(e_sample_vals) is not None else None,
            "E_target_global_mean": mean(e_global_vals),
            "E_source": ";".join(sorted(set(e_source_grouped.get(key, [])))),
            "Efficiency_QC_warning": ";".join(sorted(set(w for item in eff_qc_grouped.get(key, []) for w in str(item).split(';') if w))),
        }

    rel_rows: List[Dict[str, Any]] = []
    for plan in analysis_plans:
        analysis = str(plan["Analysis"])
        control = str(plan["control_condition"])
        selected_condition_keys = {normalize_key(c) for c in plan["conditions"]}
        selected_conditions = [c for c in plan["conditions"] if normalize_key(c) in selected_condition_keys]

        for condition in selected_conditions:
            bioreps = sorted({k[1] for k in agg if normalize_key(k[0]) == normalize_key(condition)})
            for biorep in bioreps:
                # Use the actual condition name from agg in case only case differs.
                matching_conditions = sorted({k[0] for k in agg if normalize_key(k[0]) == normalize_key(condition)})
                if not matching_conditions:
                    continue
                actual_condition = matching_conditions[0]
                ref_key = (actual_condition, biorep, reference)
                ref_agg = agg.get(ref_key)
                if not ref_agg:
                    continue
                ref_q = ref_agg.get("mean_quantity_E_minus_Cq")
                if ref_q is None or ref_q == 0:
                    continue
                for target in targets:
                    key = (actual_condition, biorep, target)
                    ta = agg.get(key)
                    if not ta:
                        continue
                    tq = ta.get("mean_quantity_E_minus_Cq")
                    if tq is None:
                        continue
                    normalized = tq / ref_q
                    warning_parts = []
                    if ta["n_technical"] < 2:
                        warning_parts.append("target_n_tech<2")
                    if ref_agg["n_technical"] < 2:
                        warning_parts.append("reference_n_tech<2")
                    if "OK" not in ta["efficiency_statuses"]:
                        warning_parts.append("target_efficiency_not_OK")
                    if "OK" not in ref_agg["efficiency_statuses"]:
                        warning_parts.append("reference_efficiency_not_OK")
                    for w in str(ta.get("Efficiency_QC_warning", "")).split(";"):
                        if w:
                            warning_parts.append("target_" + w)
                    for w in str(ref_agg.get("Efficiency_QC_warning", "")).split(";"):
                        if w:
                            warning_parts.append("reference_" + w)

                    rel_rows.append({
                        "Analysis": analysis,
                        "Target": target,
                        "Condition": actual_condition,
                        "BiologicalRep": biorep,
                        "Unit": ta["Unit"],
                        "Reference": reference,
                        "control_condition": control,
                        "Efficiency_mode": getattr(args, "efficiency_mode", ""),
                        "relative_expression": None,  # filled below
                        "log2_relative_expression": None,
                        "normalized_quantity_target_over_ref": normalized,
                        "control_mean_normalized_quantity": None,
                        "mean_Cq_target": ta["mean_Cq"],
                        "sd_Cq_target": ta["sd_Cq"],
                        "mean_Cq_reference": ref_agg["mean_Cq"],
                        "sd_Cq_reference": ref_agg["sd_Cq"],
                        "E_target_used": ta.get("E_used_mean"),
                        "E_reference_used": ref_agg.get("E_used_mean"),
                        "E_target_global_mean": ta.get("E_target_global_mean") or target_mean_map.get(target),
                        "E_reference_global_mean": ref_agg.get("E_target_global_mean") or target_mean_map.get(reference),
                        "E_sample_target_mean": ta.get("E_sample_target_mean"),
                        "E_sample_target_sd": ta.get("E_sample_target_sd"),
                        "E_sample_target_cv_percent": ta.get("E_sample_target_cv_percent"),
                        "E_sample_reference_mean": ref_agg.get("E_sample_target_mean"),
                        "E_sample_reference_sd": ref_agg.get("E_sample_target_sd"),
                        "E_sample_reference_cv_percent": ref_agg.get("E_sample_target_cv_percent"),
                        "E_source_target": ta.get("E_source"),
                        "E_source_reference": ref_agg.get("E_source"),
                        "n_tech_target": ta["n_technical"],
                        "n_tech_reference": ref_agg["n_technical"],
                        "target_wells": ta["wells"],
                        "reference_wells": ref_agg["wells"],
                        "target_efficiency_statuses": ta["efficiency_statuses"],
                        "reference_efficiency_statuses": ref_agg["efficiency_statuses"],
                        "QC_warning": ";".join(warning_parts) if warning_parts else "",
                    })

        # Normalize each target independently to the control inside this analysis only.
        for target in targets:
            control_vals = [
                r["normalized_quantity_target_over_ref"]
                for r in rel_rows
                if r["Analysis"] == analysis
                and normalize_key(r["Target"]) == normalize_key(target)
                and normalize_key(r["Condition"]) == normalize_key(control)
            ]
            control_mean = mean(control_vals)
            for r in rel_rows:
                if r["Analysis"] == analysis and normalize_key(r["Target"]) == normalize_key(target):
                    r["control_mean_normalized_quantity"] = control_mean
                    if control_mean is not None and control_mean != 0:
                        r["relative_expression"] = r["normalized_quantity_target_over_ref"] / control_mean
                        r["log2_relative_expression"] = math.log(r["relative_expression"], 2) if r["relative_expression"] > 0 else None
                    else:
                        r["relative_expression"] = None
                        r["log2_relative_expression"] = None
                        r["QC_warning"] = (r.get("QC_warning", "") + ";no_control_mean").strip(";")

    # Clean summary by analysis, target and condition.
    summary: List[Dict[str, Any]] = []
    analyses = sorted({r["Analysis"] for r in rel_rows})
    for analysis in analyses:
        analysis_rows = [r for r in rel_rows if r["Analysis"] == analysis]
        for target in targets:
            conditions = sorted({r["Condition"] for r in analysis_rows if normalize_key(r["Target"]) == normalize_key(target)})
            for condition in conditions:
                sub = [r for r in analysis_rows if normalize_key(r["Target"]) == normalize_key(target) and r["Condition"] == condition]
                vals = [r.get("relative_expression") for r in sub]
                vals = [float(v) for v in vals if v is not None and math.isfinite(float(v))]
                logs = [r.get("log2_relative_expression") for r in sub]
                logs = [float(v) for v in logs if v is not None and math.isfinite(float(v))]
                sd = stdev(vals)
                sd_log = stdev(logs)
                warnings = sorted({w for r in sub for w in str(r.get("QC_warning", "")).split(";") if w})
                control = sub[0].get("control_condition", "") if sub else ""
                summary.append({
                    "Analysis": analysis,
                    "Target": target,
                    "Condition": condition,
                    "Reference": reference,
                    "control_condition": control,
                    "n_units": len(vals),
                    "mean_relative_expression": mean(vals),
                    "sd_relative_expression": sd,
                    "sem_relative_expression": (sd / math.sqrt(len(vals))) if sd is not None and len(vals) > 0 else None,
                    "mean_log2_relative_expression": mean(logs),
                    "sd_log2_relative_expression": sd_log,
                    "sem_log2_relative_expression": (sd_log / math.sqrt(len(logs))) if sd_log is not None and len(logs) > 0 else None,
                    "E_target_used": target_mean_map.get(target),
                    "Efficiency_percent_target_used": ((target_mean_map[target] - 1) * 100.0) if target in target_mean_map and target_mean_map.get(target) is not None else None,
                    "E_reference_used": target_mean_map.get(reference),
                    "Efficiency_percent_reference_used": ((target_mean_map[reference] - 1) * 100.0) if reference in target_mean_map and target_mean_map.get(reference) is not None else None,
                    "QC_warning": ";".join(warnings),
                })
    return rel_rows, summary



def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for incomplete beta; dependency-free p-values."""
    MAXIT = 200
    EPS = 3.0e-12
    FPMIN = 1.0e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        delh = d * c
        h *= delh
        if abs(delh - 1.0) < EPS:
            break
    return h


def _regularized_beta(a: float, b: float, x: float) -> Optional[float]:
    if not (a > 0 and b > 0) or x < 0 or x > 1:
        return None
    if x == 0:
        return 0.0
    if x == 1:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def _t_two_sided_p(t: float, df: float) -> Optional[float]:
    if not (math.isfinite(t) and math.isfinite(df) and df > 0):
        return None
    x = df / (df + t * t)
    ib = _regularized_beta(df / 2.0, 0.5, x)
    return max(0.0, min(1.0, ib)) if ib is not None else None


def _f_right_tail_p(F: float, df1: float, df2: float) -> Optional[float]:
    if not (math.isfinite(F) and F >= 0 and df1 > 0 and df2 > 0):
        return None
    x = (df1 * F) / (df1 * F + df2)
    cdf = _regularized_beta(df1 / 2.0, df2 / 2.0, x)
    if cdf is None:
        return None
    return max(0.0, min(1.0, 1.0 - cdf))


def welch_ttest_p(a: List[float], b: List[float]) -> Optional[float]:
    a = [float(x) for x in a if x is not None and math.isfinite(float(x))]
    b = [float(x) for x in b if x is not None and math.isfinite(float(x))]
    if len(a) < 2 or len(b) < 2:
        return None
    ma, mb = mean(a), mean(b)
    va = stdev(a)
    vb = stdev(b)
    if ma is None or mb is None or va is None or vb is None:
        return None
    va2 = va * va
    vb2 = vb * vb
    se2 = va2 / len(a) + vb2 / len(b)
    if se2 <= 0:
        return None
    t = (ma - mb) / math.sqrt(se2)
    df_num = se2 * se2
    df_den = (va2 / len(a)) ** 2 / (len(a) - 1) + (vb2 / len(b)) ** 2 / (len(b) - 1)
    if df_den <= 0:
        return None
    return _t_two_sided_p(t, df_num / df_den)


def one_way_anova_p(groups: List[List[float]]) -> Optional[float]:
    groups = [[float(x) for x in g if x is not None and math.isfinite(float(x))] for g in groups]
    groups = [g for g in groups if len(g) >= 1]
    if len(groups) < 2 or sum(len(g) for g in groups) <= len(groups):
        return None
    all_vals = [x for g in groups for x in g]
    grand = mean(all_vals)
    if grand is None:
        return None
    ss_between = sum(len(g) * ((mean(g) or 0.0) - grand) ** 2 for g in groups)
    ss_within = sum(sum((x - (mean(g) or 0.0)) ** 2 for x in g) for g in groups)
    df_between = len(groups) - 1
    df_within = len(all_vals) - len(groups)
    if df_between <= 0 or df_within <= 0 or ss_within <= 0:
        return None
    F = (ss_between / df_between) / (ss_within / df_within)
    return _f_right_tail_p(F, df_between, df_within)


def adjust_pvalues(pvals: List[Optional[float]], method: str = "bh") -> List[Optional[float]]:
    method = (method or "bh").lower()
    out: List[Optional[float]] = [None] * len(pvals)
    indexed = [(i, float(p)) for i, p in enumerate(pvals) if p is not None and math.isfinite(float(p))]
    m = len(indexed)
    if m == 0 or method == "none":
        for i, p in indexed:
            out[i] = p
        return out
    if method == "bonferroni":
        for i, p in indexed:
            out[i] = min(1.0, p * m)
        return out
    if method == "holm":
        sorted_pairs = sorted(indexed, key=lambda x: x[1])
        prev = 0.0
        for rank, (i, p) in enumerate(sorted_pairs, start=1):
            adj = min(1.0, (m - rank + 1) * p)
            adj = max(adj, prev)
            out[i] = adj
            prev = adj
        return out
    # Benjamini-Hochberg
    sorted_pairs = sorted(indexed, key=lambda x: x[1], reverse=True)
    prev = 1.0
    for rank_from_high, (i, p) in enumerate(sorted_pairs, start=1):
        rank = m - rank_from_high + 1
        adj = min(prev, p * m / rank)
        out[i] = min(1.0, adj)
        prev = adj
    return out


def p_to_stars(p: Optional[float]) -> str:
    if p is None:
        return "NA"
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def _clean_numeric(vals: List[Any]) -> List[float]:
    out: List[float] = []
    for v in vals:
        try:
            fv = float(v)
        except Exception:
            continue
        if math.isfinite(fv):
            out.append(fv)
    return out


def anova_status_from_groups(values_by_cond: Dict[str, List[float]], control: str = "") -> Tuple[str, str]:
    """Human-readable status for one-way ANOVA preconditions and result availability."""
    nonempty = {c: _clean_numeric(v) for c, v in values_by_cond.items() if len(_clean_numeric(v)) > 0}
    if len(nonempty) < 2:
        return "NOT_RUN", "fewer_than_2_groups_with_values"
    total_n = sum(len(v) for v in nonempty.values())
    if total_n <= len(nonempty):
        return "NOT_RUN", "no_within_group_df; need_at_least_one_group_with_n>=2"
    if all(len(v) < 2 for v in nonempty.values()):
        return "NOT_RUN", "all_groups_have_n=1; ANOVA_needs_within_group_variance"
    if all((stdev(v) in (None, 0)) for v in nonempty.values()):
        return "NOT_RUN", "zero_within_group_variance"
    return "OK", "ANOVA_on_log2_relative_expression"


def pairwise_status_for_groups(test_vals: List[float], control_vals: List[float], is_control: bool = False) -> Tuple[str, str]:
    """Human-readable status for Welch t-test versus control."""
    if is_control:
        return "CONTROL", "control_group_not_compared_to_itself"
    a = _clean_numeric(test_vals)
    b = _clean_numeric(control_vals)
    if len(a) < 2:
        return "NOT_RUN", "condition_n<2"
    if len(b) < 2:
        return "NOT_RUN", "control_n<2"
    sa = stdev(a)
    sb = stdev(b)
    if sa in (None, 0) and sb in (None, 0):
        return "NOT_RUN", "zero_variance_in_both_groups"
    return "OK", "Welch_t_test_on_log2_relative_expression"


def stats_condition_name(condition: Any) -> str:
    """Default biological grouping for mapped runs: ptfeb2/ptfeb3 -> ptfeb; k-3/k-4 -> k-."""
    s = str(condition).strip()
    m = re.match(r"^(.*?-)[0-9]+$", s)
    if m and m.group(1):
        return m.group(1)
    m = re.match(r"^(.+?)[0-9]+$", s)
    if m and m.group(1):
        return m.group(1)
    return s


def build_combined_summary_rows(rel_rows: List[Dict[str, Any]], target_mean_map: Optional[Dict[str, Optional[float]]] = None) -> List[Dict[str, Any]]:
    """Summarize already-normalized rows across analysis groups using stats_condition_name()."""
    target_mean_map = target_mean_map or {}
    out: List[Dict[str, Any]] = []
    reference = str(rel_rows[0].get("Reference", "")).strip() if rel_rows else ""
    targets = sorted({str(r.get("Target", "")).strip() for r in rel_rows if str(r.get("Target", "")).strip()})
    conditions = sorted({stats_condition_name(r.get("Condition")) for r in rel_rows if str(r.get("Condition", "")).strip()})
    controls = [stats_condition_name(r.get("control_condition")) for r in rel_rows if str(r.get("control_condition", "")).strip()]
    control = controls[0] if controls else ""
    if control in conditions:
        conditions = [control] + [c for c in conditions if c != control]
    for target in targets:
        for cond in conditions:
            sub = [r for r in rel_rows if normalize_key(r.get("Target")) == normalize_key(target) and normalize_key(stats_condition_name(r.get("Condition"))) == normalize_key(cond)]
            vals = [r.get("relative_expression") for r in sub]
            vals = [float(v) for v in vals if v is not None and math.isfinite(float(v))]
            logs = [r.get("log2_relative_expression") for r in sub]
            logs = [float(v) for v in logs if v is not None and math.isfinite(float(v))]
            if not vals:
                continue
            sd = stdev(vals)
            sd_log = stdev(logs)
            warnings = sorted({w for r in sub for w in str(r.get("QC_warning", "")).split(";") if w})
            out.append({
                "Analysis": "combined",
                "Target": target,
                "Condition": cond,
                "Reference": reference,
                "control_condition": control,
                "n_units": len(vals),
                "mean_relative_expression": mean(vals),
                "sd_relative_expression": sd,
                "sem_relative_expression": (sd / math.sqrt(len(vals))) if sd is not None and len(vals) > 0 else None,
                "mean_log2_relative_expression": mean(logs),
                "sd_log2_relative_expression": sd_log,
                "sem_log2_relative_expression": (sd_log / math.sqrt(len(logs))) if sd_log is not None and len(logs) > 0 else None,
                "E_target_used": target_mean_map.get(target),
                "Efficiency_percent_target_used": ((target_mean_map[target] - 1) * 100.0) if target in target_mean_map and target_mean_map.get(target) is not None else None,
                "E_reference_used": target_mean_map.get(reference),
                "Efficiency_percent_reference_used": ((target_mean_map[reference] - 1) * 100.0) if reference in target_mean_map and target_mean_map.get(reference) is not None else None,
                "QC_warning": ";".join(warnings),
            })
    return out


def compute_combined_relative_expression_stats(rel_rows: List[Dict[str, Any]], p_adjust: str = "bh") -> List[Dict[str, Any]]:
    """ANOVA and pairwise Welch tests after combining mapped biological runs by stats_condition_name()."""
    stats_rows: List[Dict[str, Any]] = []
    targets = sorted({str(r.get("Target", "")).strip() for r in rel_rows if str(r.get("Target", "")).strip()})
    controls = [stats_condition_name(r.get("control_condition")) for r in rel_rows if str(r.get("control_condition", "")).strip()]
    control = controls[0] if controls else ""
    for target in targets:
        trows = [r for r in rel_rows if normalize_key(r.get("Target")) == normalize_key(target)]
        conditions = sorted({stats_condition_name(r.get("Condition")) for r in trows if str(r.get("Condition", "")).strip()})
        if control in conditions:
            conditions = [control] + [c for c in conditions if c != control]
        values_by_cond = {}
        for cond in conditions:
            vals = [r.get("log2_relative_expression") for r in trows if normalize_key(stats_condition_name(r.get("Condition"))) == normalize_key(cond)]
            vals = [float(v) for v in vals if v is not None and math.isfinite(float(v))]
            values_by_cond[cond] = vals
        anova_status, anova_note = anova_status_from_groups(values_by_cond, control)
        anova_p = one_way_anova_p(list(values_by_cond.values())) if anova_status == "OK" else None
        if anova_status == "OK" and anova_p is None:
            anova_status, anova_note = "FAILED", "ANOVA_preconditions_ok_but_p_value_not_computed"
        raw = []
        conds = []
        pair_statuses: Dict[str, str] = {}
        pair_notes: Dict[str, str] = {}
        control_vals = values_by_cond.get(control, [])
        for cond in conditions:
            is_control = normalize_key(cond) == normalize_key(control)
            pair_status, pair_note = pairwise_status_for_groups(values_by_cond.get(cond, []), control_vals, is_control=is_control)
            p_raw = None if pair_status != "OK" else welch_ttest_p(values_by_cond.get(cond, []), control_vals)
            if pair_status == "OK" and p_raw is None:
                pair_status, pair_note = "FAILED", "Welch_preconditions_ok_but_p_value_not_computed"
            pair_statuses[cond] = pair_status
            pair_notes[cond] = pair_note
            raw.append(p_raw)
            conds.append(cond)
        adj = adjust_pvalues(raw, p_adjust)
        for cond, p_raw, p_adj in zip(conds, raw, adj):
            stats_rows.append({
                "Analysis": "combined",
                "Target": target,
                "Condition": cond,
                "control_condition": control,
                "n_units_for_stats": len(values_by_cond.get(cond, [])),
                "n_control_units_for_stats": len(control_vals),
                "anova_status": anova_status,
                "anova_note": anova_note,
                "anova_p": anova_p,
                "pairwise_test_status": pair_statuses.get(cond),
                "pairwise_test_note": pair_notes.get(cond),
                "p_vs_control_raw": p_raw,
                "p_vs_control_adjusted": p_adj,
                "p_adjust_method": p_adjust,
                "significance_vs_control": "control" if normalize_key(cond) == normalize_key(control) else p_to_stars(p_adj),
            })
    return stats_rows


def build_compact_output_combined(rel_rows: List[Dict[str, Any]], combined_summary_rows: List[Dict[str, Any]], stats_rows: List[Dict[str, Any]], target_eff_qc_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    summary_key = {(normalize_key(r.get("Target")), normalize_key(r.get("Condition"))): r for r in combined_summary_rows}
    stats_key = {(normalize_key(r.get("Target")), normalize_key(r.get("Condition"))): r for r in stats_rows}
    eff_key = {normalize_key(r.get("Target")): r for r in target_eff_qc_rows}
    out: List[Dict[str, Any]] = []
    for r in rel_rows:
        stats_cond = stats_condition_name(r.get("Condition"))
        key = (normalize_key(r.get("Target")), normalize_key(stats_cond))
        s = summary_key.get(key, {})
        st = stats_key.get(key, {})
        eff = eff_key.get(normalize_key(r.get("Target")), {})
        out.append({
            "Analysis": r.get("Analysis"),
            "Stats_condition": stats_cond,
            "Target": r.get("Target"),
            "Original_condition": r.get("Condition"),
            "BiologicalRep": r.get("BiologicalRep") or r.get("Analysis"),
            "Reference": r.get("Reference"),
            "Control_group": st.get("control_condition") or stats_condition_name(r.get("control_condition")),
            "Efficiency_mode": r.get("Efficiency_mode"),
            "Relative_expression": r.get("relative_expression"),
            "Log2_relative_expression": r.get("log2_relative_expression"),
            "Mean_relative_expression": s.get("mean_relative_expression"),
            "SD_relative_expression": s.get("sd_relative_expression"),
            "SEM_relative_expression": s.get("sem_relative_expression"),
            "n_biological_units": s.get("n_units"),
            "ANOVA_status": st.get("anova_status"),
            "ANOVA_note": st.get("anova_note"),
            "ANOVA_p": st.get("anova_p"),
            "Pairwise_test_status": st.get("pairwise_test_status"),
            "Pairwise_test_note": st.get("pairwise_test_note"),
            "p_vs_control_raw": st.get("p_vs_control_raw"),
            "p_vs_control_adjusted": st.get("p_vs_control_adjusted"),
            "Significance_vs_control": st.get("significance_vs_control"),
            "E_target_used": r.get("E_target_used"),
            "E_reference_used": r.get("E_reference_used"),
            "E_target_global_mean": r.get("E_target_global_mean"),
            "E_reference_global_mean": r.get("E_reference_global_mean"),
            "E_sample_target_mean": r.get("E_sample_target_mean"),
            "E_sample_target_sd": r.get("E_sample_target_sd"),
            "E_sample_target_cv_percent": r.get("E_sample_target_cv_percent"),
            "E_sample_reference_mean": r.get("E_sample_reference_mean"),
            "E_sample_reference_sd": r.get("E_sample_reference_sd"),
            "E_sample_reference_cv_percent": r.get("E_sample_reference_cv_percent"),
            "E_source_target": r.get("E_source_target"),
            "E_source_reference": r.get("E_source_reference"),
            "E_target_global_mean": r.get("E_target_global_mean"),
            "E_reference_global_mean": r.get("E_reference_global_mean"),
            "E_sample_target_mean": r.get("E_sample_target_mean"),
            "E_sample_target_sd": r.get("E_sample_target_sd"),
            "E_sample_target_cv_percent": r.get("E_sample_target_cv_percent"),
            "E_sample_reference_mean": r.get("E_sample_reference_mean"),
            "E_sample_reference_sd": r.get("E_sample_reference_sd"),
            "E_sample_reference_cv_percent": r.get("E_sample_reference_cv_percent"),
            "E_source_target": r.get("E_source_target"),
            "E_source_reference": r.get("E_source_reference"),
            "E_target_mean": s.get("E_target_used"),
            "Efficiency_percent_target": s.get("Efficiency_percent_target_used"),
            "E_reference_mean": s.get("E_reference_used"),
            "Efficiency_percent_reference": s.get("Efficiency_percent_reference_used"),
            "n_tech_target": r.get("n_tech_target"),
            "n_tech_reference": r.get("n_tech_reference"),
            "mean_Cq_target": r.get("mean_Cq_target"),
            "mean_Cq_reference": r.get("mean_Cq_reference"),
            "target_wells": r.get("target_wells"),
            "reference_wells": r.get("reference_wells"),
            "Efficiency_QC": eff.get("outlier_method"),
            "QC_warning": r.get("QC_warning") or s.get("QC_warning") or "",
        })
    return out




def build_rstudio_input_rows(rel_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Minimal long-format table for external statistics in R/RStudio/GraphPad.

    One row = one biological unit/sample x one target gene.

    Columns are intentionally limited to what is needed for ANOVA:
      sample                unique biological unit name, e.g. empty1, k-2
      condition             collapsed group used for statistics, e.g. empty, k-, ptfeb
      target                target gene
      relative_expression   final efficiency-corrected fold expression
      log2_relative_expression  log2 of final expression; usually preferable for ANOVA

    QC, Cq, efficiency and well-level information remain in the compact/QC files,
    not in this statistics input table.
    """
    out: List[Dict[str, Any]] = []
    for r in rel_rows:
        condition_original = str(r.get("Condition", "")).strip()
        stats_condition = stats_condition_name(condition_original)
        bio_rep = str(r.get("BiologicalRep", "")).strip()
        if not bio_rep:
            bio_rep = str(r.get("Analysis", "")).strip()
        sample = str(r.get("Unit") or "").strip()
        if not sample:
            sample = f"{condition_original}|{bio_rep}" if bio_rep else condition_original

        out.append({
            "sample": sample,
            "condition": stats_condition,
            "target": r.get("Target"),
            "relative_expression": r.get("relative_expression"),
            "log2_relative_expression": r.get("log2_relative_expression"),
        })

    out.sort(key=lambda x: (str(x.get("target")), str(x.get("condition")), str(x.get("sample"))))
    return out


def build_sample_target_efficiency_qc_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Compact QC table: one row = one sample x target efficiency estimate."""
    grouped: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = {}
    for r in rows:
        target = str(r.get("Target", "")).strip()
        condition = str(r.get("Condition", "")).strip()
        biorep = str(r.get("BiologicalRep", "")).strip()
        if target and condition:
            grouped.setdefault((condition, biorep, target), []).append(r)
    out: List[Dict[str, Any]] = []
    for (condition, biorep, target), sub in sorted(grouped.items(), key=lambda x: (x[0][2].lower(), x[0][0].lower(), x[0][1].lower())):
        e_vals = []
        wells = []
        statuses = []
        for r in sub:
            e = r.get("E_amp_factor_per_well")
            if e is not None and math.isfinite(float(e)):
                e_vals.append(float(e))
            wells.append(str(r.get("Well", "")))
            statuses.append(str(r.get("Efficiency_status", "")))
        em = mean(e_vals)
        esd = stdev(e_vals)
        cv = (100.0 * esd / em) if em not in (None, 0) and esd is not None else None
        global_e = sub[0].get("E_target_global_mean")
        ratio = (em / global_e) if em is not None and global_e not in (None, 0) else None
        delta = (em - global_e) if em is not None and global_e is not None else None
        qc = []
        if cv is not None and cv > 10.0:
            qc.append("sample_target_E_CV_gt_10pct")
        if ratio is not None and abs(ratio - 1.0) > 0.10:
            qc.append("sample_target_E_differs_from_global_gt_10pct")
        if any(s and s != "OK" for s in statuses):
            qc.append("one_or_more_curve_E_warnings")
        out.append({
            "condition": condition,
            "stats_condition": stats_condition_name(condition),
            "biological_rep": biorep,
            "target": target,
            "E_sample_target_mean": em,
            "Efficiency_percent_sample_target_mean": ((em - 1.0) * 100.0) if em is not None else None,
            "E_sample_target_sd": esd,
            "E_sample_target_cv_percent": cv,
            "n_E_technical": len(e_vals),
            "E_target_global_mean": global_e,
            "Efficiency_percent_target_global_mean": ((global_e - 1.0) * 100.0) if global_e is not None else None,
            "E_sample_vs_target_mean_ratio": ratio,
            "E_sample_vs_target_mean_delta": delta,
            "efficiency_statuses": ";".join(sorted(set(statuses))),
            "wells": ";".join(wells),
            "QC_warning": ";".join(qc),
        })
    return out


def compute_relative_expression_stats(rel_rows: List[Dict[str, Any]], p_adjust: str = "bh") -> List[Dict[str, Any]]:
    """ANOVA by target/analysis and pairwise Welch tests against control on log2(relative_expression)."""
    stats_rows: List[Dict[str, Any]] = []
    analyses = sorted({str(r.get("Analysis", "analysis_1")).strip() or "analysis_1" for r in rel_rows})
    for analysis in analyses:
        arows = [r for r in rel_rows if (str(r.get("Analysis", "analysis_1")).strip() or "analysis_1") == analysis]
        targets = sorted({str(r.get("Target", "")).strip() for r in arows if str(r.get("Target", "")).strip()})
        for target in targets:
            trows = [r for r in arows if normalize_key(r.get("Target")) == normalize_key(target)]
            if not trows:
                continue
            control = str(trows[0].get("control_condition", "")).strip()
            conditions = sorted({str(r.get("Condition", "")).strip() for r in trows if str(r.get("Condition", "")).strip()})
            if control in conditions:
                conditions = [control] + [c for c in conditions if c != control]
            values_by_cond = {}
            for cond in conditions:
                vals = [r.get("log2_relative_expression") for r in trows if normalize_key(r.get("Condition")) == normalize_key(cond)]
                vals = [float(v) for v in vals if v is not None and math.isfinite(float(v))]
                values_by_cond[cond] = vals
            anova_status, anova_note = anova_status_from_groups(values_by_cond, control)
            anova_p = one_way_anova_p(list(values_by_cond.values())) if anova_status == "OK" else None
            if anova_status == "OK" and anova_p is None:
                anova_status, anova_note = "FAILED", "ANOVA_preconditions_ok_but_p_value_not_computed"
            raw_p = []
            row_indices = []
            control_vals = values_by_cond.get(control, [])
            pair_statuses: Dict[str, str] = {}
            pair_notes: Dict[str, str] = {}
            for cond in conditions:
                vals = values_by_cond.get(cond, [])
                is_control = normalize_key(cond) == normalize_key(control)
                pair_status, pair_note = pairwise_status_for_groups(vals, control_vals, is_control=is_control)
                p_raw = None if pair_status != "OK" else welch_ttest_p(vals, control_vals)
                if pair_status == "OK" and p_raw is None:
                    pair_status, pair_note = "FAILED", "Welch_preconditions_ok_but_p_value_not_computed"
                pair_statuses[cond] = pair_status
                pair_notes[cond] = pair_note
                row_indices.append((cond, p_raw))
                raw_p.append(p_raw)
            adj_p = adjust_pvalues(raw_p, p_adjust)
            for (cond, p_raw), p_adj in zip(row_indices, adj_p):
                stats_rows.append({
                    "Analysis": analysis,
                    "Target": target,
                    "Condition": cond,
                    "control_condition": control,
                    "n_units_for_stats": len(values_by_cond.get(cond, [])),
                    "n_control_units_for_stats": len(control_vals),
                    "anova_status": anova_status,
                    "anova_note": anova_note,
                    "anova_p": anova_p,
                    "pairwise_test_status": pair_statuses.get(cond),
                    "pairwise_test_note": pair_notes.get(cond),
                    "p_vs_control_raw": p_raw,
                    "p_vs_control_adjusted": p_adj,
                    "p_adjust_method": p_adjust,
                    "significance_vs_control": "control" if normalize_key(cond) == normalize_key(control) else p_to_stars(p_adj),
                })
    return stats_rows


def build_compact_output(rel_rows: List[Dict[str, Any]], summary_rows: List[Dict[str, Any]], stats_rows: List[Dict[str, Any]], target_eff_qc_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    summary_key = {(r.get("Analysis"), normalize_key(r.get("Target")), normalize_key(r.get("Condition"))): r for r in summary_rows}
    stats_key = {(r.get("Analysis"), normalize_key(r.get("Target")), normalize_key(r.get("Condition"))): r for r in stats_rows}
    eff_key = {normalize_key(r.get("Target")): r for r in target_eff_qc_rows}
    out: List[Dict[str, Any]] = []
    for r in rel_rows:
        key = (r.get("Analysis"), normalize_key(r.get("Target")), normalize_key(r.get("Condition")))
        s = summary_key.get(key, {})
        st = stats_key.get(key, {})
        eff = eff_key.get(normalize_key(r.get("Target")), {})
        out.append({
            "Analysis": r.get("Analysis"),
            "Target": r.get("Target"),
            "Condition": r.get("Condition"),
            "BiologicalRep": r.get("BiologicalRep"),
            "Reference": r.get("Reference"),
            "Control": r.get("control_condition"),
            "Relative_expression": r.get("relative_expression"),
            "Log2_relative_expression": r.get("log2_relative_expression"),
            "Mean_relative_expression": s.get("mean_relative_expression"),
            "SD_relative_expression": s.get("sd_relative_expression"),
            "SEM_relative_expression": s.get("sem_relative_expression"),
            "n_biological_units": s.get("n_units"),
            "ANOVA_status": st.get("anova_status"),
            "ANOVA_note": st.get("anova_note"),
            "ANOVA_p": st.get("anova_p"),
            "Pairwise_test_status": st.get("pairwise_test_status"),
            "Pairwise_test_note": st.get("pairwise_test_note"),
            "p_vs_control_raw": st.get("p_vs_control_raw"),
            "p_vs_control_adjusted": st.get("p_vs_control_adjusted"),
            "Significance_vs_control": st.get("significance_vs_control"),
            "E_target_mean": s.get("E_target_used"),
            "Efficiency_percent_target": s.get("Efficiency_percent_target_used"),
            "E_reference_mean": s.get("E_reference_used"),
            "Efficiency_percent_reference": s.get("Efficiency_percent_reference_used"),
            "n_tech_target": r.get("n_tech_target"),
            "n_tech_reference": r.get("n_tech_reference"),
            "mean_Cq_target": r.get("mean_Cq_target"),
            "mean_Cq_reference": r.get("mean_Cq_reference"),
            "target_wells": r.get("target_wells"),
            "reference_wells": r.get("reference_wells"),
            "Efficiency_QC": eff.get("outlier_method"),
            "QC_warning": r.get("QC_warning") or s.get("QC_warning") or "",
        })
    return out


def plot_publication_with_stats(summary_rows: List[Dict[str, Any]], stats_rows: List[Dict[str, Any]], out_path: str, metric: str = "relative", error: str = "sd", title: Optional[str] = None) -> bool:
    """Publication-style grouped bar plot WITHOUT statistical labels."""
    rows = [r for r in summary_rows if str(r.get("Target", "")).strip() and str(r.get("Condition", "")).strip()]
    if not rows:
        return False
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as e:
        print(f"WARNING: plot skipped: matplotlib/numpy is not available: {e}", file=sys.stderr)
        return False

    mean_col = "mean_log2_relative_expression" if metric == "log2" else "mean_relative_expression"
    sd_col = "sd_log2_relative_expression" if metric == "log2" else "sd_relative_expression"
    sem_col = "sem_log2_relative_expression" if metric == "log2" else "sem_relative_expression"
    err_col = sem_col if error == "sem" else (sd_col if error == "sd" else None)
    y_label = "log2 relative expression" if metric == "log2" else "Relative expression"
    baseline = 0 if metric == "log2" else 1

    targets = sorted({str(r.get("Target", "")).strip() for r in rows})
    conditions = sorted({str(r.get("Condition", "")).strip() for r in rows})
    controls = [str(r.get("control_condition", "")).strip() for r in rows if str(r.get("control_condition", "")).strip()]
    control = controls[0] if controls else None
    if control in conditions:
        conditions = [control] + [c for c in conditions if c != control]

    stat_key: Dict[Tuple[str, str], Dict[str, Any]] = {}  # no statistical labels on the plot
    data = {(str(r.get("Target", "")).strip(), str(r.get("Condition", "")).strip()): r for r in rows}

    n_targets = len(targets)
    n_conditions = len(conditions)
    x = np.arange(n_targets, dtype=float)
    total_width = 0.78
    bar_width = total_width / max(n_conditions, 1)
    fig_width = max(7.0, 1.1 * n_targets * max(n_conditions, 2))
    fig, ax = plt.subplots(figsize=(fig_width, 5.0), dpi=600)

    max_y = baseline
    for j, cond in enumerate(conditions):
        offsets = x - total_width / 2 + bar_width / 2 + j * bar_width
        means, errors = [], []
        for target in targets:
            row = data.get((target, cond), {})
            val = to_float(row.get(mean_col))
            err = to_float(row.get(err_col)) if err_col else None
            means.append(val if val is not None else float("nan"))
            errors.append(err if err is not None else 0.0)
            if val is not None:
                max_y = max(max_y, val + (err or 0.0))
        ax.bar(offsets, means, width=bar_width, label=cond, yerr=errors if err_col else None, capsize=3 if err_col else 0, linewidth=0.8, edgecolor="black")
        for xpos, target, val, err in zip(offsets, targets, means, errors):
            if not math.isfinite(val):
                continue
            st = stat_key.get((normalize_key(target), normalize_key(cond)), {})
            stars = st.get("significance_vs_control", "")
            if stars and stars not in {"control", "NA", "ns"}:
                y = val + err + (0.04 * max(1.0, max_y))
                ax.text(xpos, y, stars, ha="center", va="bottom", fontsize=10)
                max_y = max(max_y, y)

    ax.axhline(baseline, linewidth=0.9, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(targets, fontsize=10)
    ax.set_ylabel(y_label, fontsize=11)
    ax.set_title(title or "Efficiency-corrected qPCR expression", fontsize=12)
    ax.legend(title="Condition", frameon=False, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=9)
    if metric == "relative":
        ax.set_ylim(bottom=0, top=max_y * 1.20 if max_y > 0 else 1.5)
    else:
        ax.set_ylim(top=max_y * 1.20 if max_y > 0 else max_y + 1)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return True

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Bio-Rad qPCR pipeline: robust LinReg-style efficiency, analysis-map normalization, compact publication output."
    )
    ap.add_argument("--amp", required=True, nargs="+", help="One or more Bio-Rad Quantification Amplification Results.xlsx files. If several are given, --cq must contain the same number of files in the same order.")
    ap.add_argument("--cq", required=True, nargs="+", help="One or more Bio-Rad Quantification Cq Results.xlsx files, in the same order as --amp.")
    ap.add_argument("--sheet", default="auto", help="Amplification sheet name. Use auto/all to read all amplification sheets except Run Information. Default: auto")
    ap.add_argument("--cq-sheet", default=None, help="Cq sheet name. Default: first sheet, often '0'")
    ap.add_argument("--out", default=None, help="Output prefix. Default: based on amplification filename")

    ap.add_argument("--min-window", type=int, default=5, help="Minimum log-linear window size. Default: 5")
    ap.add_argument("--max-window", type=int, default=7, help="Maximum log-linear window size. Default: 7")
    ap.add_argument("--min-e", type=float, default=1.20, help="Low-E warning threshold. Low E is annotated, not excluded. Default: 1.20")
    ap.add_argument("--max-e", type=float, default=2.40, help="High-E warning threshold. High E is annotated, not excluded. Default: 2.40")
    ap.add_argument("--good-r2", type=float, default=0.990, help="R2 warning threshold. Low R2 is annotated, not excluded. Default: 0.990")

    ap.add_argument("--well-col", default="Well", help="Well column in Cq file. Default: Well")
    ap.add_argument("--cq-col", default="Cq", help="Cq column in Cq file. Default: Cq")
    ap.add_argument("--target-col", default="Target", help="Target/gene column in Cq file. Default: Target")
    ap.add_argument("--target-fallback-col", default=None, help="Use this column if target-col is blank, e.g. Sample")
    ap.add_argument("--sample-col", default="Sample", help="Sample column in Cq file. Default: Sample")
    ap.add_argument("--condition-col", default="Sample", help="Experimental condition column. Default: Sample")
    ap.add_argument("--condition-regex", default=None, help="Optional regex to extract condition from condition-col.")
    ap.add_argument("--condition-regex-group", type=int, default=1, help="Regex capture group to use for --condition-regex. Default: 1")
    ap.add_argument("--biological-rep-col", default=None, help="Biological replicate column, e.g. Biological Set Name or SampleID")
    ap.add_argument("--analysis-col", default=None, help="Optional column in the Cq/layout file that marks independent experiments/analysis groups.")
    ap.add_argument("--content-col", default="Content", help="Content column. Default: Content")
    ap.add_argument("--fluor-col", default="Fluor", help="Fluor column. Default: Fluor")

    ap.add_argument("--reference", default=None, help="Reference gene/target name, e.g. PPIA. Required for relative expression.")
    ap.add_argument("--targets", nargs="*", default=None, help="Target genes. If omitted, all non-reference targets are used.")
    ap.add_argument("--control-condition", default=None, help="Control condition name for fold-change normalization. Used when --analysis-map is not provided.")
    ap.add_argument("--analysis-map", default=None, help="CSV/TSV file defining independent analyses. Required columns: Analysis, Condition, ControlCondition.")
    ap.add_argument("--control-map", default=None, help="Controls for --analysis-col mode, format: analysis1=control1;analysis2=control2")
    ap.add_argument("--efficiency-mode", choices=["rdml_like", "target_mean", "per_well", "sample_target_mean", "fixed", "fixed_target"], default="rdml_like", help="Default: rdml_like. Use fixed with --fixed-efficiency 2.0, or fixed_target with --fixed-target-efficiency ppia=1.542,atg9b=1.739")
    ap.add_argument("--fixed-efficiency", type=float, default=2.0, help="Amplification factor used only with --efficiency-mode fixed. Example: 2.0 means 100%% efficiency. Default: 2.0")
    ap.add_argument("--fixed-target-efficiency", default=None, help="Target-specific amplification factors used only with --efficiency-mode fixed_target. Example: ppia=1.542,atg9b=1.739,lamp1=1.735,p62=1.698. Percent values like ppia=54.2%% are also accepted.")
    ap.add_argument("--no-calc-cq-from-curves", action="store_true", help="Do not calculate missing Cq from amplification curves. Recommended when Bio-Rad Cq is available.")

    ap.add_argument("--efficiency-outlier-mode", choices=["iqr", "mad", "none"], default="none", help="Optional outlier filter for individual efficiencies. Default: none.")
    ap.add_argument("--efficiency-iqr-k", type=float, default=1.5, help="IQR multiplier for optional efficiency outlier filtering. Default: 1.5")
    ap.add_argument("--efficiency-mad-k", type=float, default=3.5, help="MAD multiplier for optional efficiency outlier filtering. Default: 3.5")
    ap.add_argument("--p-adjust", choices=["bh", "holm", "bonferroni", "none"], default="bh", help="Multiple-testing correction for pairwise tests vs control. Default: bh")
    ap.add_argument("--write-intermediate", action="store_true", help="Write old intermediate CSV files. Default: off; main output is one compact CSV.")
    ap.add_argument("--debug", action="store_true", help="Only affects intermediate diagnostic output when --write-intermediate is used.")

    ap.add_argument("--no-plot", action="store_true", help="Do not save a relative-expression plot.")
    ap.add_argument("--plot-metric", choices=["relative", "log2"], default="relative", help="Plot relative or log2 expression. Default: relative")
    ap.add_argument("--plot-error", choices=["sem", "sd", "none"], default="sd", help="Error bars for the plot. Default: sd")
    ap.add_argument("--plot-title", default=None, help="Optional custom plot title.")
    ap.add_argument("--plot-file", default=None, help="Optional output image path. Default: <out>_publication_plot_<analysis>.png")

    args = ap.parse_args()
    print("Script version: linregpcr_strict_window_multiinput_target_casefix_manual_efficiency_annotated_per_well_2026-06-14")

    if len(args.amp) != len(args.cq):
        print(f"ERROR: --amp and --cq must contain the same number of files. Got {len(args.amp)} amp and {len(args.cq)} cq.", file=sys.stderr)
        return 2
    for amp_path in args.amp:
        if not os.path.exists(amp_path):
            print(f"ERROR: amplification file not found: {amp_path}", file=sys.stderr)
            return 2
    for cq_path in args.cq:
        if not os.path.exists(cq_path):
            print(f"ERROR: Cq file not found: {cq_path}", file=sys.stderr)
            return 2
    if args.reference and not args.control_condition and not args.analysis_map:
        print("ERROR: when --reference is provided, use either --control-condition or --analysis-map", file=sys.stderr)
        return 2
    if args.analysis_map and args.control_map:
        print("WARNING: --control-map is ignored when --analysis-map is used", file=sys.stderr)

    prefix = safe_out_prefix(args.amp[0], args.out)

    # 1-2. Efficiency per well and Cq merge. Multiple qPCR files are processed
    # as independent Bio-Rad exports, then combined before relative expression/statistics.
    all_eff_rows: List[Dict[str, Any]] = []
    all_merged_rows: List[Dict[str, Any]] = []
    total_file_cq_count = 0
    total_calc_cq_count = 0

    for dataset_index, (amp_path, cq_path) in enumerate(zip(args.amp, args.cq), start=1):
        ds_args = argparse.Namespace(**vars(args))
        ds_args.amp = amp_path
        ds_args.cq = cq_path
        dataset_name = os.path.splitext(os.path.basename(amp_path))[0]
        print(f"\nDataset {dataset_index}: {dataset_name}")

        eff_rows = compute_per_well_efficiency(ds_args)
        for r in eff_rows:
            r["Dataset"] = dataset_index
            r["Amp_File"] = os.path.basename(amp_path)
            r["Cq_File"] = os.path.basename(cq_path)
        n = len(eff_rows)
        n_ok = sum(1 for r in eff_rows if r.get("status") == "OK")
        n_warn = sum(1 for r in eff_rows if r.get("status") == "QC_WARNING")
        n_no = sum(1 for r in eff_rows if r.get("E_amp_factor") is None)
        print(f"Wells in amplification file: {n}; OK: {n_ok}; QC warnings kept: {n_warn}; no E: {n_no}")

        merged_rows, _header = build_merged_cq_efficiency(ds_args, eff_rows)
        for r in merged_rows:
            r["Dataset"] = dataset_index
            r["Amp_File"] = os.path.basename(amp_path)
            r["Cq_File"] = os.path.basename(cq_path)
        file_cq_count = sum(1 for r in merged_rows if r.get("Cq") is not None)
        calculate_missing_cq_from_amp_fit(merged_rows, enabled=(not args.no_calc_cq_from_curves))
        calc_cq_count = sum(1 for r in merged_rows if r.get("Cq_source") == "calculated_from_amp_curve")
        print(f"Cq rows: {len(merged_rows)}; Cq from file: {file_cq_count}; Cq calculated: {calc_cq_count}; valid Cq total: {sum(1 for r in merged_rows if r.get('Cq') is not None)}")

        all_eff_rows.extend(eff_rows)
        all_merged_rows.extend(merged_rows)
        total_file_cq_count += file_cq_count
        total_calc_cq_count += calc_cq_count

    eff_rows = all_eff_rows
    merged_rows = all_merged_rows

    # IMPORTANT: Bio-Rad exports from different runs may spell the same target
    # differently (for example PPIA vs ppia, Atg9b vs atg9b, LAMP1 vs lamp1).
    # Earlier versions resolved the reference name case-insensitively once, but then
    # aggregated technical quantities with case-sensitive Target keys. As a result,
    # a whole qPCR run could be read and included in efficiency QC, but silently fail
    # to enter relative expression if its reference target spelling differed from
    # the first run. Normalize all Target names before efficiency assignment and
    # relative-expression aggregation. Keep the original spelling for diagnostics.
    target_spellings: Dict[str, set] = {}
    for r in merged_rows:
        original_target = str(r.get("Target", "")).strip()
        if original_target:
            canonical_target = normalize_key(original_target)
            r["Target_original"] = original_target
            r["Target"] = canonical_target
            target_spellings.setdefault(canonical_target, set()).add(original_target)
    if args.reference:
        args.reference = normalize_key(args.reference)
    if args.targets:
        args.targets = [normalize_key(t) for t in args.targets]
    changed_spellings = {k: sorted(v) for k, v in target_spellings.items() if len(v) > 1}
    if changed_spellings:
        print("Target names normalized across qPCR files:")
        for canonical, variants in sorted(changed_spellings.items()):
            print(f"  {canonical}: " + ", ".join(variants))

    # Add gene/sample annotations from the Cq/layout table back to the per-well
    # efficiency rows. This does not change calculations; it only makes the
    # *_per_well_efficiency.csv QC table interpretable by target and sample.
    annotate_efficiency_rows_from_cq(eff_rows, merged_rows)

    target_mean_map, target_eff_qc_rows = assign_efficiency_by_mode(merged_rows, args.efficiency_mode, args)
    for r in merged_rows:
        target = str(r.get("Target", "")).strip()
        r["E_target_mean"] = target_mean_map.get(target)
        r["Efficiency_percent_target_mean"] = ((target_mean_map[target] - 1) * 100) if target in target_mean_map and target_mean_map[target] is not None else None

    print(f"\nCombined datasets: {len(args.amp)}")
    print(f"Combined Cq rows: {len(merged_rows)}; Cq from file: {total_file_cq_count}; Cq calculated: {total_calc_cq_count}; valid Cq total: {sum(1 for r in merged_rows if r.get('Cq') is not None)}")

    if args.write_intermediate:
        eff_path = prefix + "_per_well_efficiency.csv"
        eff_fields = [
            "Dataset", "Amp_File", "Cq_File", "Well", "Amp_Sheet",
            "Target", "Target_original", "Sample", "Condition", "BiologicalRep", "Analysis", "Content", "Fluor",
            "Cq", "Cq_source", "Cq_missing", "Cq_annotation_found",
            "status", "E_amp_factor", "Efficiency_percent", "r2", "slope", "intercept",
            "cycle_start", "cycle_end", "window_cycles", "window_n", "residual_sd", "baseline", "baseline_cycles",
            "noise", "max_delta", "threshold_low", "threshold_high", "derivative_peak_cycle", "derivative_peak_slope",
            "slope_quality_vs_peak", "self_check_status", "self_check_details", "baseline_slope_balance",
        ]
        write_csv(eff_path, eff_rows, eff_fields)
        print(f"Saved: {eff_path}")

        target_eff_path = prefix + "_qc_target_efficiency.csv"
        target_eff_fields = [
            "Target", "E_target_mean", "Efficiency_percent_target_mean",
            "n_efficiency_OK_raw", "n_efficiency_kept_after_outlier_filter",
            "n_efficiency_rejected_bad_curve", "n_efficiency_rejected_outlier",
            "sd_E_kept", "cv_E_percent_kept", "outlier_method", "outlier_low", "outlier_high", "kept_wells", "outlier_wells",
        ]
        write_csv(target_eff_path, target_eff_qc_rows, target_eff_fields)
        print(f"Saved: {target_eff_path}")

        if args.debug:
            merged_path = prefix + "_debug_cq_merged_with_efficiency.csv"
            merged_fields = [
                "Dataset", "Amp_File", "Cq_File", "Well", "Amp_Sheet", "Fluor", "Target", "Sample", "Condition", "BiologicalRep", "Analysis", "Content",
                "Cq", "Cq_original", "Cq_source", "Cq_missing", "Cq_auto_log10_threshold", "Cq_auto_method",
                "E_amp_factor_per_well", "Efficiency_percent_per_well", "Efficiency_status", "Efficiency_self_check",
                "Eff_slope", "Eff_intercept", "Eff_cycle_start", "Eff_cycle_end", "Eff_window_cycles",
                "E_target_mean", "Efficiency_percent_target_mean", "E_used", "E_source", "Quantity_E_minus_Cq",
            ]
            write_csv(merged_path, merged_rows, merged_fields)
            print(f"Saved: {merged_path}")

    # 3. Relative expression + statistics + compact output.
    if args.reference:
        rel_rows, summary_rows_per_analysis = aggregate_relative_expression(args, merged_rows, target_mean_map)
        # For statistics, combine biological runs after each run has been normalized to its own control.
        # Example: ptfeb2 + ptfeb3 -> ptfeb; k-3 + k-4 -> k-.
        summary_rows = build_combined_summary_rows(rel_rows, target_mean_map)
        stats_rows: List[Dict[str, Any]] = []  # statistics are intentionally not calculated in this version; use R/RStudio
        compact_rows = build_compact_output_combined(rel_rows, summary_rows, stats_rows, target_eff_qc_rows)
        compact_path = prefix + "_compact_results.csv"
        compact_fields = [
            "Analysis", "Stats_condition", "Target", "Original_condition", "BiologicalRep", "Reference", "Control_group", "Efficiency_mode",
            "Relative_expression", "Log2_relative_expression", "Mean_relative_expression", "SD_relative_expression", "SEM_relative_expression", "n_biological_units",
            "E_target_used", "E_reference_used", "E_target_global_mean", "E_reference_global_mean",
            "E_sample_target_mean", "E_sample_target_sd", "E_sample_target_cv_percent",
            "E_sample_reference_mean", "E_sample_reference_sd", "E_sample_reference_cv_percent",
            "E_source_target", "E_source_reference",
            "E_target_mean", "Efficiency_percent_target", "E_reference_mean", "Efficiency_percent_reference",
            "n_tech_target", "n_tech_reference", "mean_Cq_target", "mean_Cq_reference", "target_wells", "reference_wells", "Efficiency_QC", "QC_warning",
        ]
        write_csv(compact_path, compact_rows, compact_fields)
        print(f"Relative-expression rows: {len(rel_rows)}")
        if len(rel_rows) == 0:
            print("WARNING: no relative expression rows were produced. Check Cq/reference/condition/analysis_map names.", file=sys.stderr)
        print(f"Saved compact result: {compact_path}")

        r_input_rows = build_rstudio_input_rows(rel_rows)
        r_input_path = prefix + "_rstats_input.csv"
        r_input_fields = [
            "sample", "condition", "target", "relative_expression", "log2_relative_expression",
        ]
        write_csv(r_input_path, r_input_rows, r_input_fields)
        print(f"Saved R/RStudio input table: {r_input_path}")

        sample_eff_path = prefix + "_sample_target_efficiency_qc.csv"
        sample_eff_rows = build_sample_target_efficiency_qc_rows(merged_rows)
        sample_eff_fields = [
            "condition", "stats_condition", "biological_rep", "target",
            "E_sample_target_mean", "Efficiency_percent_sample_target_mean",
            "E_sample_target_sd", "E_sample_target_cv_percent", "n_E_technical",
            "E_target_global_mean", "Efficiency_percent_target_global_mean",
            "E_sample_vs_target_mean_ratio", "E_sample_vs_target_mean_delta",
            "efficiency_statuses", "wells", "QC_warning",
        ]
        write_csv(sample_eff_path, sample_eff_rows, sample_eff_fields)
        print(f"Saved sample-target efficiency QC table: {sample_eff_path}")

        if args.write_intermediate:
            write_csv(prefix + "_final_relative_expression_by_unit.csv", rel_rows, [
                "Analysis", "Target", "Condition", "BiologicalRep", "Unit", "Reference", "control_condition",
                "relative_expression", "log2_relative_expression", "normalized_quantity_target_over_ref",
                "mean_Cq_target", "sd_Cq_target", "mean_Cq_reference", "sd_Cq_reference",
                "E_target_used", "E_reference_used", "n_tech_target", "n_tech_reference",
                "target_wells", "reference_wells", "QC_warning",
            ])
            write_csv(prefix + "_final_summary_by_condition.csv", summary_rows, [
                "Analysis", "Target", "Condition", "Reference", "control_condition", "n_units",
                "mean_relative_expression", "sd_relative_expression", "sem_relative_expression",
                "mean_log2_relative_expression", "sd_log2_relative_expression", "sem_log2_relative_expression",
                "E_target_used", "Efficiency_percent_target_used", "E_reference_used", "Efficiency_percent_reference_used", "QC_warning",
            ])
            write_csv(prefix + "_stats.csv", stats_rows, [
                "Analysis", "Target", "Condition", "control_condition", "n_units_for_stats", "n_control_units_for_stats", "anova_status", "anova_note", "anova_p", "pairwise_test_status", "pairwise_test_note", "p_vs_control_raw", "p_vs_control_adjusted", "p_adjust_method", "significance_vs_control",
            ])

        if not args.no_plot:
            plot_path = args.plot_file or (prefix + "_publication_plot.png")
            plotted = plot_publication_with_stats(summary_rows, stats_rows, plot_path, metric=args.plot_metric, error=args.plot_error, title=args.plot_title)
            if plotted:
                print(f"Saved plot: {plot_path}")
    else:
        compact_path = prefix + "_compact_cq_efficiency.csv"
        write_csv_dynamic(compact_path, merged_rows)
        print(f"No --reference provided. Saved compact Cq+efficiency table: {compact_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
