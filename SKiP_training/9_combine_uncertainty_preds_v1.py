#Code to combine aleatoric and epistemic variances into a single estimate of prediction uncertainty

import os
import sys
import csv
import pandas as pd
import numpy as np
from scipy.special import erfinv
from typing import Tuple, Union, Literal

infile = sys.argv[1]
outfile = sys.argv[2]

#XXX: python 9_combine_uncertainty_preds_v1.py "./Sample_output/kcat_X_test_uncertainty_sample.csv" "./Sample_output/kcat_X_test_uncertainty_combined_sample.csv"

#XXX: AI Usage Disclaimer: Code generated with help from Microsoft CoPilot
def combine_uncertainties_from_mapie(
    y_pred: np.ndarray,
    y_lower: np.ndarray,
    y_upper: np.ndarray,
    alpha: float,
    epistemic_var: Union[float, np.ndarray],
    *,
    asymmetry: Literal["auto", "symmetric", "average"] = "auto",
    clip_nonnegative: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

    """
    Convert a MAPIE prediction interval (for a single alpha) into an
    aleatoric variance under a Normal approximation, then combine with
    an epistemic variance via variance addition.

    Parameters
    ----------
    y_pred : (n,) array
        MAPIE point predictions for the samples.
    y_lower : (n,) array
        Lower bounds of the MAPIE prediction interval at the given alpha.
    y_upper : (n,) array
        Upper bounds of the MAPIE prediction interval at the given alpha.
    alpha : float
        Miscoverage level used to compute the MAPIE interval (e.g., 0.1 for 90% PI).
        The function converts the half-width to sigma using z_{1 - alpha/2}.
    epistemic_var : float or (n,) array
        Epistemic variance per sample (or a single scalar to broadcast). Must be >= 0.
        Typical sources: ensembles, MC-dropout, Bayesian model variance.
    asymmetry : {"auto", "symmetric", "average"}, default "auto"
        How to handle asymmetric intervals:
        - "symmetric": assume PI ≈ y_pred ± w, i.e., use half-width (U-L)/2.
        - "average": compute two one-sided widths and average their sigma estimates:
            sigma = 0.5 * ( (U - y_pred)/z + (y_pred - L)/z )
        - "auto": use "average" if the interval is visibly asymmetric for a sample,
                  else "symmetric". Asymmetry is detected by relative difference
                  between upper and lower half-widths (> 10% by default).
    clip_nonnegative : bool, default True
        If True, clips any tiny negative variances to zero after numerical noise.

    Returns
    -------
    sigma_aleatoric : (n,) array
        Estimated aleatoric standard deviation (Normal approximation).
    var_aleatoric : (n,) array
        Estimated aleatoric variance.
    var_total : (n,) array
        Total variance = epistemic_var + var_aleatoric.
    sigma_total : (n,) array
        Total standard deviation = sqrt(var_total).

    Notes
    -----
    - This function *assumes* the PI can be mapped to a Normal central interval:
         PI_{1-alpha} ≈ [ y_pred - z * sigma,  y_pred + z * sigma ],
      where z = Phi^{-1}(1 - alpha/2).
    - With alpha=0.1 (90% PI), z ≈ 1.645; with alpha≈0.3173 (~68% PI), z≈1.
    - If MAPIE intervals are strongly asymmetric or heavy-tailed, prefer
      reporting the intervals directly or compute PIs at multiple alpha values
      to check scaling before compressing into a single σ.
    """
    
    y_pred = np.asarray(y_pred).ravel()
    y_lower = np.asarray(y_lower).ravel()
    y_upper = np.asarray(y_upper).ravel()

    if y_pred.shape != y_lower.shape or y_pred.shape != y_upper.shape:
        raise ValueError("y_pred, y_lower, and y_upper must have the same shape (n,).")

    # Broadcast epistemic variance
    epistemic_var = np.asarray(epistemic_var)
    if epistemic_var.ndim == 0:
        epistemic_var = np.full_like(y_pred, float(epistemic_var))
    elif epistemic_var.shape != y_pred.shape:
        raise ValueError("epistemic_var must be a scalar or shape-matched (n,) array.")

    if np.any(epistemic_var < 0):
        raise ValueError("epistemic_var must be nonnegative.")

    # Guard rails on alpha
    if not (0 < alpha < 1):
        raise ValueError("alpha must be in (0, 1).")

    # Normal z for central (1 - alpha) interval
    # For robustness (no SciPy dependency), approximate inverse CDF via numpy + rational approx
    # but here we use a small lookup-based Newton refinement for accuracy.
    def norm_ppf(p: float) -> float:
        # Wichura-like approximation adapted for typical central probabilities
        # Good enough for alpha in [1e-4, 0.5]
        # For simplicity, use numpy from erf inverse
        from math import sqrt
        return sqrt(2) * erfinv(2 * p - 1)

    z = float(norm_ppf(1.0 - alpha / 2.0))
    if not np.isfinite(z) or z <= 0:
        raise ValueError(f"Computed invalid z for alpha={alpha}. Check alpha value.")

    # Compute half-widths
    hw_upper = y_upper - y_pred   # may be negative if inputs are swapped
    hw_lower = y_pred - y_lower

    # Ensure bounds are ordered (swap if necessary)
    swapped = hw_upper < 0
    if np.any(swapped):
        # Swap to enforce lower <= pred <= upper
        yu = y_upper.copy()
        yl = y_lower.copy()
        y_upper[swapped] = yl[swapped]
        y_lower[swapped] = yu[swapped]
        hw_upper = y_upper - y_pred
        hw_lower = y_pred - y_lower

    # Detect asymmetry per sample (10% relative tolerance)
    rel_diff = np.abs(hw_upper - hw_lower) / np.maximum(1e-12, (hw_upper + hw_lower))
    is_asym = rel_diff > 0.10

    if asymmetry == "symmetric":
        half_width = 0.5 * (y_upper - y_lower)
        sigma_aleatoric = half_width / z
    elif asymmetry == "average":
        sigma_aleatoric = 0.5 * (hw_upper / z + hw_lower / z)
    elif asymmetry == "auto":
        sigma_aleatoric = np.where(
            is_asym,
            0.5 * (hw_upper / z + hw_lower / z),    # average when asymmetric
            0.5 * (y_upper - y_lower) / z,          # symmetric when balanced
        )
    else:
        raise ValueError("asymmetry must be one of {'auto','symmetric','average'}.")

    # Prevent negative/NaN (can occur if bounds are malformed)
    sigma_aleatoric = np.where(np.isfinite(sigma_aleatoric), sigma_aleatoric, np.nan)
    var_aleatoric = sigma_aleatoric ** 2

    # Combine with epistemic via variance addition
    var_total = var_aleatoric + epistemic_var

    if clip_nonnegative:
        var_aleatoric = np.clip(var_aleatoric, 0.0, None)
        var_total = np.clip(var_total, 0.0, None)

    sigma_total = np.sqrt(var_total)

    return sigma_aleatoric, var_aleatoric, var_total, sigma_total

df = pd.read_csv(infile, sep="\t", header=0)
#print(df.columns)  #['Protein_file', 'Ligand_file', 'ES_complex', 'Actual_pKm', 'Epistemic_error', 'MAPIE_mean_pred', 'MAPIE_CI_lower', 'MAPIE_CI_upper']

aleatoric_sd = []
aleatoric_variance = []
total_uncertainty = []
total_uncertainty_sd = []
for i, row in df.iterrows():
	sigma_aleatoric, var_aleatoric, var_total, sigma_total = combine_uncertainties_from_mapie(row["MAPIE_mean_pred"], row["MAPIE_CI_lower"], row["MAPIE_CI_upper"], 0.1, row["Epistemic_error"])
	#print(sigma_aleatoric, var_aleatoric, var_total, sigma_total)
	aleatoric_sd.append(sigma_aleatoric[0])
	aleatoric_variance.append(var_aleatoric[0])
	total_uncertainty.append(var_total[0])
	total_uncertainty_sd.append(sigma_total[0])
	#break

df["Aleatoric_SD"] = aleatoric_sd
df["Aleatoric_var"] = aleatoric_variance
df["Total_uncertainty"] = total_uncertainty
df["Total_uncertainty_SD"] = total_uncertainty_sd

df.to_csv(outfile, sep="\t", header=True, index=False)
print(min(total_uncertainty_sd), max(total_uncertainty_sd), np.mean(total_uncertainty_sd), np.median(total_uncertainty_sd))

























