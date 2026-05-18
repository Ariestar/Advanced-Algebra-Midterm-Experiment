import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.config import FIGURES_DIR, RESULTS_DIR, RNG_SEED, ensure_output_dirs
from src.metrics import amplification_factor, relative_l2_error
from src.operators import add_relative_noise, make_ill_conditioned_matrix, make_test_signal
from src.solvers import solve_pinv, solve_tikhonov_from_svd


def run_trial(n: int, alpha: float, rho: float, lam: float, seed: int) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    x = make_test_signal(n)
    a, _ = make_ill_conditioned_matrix(n, alpha, rng)
    b = a @ x
    b_tilde, noise = add_relative_noise(b, rho, rng)
    u, s, vt = np.linalg.svd(a, full_matrices=False)

    x_pinv = solve_pinv(a, b_tilde)
    x_tikh = solve_tikhonov_from_svd(u, s, vt, b_tilde, lam)

    return {
        "pinv_error": relative_l2_error(x_pinv, x),
        "tikhonov_error": relative_l2_error(x_tikh, x),
        "pinv_amplification": amplification_factor(x_pinv, x, noise, b),
        "tikhonov_amplification": amplification_factor(x_tikh, x, noise, b),
    }


def plot_error_ratio_heatmap(path, summary: pd.DataFrame, rhos: list[float], kappas: list[float]) -> None:
    ratio = np.empty((len(rhos), len(kappas)))
    for i, rho in enumerate(rhos):
        for j, kappa in enumerate(kappas):
            row = summary[(summary["rho"] == rho) & (summary["condition_number"] == kappa)].iloc[0]
            ratio[i, j] = row["pinv_error_mean"] / row["tikhonov_error_mean"]

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    image = ax.imshow(np.log10(ratio), cmap="viridis", aspect="auto")
    ax.set_xticks(np.arange(len(kappas)), [f"$10^{int(np.log10(k))}$" for k in kappas])
    ax.set_yticks(np.arange(len(rhos)), [f"$10^{{{int(np.log10(r))}}}$" for r in rhos])
    ax.set_xlabel(r"condition number $\kappa(A)$")
    ax.set_ylabel(r"relative noise level $\rho$")
    ax.set_title(r"log$_{10}$(pseudoinverse error / Tikhonov error)")
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label(r"log$_{10}$ error ratio")

    for i in range(len(rhos)):
        for j in range(len(kappas)):
            ax.text(j, i, f"{ratio[i, j]:.1f}x", ha="center", va="center", color="white", fontsize=8)

    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    ensure_output_dirs()
    n = 256
    lam = 1e-6
    alphas = [2, 4, 6, 8]
    rhos = [1e-4, 1e-3, 1e-2]
    repeats = 20

    rows = []
    for alpha in alphas:
        for rho in rhos:
            for trial in range(repeats):
                seed = RNG_SEED + 3000 * int(alpha) + 100 * int(round(-np.log10(rho))) + trial
                result = run_trial(n, alpha, rho, lam, seed)
                rows.append(
                    {
                        "alpha": alpha,
                        "condition_number": 10**alpha,
                        "rho": rho,
                        "trial": trial,
                        "lambda": lam,
                        **result,
                    }
                )

    df = pd.DataFrame(rows)
    csv_path = RESULTS_DIR / "exp_1d_noise_grid.csv"
    df.to_csv(csv_path, index=False)

    summary = (
        df.groupby(["rho", "condition_number"])
        .agg(
            pinv_error_mean=("pinv_error", "mean"),
            pinv_error_std=("pinv_error", "std"),
            tikhonov_error_mean=("tikhonov_error", "mean"),
            tikhonov_error_std=("tikhonov_error", "std"),
            pinv_amplification_mean=("pinv_amplification", "mean"),
            tikhonov_amplification_mean=("tikhonov_amplification", "mean"),
        )
        .reset_index()
    )
    summary_path = RESULTS_DIR / "exp_1d_noise_grid_summary.csv"
    summary.to_csv(summary_path, index=False)

    plot_error_ratio_heatmap(
        FIGURES_DIR / "signal" / "fig_1d_noise_kappa_heatmap.png",
        summary,
        rhos,
        [10**alpha for alpha in alphas],
    )

    metadata = {
        "seed": RNG_SEED,
        "n": n,
        "lambda": lam,
        "alphas": alphas,
        "rhos": rhos,
        "repeats": repeats,
        "csv": str(csv_path),
        "summary_csv": str(summary_path),
        "figure": str(FIGURES_DIR / "signal" / "fig_1d_noise_kappa_heatmap.png"),
    }
    (RESULTS_DIR / "exp_1d_noise_grid.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
