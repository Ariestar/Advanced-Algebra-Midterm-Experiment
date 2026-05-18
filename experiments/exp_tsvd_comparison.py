import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.config import FIGURES_DIR, RESULTS_DIR, RNG_SEED, ensure_output_dirs
from src.metrics import relative_l2_error
from src.operators import add_relative_noise, make_ill_conditioned_matrix, make_test_signal
from src.solvers import solve_tikhonov_from_svd, solve_truncated_svd


def main() -> None:
    ensure_output_dirs()
    n = 256
    alpha = 8
    rho = 1e-3
    lam = 1e-6
    ranks = [5, 10, 20, 40, 80, 128, 192, 256]
    rng = np.random.default_rng(RNG_SEED + 444)

    x = make_test_signal(n)
    a, _ = make_ill_conditioned_matrix(n, alpha, rng)
    b = a @ x
    b_tilde, _ = add_relative_noise(b, rho, rng)
    u, s, vt = np.linalg.svd(a, full_matrices=False)

    rows = []
    for rank in ranks:
        x_tsvd = solve_truncated_svd(u, s, vt, b_tilde, rank)
        rows.append(
            {
                "method": "truncated_svd",
                "rank": rank,
                "lambda": np.nan,
                "relative_error": relative_l2_error(x_tsvd, x),
            }
        )

    x_tikh = solve_tikhonov_from_svd(u, s, vt, b_tilde, lam)
    rows.append(
        {
            "method": "tikhonov",
            "rank": np.nan,
            "lambda": lam,
            "relative_error": relative_l2_error(x_tikh, x),
        }
    )

    df = pd.DataFrame(rows)
    csv_path = RESULTS_DIR / "exp_tsvd_comparison.csv"
    df.to_csv(csv_path, index=False)

    tsvd_df = df[df["method"] == "truncated_svd"]
    tikh_error = float(df[df["method"] == "tikhonov"]["relative_error"].iloc[0])
    best_tsvd = tsvd_df.loc[tsvd_df["relative_error"].idxmin()]

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(tsvd_df["rank"], tsvd_df["relative_error"], marker="o", label="truncated SVD")
    ax.axhline(tikh_error, color="#F58518", linestyle="--", label=f"Tikhonov $\\lambda={lam:g}$")
    ax.scatter([best_tsvd["rank"]], [best_tsvd["relative_error"]], color="#54A24B", zorder=3, label="best tested TSVD")
    ax.set_yscale("log")
    ax.set_xlabel("truncation rank k")
    ax.set_ylabel("relative L2 error")
    ax.set_title("Truncated SVD and Tikhonov comparison")
    ax.grid(alpha=0.25, which="both")
    ax.legend()
    fig.tight_layout()
    figure_path = FIGURES_DIR / "signal" / "fig_tsvd_comparison.png"
    fig.savefig(figure_path, dpi=220)
    plt.close(fig)

    metadata = {
        "seed": RNG_SEED + 444,
        "n": n,
        "alpha": alpha,
        "condition_number": 10**alpha,
        "rho": rho,
        "lambda": lam,
        "ranks": ranks,
        "best_tsvd_rank": int(best_tsvd["rank"]),
        "best_tsvd_relative_error": float(best_tsvd["relative_error"]),
        "tikhonov_relative_error": tikh_error,
        "csv": str(csv_path),
        "figure": str(figure_path),
    }
    (RESULTS_DIR / "exp_tsvd_comparison.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
