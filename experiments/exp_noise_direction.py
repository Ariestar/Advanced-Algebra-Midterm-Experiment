import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.config import FIGURES_DIR, RESULTS_DIR, RNG_SEED, ensure_output_dirs
from src.metrics import relative_l2_error
from src.operators import make_ill_conditioned_matrix, make_test_signal
from src.solvers import solve_tikhonov_from_svd


def main() -> None:
    ensure_output_dirs()
    n = 256
    alpha = 8
    rho = 1e-3
    lam = 1e-6
    rng = np.random.default_rng(RNG_SEED + 333)

    x = make_test_signal(n)
    a, s = make_ill_conditioned_matrix(n, alpha, rng)
    b = a @ x
    u, s_svd, vt = np.linalg.svd(a, full_matrices=False)

    directions = [
        ("largest singular direction", 0),
        ("middle singular direction", n // 2),
        ("smallest singular direction", n - 1),
    ]

    rows = []
    for label, index in directions:
        noise = rho * np.linalg.norm(b) * u[:, index]
        b_tilde = b + noise

        x_pinv = np.linalg.pinv(a, rcond=1e-15) @ b_tilde
        x_tikh = solve_tikhonov_from_svd(u, s_svd, vt, b_tilde, lam)

        pinv_noise_error = np.linalg.pinv(a, rcond=1e-15) @ noise
        predicted_pinv_noise_error_norm = np.linalg.norm(noise) / s_svd[index]

        rows.append(
            {
                "direction": label,
                "index": index + 1,
                "singular_value": s_svd[index],
                "inverse_factor": 1 / s_svd[index],
                "tikhonov_amplification_factor": s_svd[index] / (s_svd[index] ** 2 + lam),
                "pinv_relative_error": relative_l2_error(x_pinv, x),
                "tikhonov_relative_error": relative_l2_error(x_tikh, x),
                "pinv_noise_error_norm": float(np.linalg.norm(pinv_noise_error)),
                "predicted_pinv_noise_error_norm": float(predicted_pinv_noise_error_norm),
            }
        )

    df = pd.DataFrame(rows)
    csv_path = RESULTS_DIR / "exp_noise_direction.csv"
    df.to_csv(csv_path, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    x_pos = np.arange(len(df))
    labels = ["$u_1$", "$u_{n/2}$", "$u_n$"]

    axes[0].bar(x_pos, df["inverse_factor"], color="#4C78A8")
    axes[0].set_yscale("log")
    axes[0].set_xticks(x_pos, labels)
    axes[0].set_ylabel(r"pseudoinverse factor $1/\sigma_i$")
    axes[0].set_title("Noise amplification predicted by SVD")
    axes[0].grid(alpha=0.25, axis="y", which="both")

    width = 0.38
    axes[1].bar(x_pos - width / 2, df["pinv_relative_error"], width, label="pseudoinverse")
    axes[1].bar(x_pos + width / 2, df["tikhonov_relative_error"], width, label="Tikhonov")
    axes[1].set_yscale("log")
    axes[1].set_xticks(x_pos, labels)
    axes[1].set_ylabel("relative L2 error")
    axes[1].set_title("Recovery error for equal-norm directional noise")
    axes[1].legend()
    axes[1].grid(alpha=0.25, axis="y", which="both")

    fig.tight_layout()
    figure_path = FIGURES_DIR / "signal" / "fig_noise_direction_svd.png"
    fig.savefig(figure_path, dpi=220)
    plt.close(fig)

    metadata = {
        "seed": RNG_SEED + 333,
        "n": n,
        "alpha": alpha,
        "condition_number": 10**alpha,
        "rho": rho,
        "lambda": lam,
        "csv": str(csv_path),
        "figure": str(figure_path),
    }
    (RESULTS_DIR / "exp_noise_direction.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
