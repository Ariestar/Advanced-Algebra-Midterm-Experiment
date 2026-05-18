# 病态线性算子与 Tikhonov 正则化

本仓库用于高等代数小组实验“病态算子构造与正则化”。

核心目标：

- 构造条件数极大的病态线性算子；
- 展示伪逆求解中的噪声放大；
- 使用 Tikhonov 正则化稳定恢复；
- 用奇异值分解解释实验现象。

## 运行

```powershell
uv sync
uv run python -m experiments.exp_filter_factors
uv run python -m experiments.exp_1d_condition_number
uv run python -m experiments.exp_1d_lambda_selection
uv run python -m experiments.exp_noise_direction
uv run python -m experiments.exp_1d_noise_grid
uv run python -m experiments.exp_tsvd_comparison
uv run python -m experiments.exp_2d_image_recovery
```

生成结果位于：

- `figures/`
- `results/`

研究计划见 [docs/research_plan.md](docs/research_plan.md)。
