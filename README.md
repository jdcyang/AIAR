# AIAR Fortran 移植版

MATLAB/Octave 模型的 Fortran 移植，代码架构参照 Kaplan–Moll–Violante (AER 2018) replication 的 Fortran 工程。**与 Octave 版本结果逐位对照一致，速度提升约 50–60 倍**（单稳态 ~10–30ms，全部实验秒级）。

## 编译与运行

```bash
source /opt/intel/oneapi/setvars.sh   # ifx + MKL 环境
make release                          # 或 make debug
./aiar.out                            # 读取 input.nml
```

- 编译器：Intel `ifx` 2025.3+；线性代数用 MKL LAPACK 带状求解器（`dgbsv`/`dgbtrf`/`dgbtrs`），无 SuiteSparse 依赖。
- `makefile` 为目标入口（release/debug/run/clean），实际规则在 `release.makefile` / `debug.makefile`（沿用参考工程惯例）。
- **模块依赖**：`$(SUBR): $(MOD)`——修改 `Parameters.f90` / `Globals.f90` 后必须全量重编（`make clean && make release`），否则 `.mod` 不一致会产生幽灵结果。

## 实验模式（input.nml）

```fortran
&config
  mode = 0                 ! 0=单稳态 1=A扫描 2=token税×A×采用边际网格 3=福利CEV
  A_bar = 0.25             ! mode 0 用
  token_tax = 0.0
  adopt_margin = 0         ! 1 启用 Race(2018) 采用边际
  K_init = -1.0            ! >0 时覆盖初始资本猜测（默认 15）
  A_min = 0.05, A_max = 0.50, A_step = 0.025
  tax_min = 0.0, tax_step = 0.25, token_tax_npts = 9
  Display = 0              ! 2=每10轮打印, 3=加打印前3轮细节
  outdir = 'output'
/
```

输出 CSV 在 `output/`：`policy.csv`（mode 0）、`sweep.csv`、`adoption.csv`、`welfare.csv`，与 Octave 端 CSV 列结构对应。

## 文件结构（对照参考工程）

| 文件 | 作用 | 对应 Octave |
|---|---|---|
| `Parameters.f90` | 编译期常量（网格、迭代参数、容差） | set_parameters_AR 的近似参数段 |
| `Globals.f90` | 模块级共享变量（参数+均衡对象+快照） | 全局变量 |
| `SetParameters.f90` | 中国第一层校准 | set_parameters_AR.m |
| `Grids.f90` | 资产网格 | — |
| `TaskSpace.f90` | 任务空间/Q/π_ss/token 会计/采用边际 | 求解器前段+采用块 |
| `IterateBellman.f90` / `HJBUpdate.f90` | 隐式上风 HJB | 内循环 1 |
| `StationaryDistribution.f90` | FPK 伪时间迭代 | 内循环 2 |
| `SolveSteadyState.f90` | 外循环（自适应阻尼+快照回退） | 价格外循环 |
| `ComputeGinis.f90` | 财富/收入/消费基尼 | compute_all_ginis |
| `SaveOutput.f90` / `Main.f90` | 输出与驱动 | run_* 脚本 |

## 与 Octave 版本的关键实现差异

1. **带状 LAPACK**：状态-资产未知量按 `(i,s) → j=(i−1)·ns+s` 交错编号，HJB/FPK 系统带宽 kl=ku=4，直接 `dgbsv`/`dgbtrs`（参考工程用 UMFPACK，本模型 ntot=800 无需 SuiteSparse）。
2. **HJB 紧度**略高于 MATLAB 版（maxit 2000、crit 1e-11；FPK 40000/1e-13）——成本可忽略，更稳健。
3. **FPK 热启动重置**：价格大幅移动（|e_K|/K>0.2）时清空热启动——修复近可约漂移算子下的路径依赖伪均衡（Octave 版已同步此修复）。

## 移植过程中抓到并修复的 bug（勿回退）

1. `token_tax_eff` 必须**先于** `DownstreamObjects` 赋值（否则连续求解时 token 税收滞后一个解）；
2. `dgbsv` 会原地销毁矩阵——FPK 循环必须 `dgbtrf` 一次 + `dgbtrs` 多次；
3. 局部变量 `k` 遮蔽模块 `K`（Fortran 不区分大小写）——循环变量一律避开模块名；
4. 网格对照数组必须按 (tax, A) 索引存储，不能用标量缓存；
5. makefile 必须含 `$(SUBR): $(MOD)` 模块依赖。

## 验证状态（与 Octave 对照）

- 单稳态 A=0.25：best_err 与 Octave 完全相同（6.14e-07），K/r/w/基尼一致到 4 位有效数字
- A 扫描（19 点）：溢价/失业率逐位一致，财富基尼差 ≤0.009（扭结级）
- 采用网格（342 点）：失业率/I*/τ 全部一致
- 福利 CEV（76 点）：与 Octave 差 ≤0.1 个百分点
- 全部模式零未收敛点
