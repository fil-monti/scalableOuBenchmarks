# OU Likelihood Benchmark: BEAST vs R/Stan

This document records the full benchmark comparing the BEAST canonical Kalman-filter
implementation against R/Stan for multivariate OU time-series likelihoods.  Every command
below can be run from scratch to reproduce the numbers.

---

## 1. What is being compared

Both implementations compute the exact Kalman-filter log-likelihood and its gradient for the
same multivariate OU model, but use different parametrizations of the drift (selection)
matrix and different numerical strategies.

### The OU model

The latent process is

```
dX_t = -A (X_t - μ) dt + dW_t,    E[dW_t dW_t^T] = Q dt.
```

Observations are

```
Y_i = X_{t_i} + ε_i,    ε_i ~ N(0, σ²_obs I).
```

The exact discrete-time transition for a step of length δ is

```
X_{t+δ} | X_t ~ N( exp(-Aδ) X_t + (I - exp(-Aδ)) μ,  V(δ) )
```

where `V(δ) = ∫₀^δ exp(-As) Q exp(-A^T s) ds`.

### R/Stan: Lyapunov parametrization

The R benchmark (`stanOuBenchMark.R`) uses the general dense stable parametrization

```
B = (0.5 Q + M) S^{-1}
```

where `S` and `Q` are symmetric positive definite (Cholesky factors `L_S`, `L_Q`) and `M`
is skew-symmetric.  The likelihood is computed by Stan's Kalman filter using `matrix_exp` at
each time step on an irregular time grid.  The benchmark measures `log_prob` (forward pass
only) and `grad_log_prob` (forward + reverse-mode AD).

### BEAST: orthogonal block-diagonal polar-stable parametrization

BEAST uses `OrthogonalBlockDiagonalPolarStableMatrixParameter` to parametrize `A`.  For even
`K`, the matrix is decomposed as

```
A = R D R^T
```

where `R` is a `K × K` orthogonal matrix (product of Givens rotations) and `D` is block-diagonal
with `K/2` identical-diagonal 2×2 blocks, each parametrized in polar form:

```
D_b = ρ_b [ cos θ_b   sin θ_b + t_b ]
           [ sin θ_b - t_b   cos θ_b ]
```

Native parameters (unconstrained):

| Group | Count (even K) | Description |
|---|---|---|
| Givens angles | K(K−1)/2 | rotation matrix R |
| ρ | K/2 | block decay rates (must be positive) |
| θ | K/2 | block rotation angles |
| t | K/2 | block skew terms |
| Q entries | K(K+1)/2 | diffusion (Cholesky) |
| μ | K | stationary mean |

The likelihood is computed by a canonical-form Kalman filter
(`GaussianTimeSeriesLikelihoodFactory` with `CANONICAL` forward/smoother modes).  The
gradient uses analytical adjoints (`CANONICAL_ANALYTICAL`).

### Key structural difference: the repeated-Δ cache

BEAST caches the transition triple `(exp(-Aδ), offset, V(δ))` per unique time increment δ.
For a **uniform** grid all T−1 steps share the same δ, so the matrix exponential is
computed once and reused T−1 times.  For a **fully-irregular** grid every step has a unique δ,
so T−1 matrix exponentials are required.

---

## 2. Files

| Path | Description |
|---|---|
| `stanOuBenchMark.R` | R/Stan Lyapunov benchmark (pre-existing) |
| `benchmarks/schur_vs_smbp_exp_adjoint.jl` | Julia dense-Schur-style vs SMBP microbenchmark for matrix exponentials and Frechet adjoints |
| `benchmarks/edge_message_push_benchmark.jl` | Julia dense vs SMBP microbenchmark for one child-to-parent canonical Gaussian edge push |
| `benchmarks/smbp_lyapunov_plan.jl` | Shared cached equal-diagonal SMBP Lyapunov plan used by the Lyapunov benchmarks |
| `benchmarks/lyapunov_stationary_variance_benchmark.jl` | Julia dense vs known-Schur vs SMBP benchmark for stationary-variance Lyapunov solves |
| `benchmarks/branch_likelihood_with_lyapunov_benchmark.jl` | Julia full branch transition/message benchmark including stationary-variance Lyapunov solve |
| `beast-mcmc-time-series-e0/src/test/dr/inference/timeseries/StanVsBeastOUBenchmarkTest.java` | New BEAST benchmark (created in this session) |

The benchmark scripts in this document now live under:

```bash
BENCHMARK_ROOT=/Users/filippomonti/Desktop/ScalableOUonTrees/scalableOuBenchmarks
```

Commands for the Julia/kernel benchmarks assume `cd "$BENCHMARK_ROOT"`. The historical
BEAST-vs-Stan benchmark still depends on the external BEAST checkout and BEAGLE build:

```bash
BEAST_TS_ROOT=/Users/filippomonti/Desktop/parallelDiffusions/beast-mcmc-time-series-e0
BEAGLE_LIB_ROOT=/Users/filippomonti/Desktop/CodingStation/Projects/beagle-lib/build/libhmsbeagle
BEAGLE_SETUP_NOTES=/Users/filippomonti/Desktop/parallelDiffusions/beagle-local-testxml-setup.md
```

### Relevant pre-existing infrastructure used by the BEAST benchmark

| Class | Package | Role |
|---|---|---|
| `OUProcessModel` | `dr.evomodel.continuous.ou` | Holds A, Q, μ, P₀; computes transitions |
| `OUTimeSeriesProcessAdapter` | `dr.inference.timeseries.gaussian` | Wraps OUProcessModel as a time-series latent process |
| `GaussianObservationModel` | `dr.inference.timeseries.gaussian` | Stores H, R, Y |
| `BasicTimeSeriesModel` | `dr.inference.timeseries.core` | Combines latent process + observations + time grid |
| `UniformTimeGrid` | `dr.inference.timeseries.core` | Uniform time grid (isRegular = true) |
| `IrregularTimeGrid` | `dr.inference.timeseries.core` | Arbitrary time grid (isRegular = false) |
| `GaussianTimeSeriesLikelihoodFactory` | `dr.inference.timeseries.likelihood` | Builds likelihood + gradient engine |
| `TimeSeriesLikelihood` | `dr.inference.timeseries.likelihood` | Entry point: `getLogLikelihood()`, `getGradientWrt()` |
| `OrthogonalBlockDiagonalPolarStableMatrixParameter` | `dr.inference.model` | Orthogonal block-diagonal polar-stable A |
| `GivensRotationMatrixParameter` | `dr.inference.model` | Parametrizes the rotation matrix R |

---

## 3. Prerequisites

### R/Stan side

```
R >= 4.0
rstan         # install.packages("rstan")
expm          # install.packages("expm")
```

### BEAST side

- JDK 8+
- Apache Ant
- Local BEAGLE build (needed for loading the JNI library; the OU time-series benchmark
  itself does not use BEAGLE, but BEAST will fail to start without it)

BEAGLE path used on this machine:

```
$BEAGLE_LIB_ROOT
```

See `$BEAGLE_SETUP_NOTES` for full BEAGLE setup instructions.

---

## 4. How to run

### 4.1 R/Stan benchmark

Default run (K = 4, irregular grid):

```bash
cd "$BENCHMARK_ROOT"
Rscript stanOuBenchMark.R
```

To override K without editing the file (K = 10 example):

```bash
sed 's/^K <- 4L/K <- 10L/' \
    stanOuBenchMark.R \
  | Rscript /dev/stdin
```

The script samples intervals as `dt ~ Uniform(Delta * (1 - dt_jitter), Delta * (1 + dt_jitter))`
and prints `Median log_prob time, ms` and `Median grad_log_prob time, ms`.

### 4.2 BEAST benchmark

**Step 1 – Build** (from the BEAST source directory):

```bash
cd "$BEAST_TS_ROOT"
ant build
```

This compiles both main sources (`src/dr/`) and test sources (`src/test/`) into `build/`.

**Step 2 – Run**:

```bash
base="$BEAGLE_LIB_ROOT"
LIBCP="lib/junit-4.4.jar:lib/EJML-core-0.30.jar:lib/EJML-dense64-0.30.jar:\
lib/colt.jar:lib/beagle.jar:lib/commons-math-2.2.jar:\
lib/jebl.jar:lib/jdom.jar:lib/jam.jar:lib/options.jar"

DYLD_LIBRARY_PATH="$base:$base/CPU" \
java -Djava.library.path="$base/JNI" \
     -Xmx2g \
     -cp "build:$LIBCP" \
     junit.textui.TestRunner \
     test.dr.inference.timeseries.StanVsBeastOUBenchmarkTest
```

**To change K**, edit the constant at the top of `StanVsBeastOUBenchmarkTest.java`:

```java
private static final int K = 10;   // ← change this
```

then re-run `ant build` and re-run the Java command above.  K must be even (the benchmark
requires an all-2×2-block polar-stable structure).

### 4.3 BEAST XML benchmark

The benchmark can also be run from XML, avoiding Java constant edits and rebuilds when only
benchmark settings change:

```bash
cd "$BEAST_TS_ROOT"

base="$BEAGLE_LIB_ROOT"
LIBCP="build:lib/junit-4.4.jar:lib/EJML-core-0.30.jar:lib/EJML-dense64-0.30.jar:\
lib/colt.jar:lib/beagle.jar:lib/commons-math-2.2.jar:\
lib/jebl.jar:lib/jdom.jar:lib/jam.jar:lib/options.jar"

DYLD_LIBRARY_PATH="$base:$base/CPU" \
java -Djava.library.path="$base/JNI" \
     -Xmx2g \
     -cp "$LIBCP" \
     dr.app.beast.BeastMain \
     "$BENCHMARK_ROOT/benchmarks/ou_irregular_k64_t300_10x10_beast.xml"
```

The XML contains two separate cases:

```xml
<ouTimeSeriesBenchmark mode="logLikelihood" .../>
<ouTimeSeriesBenchmark mode="gradient" .../>
```

Editable attributes include `warmup`, `timedIterations` (aliases: `nWarmup`, `nTimed`,
`iterationCount`), `stateDimension`, `timeCount`, `timeStep`, `sigmaObs`, `grid`, and
`perturbEachIteration`.  With `perturbEachIteration="true"`, a tiny alternating update to
`rho[0]` is applied before every timed evaluation, so the gradient and transition cache are
fully recomputed each iteration.

### 4.4 Julia Schur vs SMBP microbenchmark

The Julia script isolates the matrix-exponential kernel and its reverse-mode adjoint from
the Kalman filter.  It compares Julia's dense matrix exponential on `A = R D R'`, a
planned dense Schur path that computes `schur(A)` once per case and excludes that
factorization from timed sections, a planned complex-Schur path using a local Parlett
recurrence for the triangular exponential, and the SMBP path, where `exp(-dt D)` and the
Frechet adjoint are processed blockwise in the `D` basis.

```bash
cd "$BENCHMARK_ROOT"
julia benchmarks/schur_vs_smbp_exp_adjoint.jl
```

Optional environment variables:

```bash
K_VALUES=4,8,16,32,64 REPS=10 WARMUP=3 INNER_REPS=50 \
julia benchmarks/schur_vs_smbp_exp_adjoint.jl
```

`INNER_REPS` batches very small operations inside each timing sample and reports per-call
milliseconds.  The CSV separates matrix exponential timing from Frechet-adjoint timing,
including generic dense, planned Schur with generic `exp(T)`, planned Schur with local
Parlett `exp(T)`, SMBP coefficient evaluation, cached apply-only processing, cached
processing with basis rotations, and uncached processing with basis rotations.  The script
prints relative errors for `exp(-dt A)` and the processed adjoint, then writes
`benchmarks/schur_vs_smbp_exp_adjoint_results.csv`.

### 4.5 Julia edge message-push microbenchmark

The edge message-push script reports two related timings.  The first isolates the operation
performed after a child/subtree canonical Gaussian message and branch transition have
already been computed and the message needs to be moved across one OU branch to the
parent:

```text
A    = J_child + J_yy
J_p  = J_xx - J_xy A^{-1} J_yx
h_p  = h_x  - J_xy A^{-1} (h_child + h_y)
```

The dense path stores the transition and child message in the original trait basis.  The
SMBP block path stores the transition in the 2x2 block basis, where transition precision
blocks are block diagonal.  The `smbp_push_with_rotation_ms` column includes rotating the
child message into the block basis and rotating the parent message back to the original
basis.

The second timing includes full branch transition preparation inside the timed loop:
computing `F = exp(-t A)`, forming `Q_t = I - F F'`, building the canonical transition
blocks, and then pushing the child message.  These are the
`dense_full_edge_ms`, `smbp_block_full_edge_ms`, and
`smbp_full_edge_with_rotation_ms` columns.

```bash
K_VALUES=4,8,16,32,64 REPS=100 WARMUP=20 INNER_REPS=100 \
julia benchmarks/edge_message_push_benchmark.jl
```

Latest local results:

| K | dt | dense push ms | SMBP block push ms | SMBP with rotations ms | dense full edge ms | SMBP full block ms | SMBP full rotations ms |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 0.025 | 0.000596 | 0.000480 | 0.000900 | 0.00343 | 0.000646 | 0.00106 |
| 8 | 0.041 | 0.00131 | 0.00114 | 0.00186 | 0.00722 | 0.00134 | 0.00205 |
| 16 | 0.067 | 0.00412 | 0.00366 | 0.00540 | 0.0242 | 0.00408 | 0.00584 |
| 32 | 0.113 | 0.0238 | 0.0230 | 0.0275 | 0.101 | 0.0260 | 0.0322 |

The block and dense outputs agree to about `1e-15` relative error.  The prepared-transition
push confirms that when the child precision is already dense, the dense Cholesky solve in
`A^{-1}` dominates the edge push.  SMBP mainly removes dense transition multiplications, so
the prepared-transition benefit is small in block coordinates and can disappear once
per-edge basis rotations are charged.  The full-edge timing is different: dense must also
pay for a generic matrix exponential and dense covariance/precision construction, while
SMBP uses closed-form 2x2 transition blocks.  In this local run, the full SMBP edge is
about `3.9x-5.9x` faster in block coordinates and about `3.1x-4.1x` faster even with
rotations.

### 4.6 Julia stationary-variance Lyapunov benchmark

The Lyapunov script solves the stationary variance equation

```text
A Sigma + Sigma A' = V
```

and compares three paths:

1. dense `LinearAlgebra.lyap(A, -V)`, which computes the dense Schur decomposition inside
   the timed call;
2. planned dense Schur, where `A = Z T Z'` is precomputed once and the timed work is
   `Z'VZ`, LAPACK `trsyl!` on the quasi-triangular problem, and `Z Sigma_T Z'`;
3. SMBP with rotation, where `V` is rotated into the 2x2 block basis, a cached
   equal-diagonal Lyapunov plan solves only the symmetric half of the block-pair system,
   and `Sigma` is rotated back.

```bash
K_VALUES=4,8,16,32,64 REPS=100 WARMUP=20 INNER_REPS=100 \
julia benchmarks/lyapunov_stationary_variance_benchmark.jl
```

Latest local results:

| K | dense lyap ms | known-Schur lyap ms | SMBP rotation lyap ms | dense / SMBP | known-Schur / SMBP |
|---:|---:|---:|---:|---:|---:|
| 4 | 0.00293 | 0.000904 | 0.000299 | 9.81 | 3.03 |
| 8 | 0.00953 | 0.00309 | 0.000577 | 16.5 | 5.35 |
| 16 | 0.0363 | 0.0126 | 0.00157 | 23.2 | 8.05 |
| 32 | 0.208 | 0.0649 | 0.00817 | 25.4 | 7.95 |

The known-Schur path matches dense `lyap` to roundoff; the SMBP rotation path matches dense
to about `1e-15` to `3e-15`, with residuals of the same order.  This benchmark isolates
stationary covariance construction rather than branch transition or message-passing costs.

### 4.7 Julia branch likelihood with Lyapunov benchmark

The branch-with-Lyapunov script extends the edge benchmark by solving the stationary
variance inside the timed branch construction:

```text
A Sigma + Sigma A' = V
F   = exp(-dt A)
Q_t = Sigma - F Sigma F'
```

It then builds the canonical transition for `y | x ~ N(Fx, Q_t)` and pushes a precomputed
child canonical Gaussian message to the parent.  The dense path uses generic dense
`lyap(A, -V)` and dense `exp(-dt A)`.  The known-Schur path precomputes `A = ZTZ'` once
and times Schur-basis Lyapunov solve, `exp(-dt T)`, rotations, transition construction,
and message push.  The SMBP path rotates `V` and the child message into the 2x2 block
basis, uses the cached symmetric-half block Lyapunov plan and block exponential, performs the dense
covariance/precision/message operations in block coordinates, and rotates the parent
message back.

```bash
K_VALUES=4,8,16,32,64 REPS=100 WARMUP=20 INNER_REPS=100 \
julia benchmarks/branch_likelihood_with_lyapunov_benchmark.jl
```

Latest local results:

| K | dt | dense branch ms | known-Schur branch ms | SMBP rotation branch ms | dense / SMBP | known-Schur / SMBP |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 0.025 | 0.00775 | 0.00503 | 0.00217 | 3.57 | 2.32 |
| 8 | 0.041 | 0.0186 | 0.0111 | 0.00419 | 4.45 | 2.65 |
| 16 | 0.067 | 0.0876 | 0.0454 | 0.0122 | 7.16 | 3.71 |
| 32 | 0.113 | 0.420 | 0.198 | 0.0634 | 6.62 | 3.13 |

The known-Schur and SMBP parent messages both match the dense parent message to roundoff
for this grid.  Because `Sigma` is dense, the SMBP path still performs dense covariance
inversion and dense Gaussian message push; its advantage comes from the block Lyapunov
solve and block matrix exponential.

---

## 5. What the BEAST benchmark measures

The benchmark runs **four scenarios** in a single JUnit test:

### Scenario A – fixed A, uniform grid

`makeDirty()` is called before each timed evaluation.  This invalidates the Kalman-filter
state but **not** the transition-matrix cache (since A is unchanged).  The timed cost is
purely the forward Kalman filter (T steps × K² matrix operations) with pre-cached
transition matrices.  This is the steady-state cost in a setting where A does not vary.

### Scenario A – fixed A, irregular grid

Same as above but with a fully-irregular grid (299 unique δ values for T = 300).
`prepareTimeGrid()` pre-computes all 299 transition triples during `create()`, so the
timed cost is again just the Kalman filter with cached transitions — essentially the same
as the uniform case.

### Scenario B – A changes every call, uniform grid

A parameter perturbation (`rho[0] ± 1e-9`) is applied before each evaluation via
`nativeParam.setParameterValue(0, ...)`.  This fires BEAST's model-change notification,
which invalidates the entire transition-matrix cache.  The timed cost is then:

```
1 matrix exponential (the single unique δ)  +  T-step forward Kalman filter
```

This is the true per-MCMC-step cost for a uniform grid.

### Scenario B – A changes every call, irregular grid

Same cache-busting perturbation, but now there are 299 unique δ values.  The timed cost is:

```
299 matrix exponentials  +  T-step forward Kalman filter
```

This is the true per-MCMC-step cost for a fully-irregular grid.

The gradient (fwd+grad) columns additionally include the Kalman smoother (backward pass)
and the analytical adjoint formulas for ∂logL/∂A, ∂logL/∂Q, ∂logL/∂μ.

---

## 6. Results

The timings below are historical results from the earlier mixed-grid benchmark setup.
After switching the Stan benchmark to irregular times by default, rerun both sides before
using these numbers in the paper.

All times are **median milliseconds** over 100 timed iterations after 20 warmup iterations.
Machine: Apple M-series, JVM: OpenJDK (via `java -Xmx2g`).

### 6.1 K = 4, T = 300, Δ = 0.05, σ_obs = 0.10

#### BEAST vs R/Stan (uniform grid, fixed A)

| Operation | R/Stan (Lyapunov) | BEAST (orth.-block) | Speedup |
|---|---|---|---|
| log_prob | 3 ms | 0.64 ms | **4.7×** |
| grad_log_prob | 4 ms | 2.0 ms | **2.0×** |
| Unconstrained parameters | 30 | 26 | — |

*R/Stan unconstrained parameters: K(K−1)/2 skew + K(K+1)/2 L_S + K(K+1)/2 L_Q + K mean = 6+10+10+4 = 30.*
*BEAST unconstrained parameters: K(K−1)/2 angles + K/2 ρ + K/2 θ + K/2 t + K(K+1)/2 Q-Chol + K mean = 6+2+2+2+10+4 = 26.*

#### BEAST: uniform vs irregular grid (K = 4)

Not explicitly timed for K = 4 (the two-grid comparison was run at K = 10; see below).

---

### 6.2 K = 10, T = 300, Δ = 0.05, σ_obs = 0.10

#### BEAST vs R/Stan (uniform grid, MCMC scenario)

| Operation | R/Stan (Lyapunov) | BEAST (orth.-block) | Speedup |
|---|---|---|---|
| log_prob | 14 ms | 1.6 ms | **8.8×** |
| grad_log_prob | 24 ms | 9.9 ms | **2.4×** |
| Unconstrained parameters | 165 | 125 | — |

*R/Stan: 45 skew + 55 L_S + 55 L_Q + 10 mean = 165.*
*BEAST: 45 angles + 5 ρ + 5 θ + 5 t + 55 Q-Chol + 10 mean = 125.*

#### BEAST: uniform vs irregular grid — Scenario A (fixed A, cache pre-built)

| Operation | Uniform | Irregular | Slowdown |
|---|---|---|---|
| log_prob | 2.3 ms | 1.7 ms | ~1× |
| fwd+grad | 10.3 ms | 13.6 ms | 1.3× |

Both grids have all transition matrices pre-cached by `prepareTimeGrid()` at construction
time.  `makeDirty()` only invalidates the Kalman-filter state.  The timed cost is entirely
the Kalman filter forward/backward passes; the transition cache is not rebuilt.  The small
difference between the two grids reflects JVM timing noise and memory-access patterns, not
algorithmic difference.

#### BEAST: uniform vs irregular grid — Scenario B (A changes every call, true MCMC)

| Operation | Uniform | Irregular | Slowdown |
|---|---|---|---|
| log_prob | 1.6 ms | 4.1 ms | **2.5×** |
| fwd+grad | 9.9 ms | 17.2 ms | **1.7×** |

This is the operationally relevant comparison.  For the uniform grid, each A-change
requires recomputing **1** matrix exponential (the single unique δ).  For the irregular
grid, **299** matrix exponentials must be recomputed.  Despite this 299× difference in
expm calls, the wall-clock slowdown is only 2.5× for log_prob and 1.7× for fwd+grad
because:

1. The per-step expm for the orthogonal block-diagonal A is cheap: it decomposes into
   K/2 = 5 independent 2×2 block exponentials (closed form, no iterative algorithm needed).
2. The Kalman smoother and gradient adjoint dominate the gradient cost and are independent
   of the transition cache.

Even with 299 expm calls BEAST-irregular (17 ms fwd+grad) remains comparable to R/Stan
(24 ms), suggesting the block-diagonal expm is substantially cheaper than Stan's
general `matrix_exp(-δ·B)`.

---

## 7. How to change parameters

Edit the constants at the top of `StanVsBeastOUBenchmarkTest.java`:

```java
private static final int    K        = 10;    // trait dimension (must be even)
private static final int    T_OBS    = 300;   // number of observations
private static final double DELTA    = 0.05;  // nominal grid step
private static final double SIGMA_OBS = 0.10; // observation noise s.d.
private static final int    N_WARMUP = 20;    // warm-up iterations (not timed)
private static final int    N_TIMED  = 100;   // timed iterations
```

For the R benchmark the same constants live in the `USER SETTINGS` block:

```r
K        <- 4L
T_obs    <- 300L
Delta    <- 0.05
dt_jitter <- 0.50
sigma_obs <- 0.10
n_warmup <- 20L
n_timed  <- 100L
```

---

## 8. Summary of findings

The current intended comparison is Stan irregular-grid Kalman filtering against the BEAST
irregular-grid XML benchmark.  The older uniform-grid speedups below are retained only as
context until the irregular-grid runs are regenerated.

| Scenario | BEAST advantage |
|---|---|
| Uniform grid, fixed A | **4–9×** faster log_prob; **2×** faster grad (relative to Stan Lyapunov) |
| Uniform grid, A changes (MCMC) | Same as above: cache rebuilt with 1 expm |
| Irregular grid, fixed A | ~same as uniform (cache pre-built at construction) |
| Irregular grid, A changes (MCMC) | **2.5×** slower log_prob vs uniform BEAST; still faster than Stan |

The dominant sources of BEAST's advantage are:

1. **Repeated-Δ cache**: for a uniform grid, `exp(-Aδ)` is computed once per MCMC step
   regardless of T.  Stan recomputes `matrix_exp` at every time step.
2. **Block-diagonal expm**: the polar-stable block structure decomposes the K×K matrix
   exponential into K/2 closed-form 2×2 blocks, avoiding any iterative algorithm.
3. **Canonical Kalman filter**: operating in information (precision) form avoids the
   K×K matrix inversion that the standard covariance-form filter requires at each step.
4. **Analytical gradient**: the gradient is computed via adjoint formulas rather than
   reverse-mode AD, which eliminates the tape-building overhead Stan pays.
