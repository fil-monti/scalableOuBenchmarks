# Dense OU Package Benchmarks

This folder contains first-pass R benchmark wrappers for dense multivariate OU likelihood or fit surfaces in:

- `PCMBase` / `PCMBaseCpp`
- `mvMORPH`
- `mvSLOUCH`
- `PCMFit`, if installed from GitHub

The scripts generate the same synthetic dense OU tree data for each package grid point. They intentionally report likelihood timings when the package exposes a reusable likelihood evaluator, and fit-time proxy timings when it does not.

## Install

```sh
Rscript benchmarks/benchmarksTrees/install_packages.R
```

`PCMFit` is not available on CRAN in the same way as the other packages. To try the optional GitHub install:

```sh
INSTALL_PCMFIT=true Rscript benchmarks/benchmarksTrees/install_packages.R
```

Packages are installed into `benchmarks/benchmarksTrees/r-lib`.

## Run

Small smoke run:

```sh
K_VALUES=2 N_VALUES=8 REPS=1 WARMUP=0 Rscript benchmarks/benchmarksTrees/run_all.R
```

Larger grid:

```sh
K_VALUES=2,4,8,16 N_VALUES=50,100,250 REPS=10 WARMUP=3 Rscript benchmarks/benchmarksTrees/run_all.R
```

Each runner writes a CSV next to the scripts:

- `results_pcmbase.csv`
- `results_mvmorph.csv`
- `results_mvslouch.csv`
- `results_pcmfit.csv`

## BEAST XML

Generate matching tree-data-likelihood XML benchmark files:

```sh
K_VALUES=2,4 N_VALUES=50,100 BEAST_ITERATIONS=1000 \
Rscript benchmarks/benchmarksTrees/generate_beast_xml.R
```

The generator writes XML files into `benchmarks/benchmarksTrees/beast_xml`.
By default it emits two BEAST cases per grid point:

- `smbp`: canonical `traitDataLikelihood` with `orthogonalBlockDiagonalPolarStableMatrixParameter`
- `dense`: canonical `traitDataLikelihood` with a dense `matrixParameter` selection matrix

Use `BEAST_XML_MODES=smbp` or `BEAST_XML_MODES=dense` to emit only one family. The XML wraps the likelihood in BEAST's `<benchmarker iterationCount="...">`, which repeatedly dirties and evaluates the tree data likelihood.

Run one XML with the local BEAST checkout:

```sh
BEAST_DIR=/Users/filippomonti/Desktop/parallelDiffusions/beast-mcmc-time-series-e0
BEAGLE_DIR=/Users/filippomonti/Desktop/CodingStation/Projects/beagle-lib/build/libhmsbeagle
XML=benchmarks/benchmarksTrees/beast_xml/beast_tree_smbp_k2_n50_seed10068.xml

DYLD_LIBRARY_PATH="$BEAGLE_DIR:$BEAGLE_DIR/CPU" \
java -Djava.library.path="$BEAGLE_DIR/JNI" -Xmx2g \
  -cp "$BEAST_DIR/build:$BEAST_DIR/lib/junit-4.4.jar:$BEAST_DIR/lib/EJML-core-0.30.jar:$BEAST_DIR/lib/EJML-dense64-0.30.jar:$BEAST_DIR/lib/colt.jar:$BEAST_DIR/lib/beagle.jar:$BEAST_DIR/lib/commons-math-2.2.jar:$BEAST_DIR/lib/jebl.jar:$BEAST_DIR/lib/jdom.jar:$BEAST_DIR/lib/jam.jar:$BEAST_DIR/lib/options.jar:$BEAST_DIR/lib/mtj.jar" \
  dr.app.beast.BeastMain "$XML"
```

## Benchmark Surfaces

`PCMBase` reports both a generic dense asymmetric OU `H` model and the default Schur-transformable dense OU model type from `PCMDefaultModelTypes()[["F"]]`. Each is timed through direct `PCMLik` calls and the reusable `PCMCreateLikelihood` closure. The `PCMBaseCpp` rows use `PCMInfoCpp` metadata with the same model and data.

`mvMORPH` uses `mvOU(..., model = "OU1", optimization = "fixed")` and repeatedly calls the returned `fit$llik(par)` closure. The default decomposition is `MVMORPH_DECOMP=schur`, because that is the closest dense stable parametrization target here. Override with `MVMORPH_DECOMP=qr`, `eigen`, etc. if needed.

`mvSLOUCH` does not expose a clean public fixed-parameter likelihood closure in the installed API. The current runner benchmarks `ouchModel` fit time with `Atype=Any` and `Syytype=Any`, and marks rows with `fit_time_proxy_no_public_likelihood_closure`. Tune the short optimizer budget with `MVSLOUCH_MAXITER=2,2`.

`PCMFit` is an inference wrapper around PCMBase. The runner starts from the simulated parameter vector, disables random starts, and uses a one-start optimizer budget. Rows are skipped cleanly when `PCMFit` is not installed.

Garbage collection is outside the timed loop in all runners.
