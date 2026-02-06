# sparsedecomp

(Fast) sparse decomposition code in sage

February 2026

# What's here

*   (Not novel) code to compute a polynomial decomposition using
    standard dense techniques

*   Sampling-based Monte Carlo test for decomposability

*   Experimental code for the Monte Carlo test

# How to install and run

1.  Install [pixi](https://pixi.prefix.dev) (if you don't have it already)

        curl -fsSL https://pixi.sh/install.sh | sh

2.  Get pixi to install various packages (notably [sage](https://www.sagemath.org/))

        pixi install

3.  Run the experimental code to generate fresh data and graphs

        pixi run python experiment_run.py
