# Smooth and shape-constrained quantile distributed lag model

## Overview
This repository contains two R scripts for fitting models as described in the simulation section of the associated paper. The scripts handle different types of model fitting and consider the two error term distributions discussed in the paper.


## Scripts

1. **Main-v01.R**: This script is used for fitting unimodal models.
2. **Main-v02.R**: This script is used for fitting concave models.

## Usage

Both scripts accept command-line arguments for specifying the model and the error distributions. The allowed values for these parameters correspond to the models and error distributions mentioned in the paper.

### Command-Line Arguments

- `--model` or `-m`: Specifies the model to run. Allowed values are `A`, `B`, or `C`.
- `--error` or `-e`: Specifies the error term distribution. Allowed values are `normal` or `t`.

### Running the Scripts

To run either script, use the following command format:

```sh
Rscript script_name.R -m model -e error_distribution
```

Replace `script_name.R` with `Main-v01.R` or `Main-v02.R` depending on the model type you want to fit. Replace `model` with `A`, `B`, or `C`, and error_distribution with `normal` or `t`.

### Example
Unimodal Model with Normal Error Distribution:
```sh
Rscript Main-v01.R -m A -e normal
```

Concave Model with t Error Distribution:
```sh
Rscript Main-v02.R -m C -e t
```

## Repeating Simulations
To repeat the simulations as described in the paper, you will need to submit 1-2250 jobs via SLURM. This will allow you to comprehensively repeat the simulations under various configurations.

## Generating Figures
You can use the `Plots.ipynb` to generate the plots similar to those in the paper. However, you will need to adjust the directory names where your results are saved accordingly.