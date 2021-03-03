# cutie-analysis #

The repository 'cutie-analysis' contains all relevant scripts for simulating correlation classes as well as automated generation of configuration files, commands, and figures associated with the CUTIE (Correlations Under the Influence) package and publication.

## File structure ###

```
data/: raw data for each figure as well as spreadsheet summarizing all analyses
scripts/
    gen_sim.R: generates simulated correlations of NP (TP and TN), FN, and FP classes
    analyze_real.py: produces figures for real datasets
    analyze_simulations.py: generates power curves for simulated correlations
    gen_batch.py: generates batch jobs
    gen_commands_configs.py: generates configs for simulated data
    gen_real_commands_configs.py: generates configs for real data
    cutie_figures.R: generates publication-ready plots based on output from above
```
