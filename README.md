# cutie-analysis #

The repositoory 'cutie-analysis' contains all relevant scripts for simulating correlation classes as well as automated generation of configuration files, commands, and figures associated with the CUTIE (Correlations Under the Influence) package and publication.

## File structure ###

data/
plots/
scripts/
    new_simulations.R: generates simulated correlations of NP, FN, FP and CD classes
    analyze_simulations_real.py: produces donut-plots of real datasets
    analyze_simulations.py: generates power curves for simulated correlations
    gen_batch.py: generates batch jobs
    gen_commands_configs.py: generates configs for simulated data
    gen_real_commands_configs.py: generates configs for real data

