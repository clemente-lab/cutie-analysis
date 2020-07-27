# cutie-analysis #

The repository 'cutie-analysis' contains all relevant scripts for simulating correlation classes as well as automated generation of configuration files, commands, and figures associated with the CUTIE (Correlations Under the Influence) package and publication.

## File structure ###

```
data/
    df_dictionary.docx: description of what column names mean in dataframe .txt files
    sim_results_df.txt: analysis of simulated scatterplots across all parmameters tested
    real_results_df.txt: analysis of real-world datasets across all parmameters tested
plots/
scripts/
    gen_sim.R: generates simulated correlations of NP, FN, FP and CD classes
    analyze_real.py: produces donut-plots of real datasets
    analyze_simulations.py: generates power curves for simulated correlations
    gen_batch.py: generates batch jobs
    gen_commands_configs.py: generates configs for simulated data
    gen_real_commands_configs.py: generates configs for real data
    cutie_figures.R: generates publication-ready plots based on output from above
```
