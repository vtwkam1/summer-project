# MSc Epidemiology summer project: How might specific biological characteristics of SARS-CoV-2 variants affect the epidemic trajectory?
A simple two-strain deterministic SEIR model to investigate the impact of introducing a novel variant with shorter pre-infectious period (smaller d<sub>E</sub>, larger &sigma;) or increased transmissibility (higher R<sub>0</sub>) on final epidemic size and composition (in terms of infecting strain), peak size and peak timing.

*(Note scripts not thoroughly commented due to time constraints.)*

## [seir_model.R](/seir_model.R)
- Model specified here, defines compartments, strain parameters and ordinary differential equations (ODE) for transitions between compartments.
- Contains function to implement model, solving ODE using [deSolve package](https://cran.r-project.org/web/packages/deSolve/index.html).
- Generates output file with values in compartments at the specified time intervals.
- For logging and error-checking purposes, generates input file recording model input arguments.
- Plots model output, running [plot_seir_model.R](/plot_seir_model.R).

![Model structure](/model_structure.png)

## [run_model_loop_vaccrate_preinf.R](/run_model_loop_vaccrate_preinf.R)
- Loops through table of combinations for variant parameters and runs [seir_model.R](/seir_model.R) for each set of input arguments.
- Runs [match_growthrate.R](/match_growthrate.R) to calculate parameters for variants with the same initial growth rate but different characteristic.
- Calculates summary measures for each scenario, e.g. cumulative infections, prevalence of cases, daily new cases, peak cases, peak timing.
- Captures final situation at end of epidemic for each scenario and appends to file (`end.csv`).
- Runs [plot_individual_function.R](/plot_individual_function.R) to plot summary measures for each scenario.

## [plot_overall.R](/plot_overall.R)
- Combines results from all model runs and generates summary plots for comparison of different variant characteristics, including:
  - Final epidemic size; final epidemic size by infecting strain (stacked bar chart)
  - Overall number of infections (differs from final epidemic size as some individuals infected twice, first by resident strain then by novel variant - "re-infected")
  - Final number re-infected
  - Peak new cases and peak timing (combined and by strain)
  - Daily new cases (time-series), visualises epidemic trajectory for different scenarios in one graph
- Also generates table of combinations for use in [run_model_loop_vaccrate_preinf.R](/run_model_loop_vaccrate_preinf.R).
- Performs simple checks: 
  - Number of scenarios for each "match set" (each set contains shorter pre-infectious period variant, matched more transmissible variant, different starting vaccine coverages for each one)
  - Correct matching of more transmissible variants
  - Runs which ended prematurely and need re-running with longer run-time

## [plot_seir_model.R](/plot_seir_model.R) 
- Plots number in each compartment over run time.

## [plot_individual_function.R](/plot_individual_function.R)
- Plots summary measures over run time, e.g. percentage population infected with resident strain/novel strain/both strains, number in each compartment as a proportion of the population and daily new cases.

## [match_growthrate.R](/match_growthrate.R)
- For some comparability, relative transmissibility advantage of more transmissible variant determined by matching its growth rate at introduction to growth rate at introduction of variant with shorter pre-infectious period.
- Script contains function to calculate R<sub>0</sub> from growth rate at introduction. Other than the characteristics of interest, all other parameters are assumed the same as the resident strain.

## [open_file.R](/open_file.R)
- Allows interrogation of specific scenarios. Specify variant characteristics and filename to open.
