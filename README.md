# Cognitive maps in spatial and non-spatial domains: Similar but different

This git repo contains the data, code, and analyses used in Wu, Schulz, Garvert, Meder & Schuck (20202).

## Data
- The data is located in `experimentData/full.csv`
- The script `dataProcess.R` provides functions for importing data in well formatted dataframes

## Experiment
- `experiment/gaborPatches.R` contains the code used to generate the gabor stimuli
- `experiment/roughEnvironment.json` and `experiment/smoothEnvironment.json` contain the underlying payoff distributions used in the experiment

## Analyses
- All analyses are described in one of three R notebooks, which provide a step-by-step guide of the functions, statistics, and plots used in the final paper
- The `.Rmd` files are located in `docs/` and are separated into [behavioral analyses](https://github.com/charleywu/cognitivemaps/blob/master/docs/behavioralResultsNotebook.Rmd), [model results](https://github.com/charleywu/cognitivemaps/blob/master/docs/modelingResultsNotebook.Rmd), and [bonus round analyses](https://github.com/charleywu/cognitivemaps/blob/master/docs/bonusRoundNotebook.Rmd)
- Each of these R notebooks produces easy to read HTML files that are an accessible way to look at the precise analyses used in each statistical comparison or plot. So please check them out!
  - [behavioral results notebook](https://charleywu.github.io/cognitivemaps/behavioralResultsNotebook.nb.html)
  - [model results notebook](https://charleywu.github.io/cognitivemaps/modelingResultsNotebook.html)
  - [bonus round notebook](https://charleywu.github.io/cognitivemaps/bonusRoundNotebook.html)
   
There are also some helper files that are loaded in these notebooks or are explicitly mentioned (e.g., `modelComparisonCV.R` is used to generate the model estimates)

This is the first time I've ever invested the effort to put together interpretable R notebooks of the analyses in a paper. It wasn't a small amount of work, but I think it could be very useful for providing transparency to the scientific process. Let me know what you think!
