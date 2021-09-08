regression experiments
======================

Code and data for regression experiments.

The key notebook, and the place you probably want to start if you're interested in calculations that substantiate claims in our article, is [```InterpretTopicResults.ipynb```](https://github.com/tedunderwood/period-cohort/blob/main/regression/InterpretTopicResults.ipynb). This documents the primary experiment we ran, and many of the figures in the article are generated in this notebook.

That notebook draws its data mostly from ```/mainresults```; if you're curious how those results themselves were generated, you can consult ```topicmodel_gridsearch.Rmd,``` which generated them.

Bootstrap resampling was carried out by ```topicmodel_gridsearch_bootstrap.Rmd```. Data is in ```/bootstrapdata```, and is interpreted in [```BootstrapDeltaEstimate.ipynb```](https://github.com/tedunderwood/period-cohort/blob/main/regression/BootstrapDeltaEstimate.ipynb).

