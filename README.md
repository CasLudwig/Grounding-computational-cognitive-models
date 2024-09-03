# Grounding computational cognitive models

Repository accompanying [Ludwig, Stuchlý, & Malhotra](https://osf.io/preprints/psyarxiv/vur6t). 

- Code for generating Figure 4 (mixture of Gaussians example).
- Code for generating Figures 5 and C1 (simulated expanded judgement paradigm).

## Mixture of Gaussians

The [mixGaussiansPF](https://github.com/CasLudwig/Grounding-computational-cognitive-models/tree/main/mixGaussiansPF) folder contains all the code for generating Figure 4 in the paper. Assuming you have all the relevant libraries and Stan working together with R, you *should* simply be able to run [timeVarNormalSimulation.R](https://github.com/CasLudwig/Grounding-computational-cognitive-models/blob/main/mixGaussiansPF/timeVarNormalSimulation.R). 

## Simulated expanded judgement paradigm

The [expandedJudgementSimulations](https://github.com/CasLudwig/Grounding-computational-cognitive-models/tree/main/expandedJudgementSimulations) folder contains code and data for generating Figures 5 and C1 (Appendix C) in the paper. Running this code is a bit more involved, but we have tried to make things easier by providing a [RNotebook]() with extensive commentary. You do not need to run the code in this notebook, but it provides insight in how the (simulated) data reported in the paper were generated. The notebook contains all the instructions for recreating the figures, but in brief:

- Figure 5: run [trajectoryExpJudgement.R](https://github.com/CasLudwig/Grounding-computational-cognitive-models/blob/main/expandedJudgementSimulations/trajectoryExpJudgement.R).
- Figure C1: run [rewardRatesExpJudgement.R](https://github.com/CasLudwig/Grounding-computational-cognitive-models/blob/main/expandedJudgementSimulations/rewardRatesExpJudgement.R).

## Disclaimer

If the code does not work because your setup is different from mine (see below), there is not much I can do about this. If the code does not work properly because there are errors, please let me know and I will of course fix as soon as I can.

## Environment

Code was written in the following environment.

