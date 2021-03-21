# Evolution Outcomes France COVID19

Code used for the paper:

**Evolution of outcomes for patients hospitalised during the first 9 months of the SARS-CoV-2 pandemic in France: A retrospective national surveillance data analysis**, Noémie Lefrancq, Juliette Paireau, Nathanaël Hozé, Noémie Courtejoie, Yazdan Yazdanpanah, Lila Bouadma, Pierre-Yves Boëlle, Fanny Chereau, Henrik Salje, Simon Cauchemez, _The Lancet Regional Health - Europe_, Volume 5, 2021 (https://doi.org/10.1016/j.lanepe.2021.100087)

Our framework is able to monitor in-hospital mortality in real-time, accounting for delays and censoring. It can be used to:
- assess changes in patient outcomes over the course of an epidemic
- estimate (in real-time) detailed age- and sex- specific estimates of those outcomes

### Codes
- Stan code: >_"Model_indiv_data_changeprobsdelays_plostdischarge_final.stan"_
- R code to run the Stan code and estimate the parameters: "Run_model.R"

### Data
- Example of a simulated dataset needed to run the model: _Individual_trajectories.rds_

