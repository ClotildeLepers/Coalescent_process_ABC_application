# Coalescent_process_ABC_application

We propose here ABC inference in population models with phylogenetical information, selection, varying population size and evolving trait structure.

The Code of the repository has been used in the paper "Inference with selection, varying population size and evolving population structure: Application of ABC to a forward-backward coalescent process with interactions", by Clotilde Lepers, Sylvain Billiard, Matthieu Porte, Sylvie Méléard, Viet Chi TRAN. 
(https://arxiv.org/abs/1910.10201) To appear in Heredity. 

## Dieckmann-Doebeli model for a population structured by a trait in [-1,1]

The files 1trait_model.zip, 1trait_simuls_A_D.rar, analysis_1trait.R, creation_kingman.py, normality_Kingman.R and Neutrality_popeffect.R correspond to application 1 presented in section 3.2 of the aforementioned paper.

### Simulation of the model

The 1trait_model.zip file allows to simulate the descriptive statistics obtained using a simple Dieckmann-Doebeli model. 

### Estimation of the parameters of the model

The 1trait_simuls_A_D.rar file gives the result obtained for four sets of parameters (pseudo-data, see App. Table 1). 
The analysis_1trait.R file is used to perform abc analyzes in order to infer the parameters used to simulate the pseudo-data (section 3.2.1).

### Neutrality tests for the phylogenies simulated under our model

The creation_kingman.py file is used to simulate descriptive statistics obtained by a Kingman's coalescent. The normality_Kingman.R and Neutrality_popeffect.R files are used to analyze discrepancy with Kingman's coalescent (section 3.2.2)

## A simple model for Central-Asian Human populations structured by 2 traits

The files 2traits_model.zip, data_central_Asia.csv and analysis_2traits.R, correspond to application 2 presented in section 3.3.

The 2traits_model.zip file allows to simulate the descriptive statistics obtained using an eco-evolutionary model taking into account geographic location and social organization of the population. 

The data_central_Asia.csv file gives the descriptive statistics of Central Asia data given in Chaix et al (2007) and Heyer et al (2015) studies. 

The file analysis_2traits.R is used to perform abc analyzes in order to check the quality of the ABC estimation of parameters, and to infer the parameters leading to the Central Asia data.

The archive CentralAsiaHumanPop_data.zip contains: 1) the summary data from the Central Asia Human Population (from Chaix et al. 2007, Heyer et al. 2015) 2) summary statistics of simulations under our model, for parameters drawn in the a priori distribution. Both files are necessary to run the ABC estimation methods presented in Section 3.3 of our paper.


### Ref: 
Chaix R, Quintana-Murci L, Hegay T, Hammer M, Mobasher Z, Austerlitz F, et al. (2007).  From social to genetic structures in Central Asia. Current Biology 17: 43–48
Heyer E, Brandenburg JT, Leonardi M, Toupance B, Balaresque P, Hegay T, et al. (2015). Patrilineal populations show more male transmission of reproductive success than cognatic populations in Central Asia, which reduces their genetic diversity. American Journal of Physical Anthropology 157: 537–543.
Lepers C, Billiard S, Porte M, Méléard S, Tran VC (2020). Inference with selection, varying population size and evolving population structure: Application of ABC to a forward-backward coalescent process with interactions. To appear in Heredity.
Clotilde Lepers, Sylvain Billiard, Matthieu Porte, Sylvie Méléard, Viet Chi Tran
