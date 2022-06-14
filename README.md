# Equilibria for integer programming games

Paper associated with this algorithmic implementation:

[*M. Carvalho, A. Lodi, J. P. Pedroso, "Computing equilibria for integer programming games". 2022. European Journal of Operational Research*](https://www.sciencedirect.com/science/article/pii/S0377221722002727)

[*M. Carvalho, A. Lodi, J. P. Pedroso, "Computing Nash equilibria for integer programming games". 2020. arXiv:2012.07082*](https://arxiv.org/abs/2012.07082)

## Instances representation ##
Each player p in the set M must solve a problem of the form

![equation](https://latex.codecogs.com/gif.latex?%5Cmax%20%5C%20%5C%20c%5ETx%5Ep%20-%5Cfrac%7B1%7D%7B2%7D%28x%5Ep%29%5ETQ_p%5Epx%5Ep%20&plus;%20%5Csum_%7Bk%20%5Cin%20M%3A%20k%20%5Cneq%20p%7D%20%28x%5Ek%29%5ETQ_k%5Epx%5Ep%5C%5C%20s.t.%20%5C%20%5C%20A%5Epx%5Ep%20%5Cleq%20b%5Ep%20%5C%5C%20x_i%5Ep%20%5Cin%20%5C%7B0%2C1%5C%7D%2C%20i%3D1%2C...%2CB_p)

where

![equation](https://latex.codecogs.com/gif.latex?A%5Ep%20%5Cin%20M_%7Br_p%20%5Ctimes%20n_p%7D%2C%20n_p%20%5Cgeq%20B_p%2C%20b%5Ep%20%5Cin%20M_%7Br_p%20%5Ctimes%201%7D)

**Remark**: The implemented algorithmic approach is guaranteed to return an equilibrium if each player strategy set is bounded and non-empty. Otherwise, the algorithm may fail to stop. In particular, note that this methodology is proven to be correct (see associated paper) **even if there are continous variables**, i.e., it is not necessary for all variables to be integer.

## Description of the files

**Instances.py**: contains all methods to generate random Knapsack Game, Kidney Exchange Game and Lot-Sizing Game instances, as well as the implementations to save and read such instances. In fact, this file creates the class of Integer Programming Games (IPGs). 

**Computer_NE.py**: contains all methods to compute Nash equilibria for an IPG. *IterativeSG* is the method m-SGM wiht time limit of 1h and as input a maximum number of iterations, and *IterativeSG_NOT_DFS* is SGM.

**Computer_CE.py**: contains all methods to compute correlated equilibria for an IPG. *IterativeSG_NOT_DFS* is the method SGM wiht time limit of 1h and as input a maximum number of iterations. Note that I use the same name as the function in Compute_NE, hence these two functions must be loaded with different names.

**Initial_str.py**: contains  three different methods for initiating the sampled game: monopoly strategy (oppenents do nothing), social optimum strategies and potential function strategies (if the game is a potential game, this is a pure Nash equilibrium). It also has the functions for computing a player best reaction/response.

**Run_Generate_instances.py**: file runned to generate all instances for the paper. If you run it, it will overwrite the Instances file.

**MaxPotential_LS.py**: contains the function that computes the pure NE associated with the lot sizing game when we maximize the associated potential function.

**Compute_NE_maximal.py** and **Compute_CE_maximal.py**: this is for KEG where we know we can restrict strategies to be maximal.

Our methodos do not account for the time to build the best response model of each player. This is only significant to the kidney exchange game where the game is quite large (and the reason why there is the function **Initial_str_KEG** which speeds up model construction).

**All_Strategies.py**: file to compare m-SGM with PNS (last experiment in the arxiv paper).


We are using gurobi 9.0.0 and python 3.8.3

## Folders

**Instances/**: it contains a folder with the Knapsack Game and Lot Sizing Game instances (npy files). For the Kidney Exchange Game you must download the KEP instances in https://rdm.inesctec.pt/dataset/ii-2019-001 and run G_KEG = Game("KEG",2,n,ins,50,3) that will transform instance ins of size n of the KEP files into a game.

**Results/**: it contains the files with granular results for the Knapsack Game, Kidney Exchange Game and Lot Sizing Game (txt files). The numbers in the files correspond to their order of apperance in the results tables in the paper.

**RunFiles/**: it contains the code used to generate the result tables in the paper (they must be placed in the main folder otherwise python will not find the code to run them).

