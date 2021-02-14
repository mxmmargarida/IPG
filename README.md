# Equilibria for integer programming games

Each player p must solve a problem of the form
![equation](http://latex.codecogs.com/gif.latex?O_t%3D%5Ctext%20%7B%20Onset%20event%20at%20time%20bin%20%7D%20t)
$$
max c*x[p] -(1/2)x[p]Q[p][p]x[p] + sum_{j \neq i} x[j]Q[p][j]x[p]
s.t. Ax <=b
x[:n_I] binary and x[:n_C] >=0 continuous
$$

Implementation of the methodology described in 
M. Carvalho, A. Lodi, J. P. Pedroso, "Computing Nash equilibria for integer programming games". arxiv paper

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

