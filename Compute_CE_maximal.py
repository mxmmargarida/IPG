from Instances import *
from Initial_str_KEG import *
# control time
from time import time
from copy import deepcopy

import numpy as np

import gurobipy as grb

from Maximal_strategy import *

# under unix: to limit time
# import signal
#
# def signal_handler(signum, frame):
#     raise Exception("Timed out!")


###########################################################
# SGM
def IterativeSG_NOT_DFS(G,max_iter,opt_solver=1, S=[]):
    r"""Create instances in a standard format.

    Parameters:
    ----------
    G: Game class (see Instances file)
    max_iter: maximum number of sampled games to be solved
    opt_solver: 0 if cplex is used and 1 otherwise (use gurobi); in the paper it is always 1.
    S: intial set of strategies (optinal)
    Returns:
    -------
    ce: array of probabilities indicating the correlated equilibrium for the strategies in the last sampled game S
    Profits: List of profits for each player under ce
    S: final set of strategies
    count: number of iterations, i.e., sampled games solved
    cpu_time: computational time
    """
    # STEP 0 - INITIALIZATION
    # initialize set of strategies
    if S == []:
        S, U_p, Best_m = InitialStrategies(G,opt_solver)
        #S, U_p, Best_m = InitialStrategiesII(G,opt_solver)
    else:
        U_p, S = IndUtilities(G.m(), G.c(), G.Q(), [[] for _ in range(G.m())], [[] for _ in range(G.m())], S)
        Best_m = CreateModels(G.m(), G.n_I(), G.n_C(), G.n_constr(), G.c(), G.Q(), G.A(), G.b())
    S_new = [[] for p in range(G.m())]
    if [[]] in S:
        print("ERROR: There is a player without feasible strategies")
        return [],[],S,0,0,False
    Numb_stra = [len(S[p]) for p in range(G.m())]
    U_depend = [[{} for k in range(G.m())] for p in range(G.m())]
    # set mandatory strategy in the support
    Numb_stra_S = [0]*G.m()
    # STEP 2 - COMPUTE Correlated EQUILIBRIA OF RESTRICTED GAME
    count = 1
    U_depend = Utilities_Poymatrix(G.m(),G.Q(),U_depend,S_new,S,Numb_stra_S)
    list_best = list(range(G.m()))
    time_aux = time()
    ce = []
    ce_previous = ce[:]
    while True and count <= max_iter and time()-time_aux<=3600:
        print("\n\n Processing node ... ", count)
        print("Computing correlated equilibrium.... \n")
        ce_previous = ce[:]
        #signal.signal(signal.SIGALRM, signal_handler)
        #signal.alarm(3600-int(time()-time_aux))   #  seconds
        try:
            ce, Profits, Profits_pure_strategy = ComputeCE(U_depend,U_p,G.m(),G.n_I(),G.n_C(),Numb_stra,opt_solver)
            ### modify ###
            #return ne,Profits,S,count,time()-time_aux
        #except Exception, msg: python2.7
        except Exception:
            print("Time limit exceeded")
            return ce_previous, [], S,count,time()-time_aux, False
        print("Correlated equilibrium computed sucessfully")
        aux = True # no player has incentive to deviate
        S_new = [[] for p in range(G.m())] # set of new strategies to be considered
        # artificial profile to avoid changing BestReactionGurobi
        #Profile = [np.identity(G.n_I()[p]+G.n_C()[p]) for p in range(G.m())] # this will change nothing as we are multiplying by an identity matrix
        aux_p = 0
        it = np.nditer(np.ones(tuple(Numb_stra)),flags=['multi_index'])
        aux_it = [it.multi_index for _ in it]
        # we just find one violated inequality at time: this might not be the most efficient
        while aux and aux_p<G.m(): # FEED BEST RESPONSES WITH CE solution
            p = list_best[aux_p]
            # for each strategy of player p verify if the correlated equilibria constraint is violated
            for bar_xp in range(Numb_stra[p]):
                # determine temporary c[p] which depends on ce and bar_xp
                coef_p = sum(ce[s] for s in aux_it if s[p]==bar_xp)
                c_tmp = G.c()[p]*coef_p
                # determine temporary Q[p][p] which depends on ce and bar_xp
                #Q_tmp = deepcopy(G.Q()[p])
                Q_tmp = [[] for _ in range(G.m())]
                Q_tmp[p] = G.Q()[p][p]*coef_p
                #Q_tmp[p] = Q_tmp[p]*coef_p
                # determine temporary Q[p][k] which depends on ce and bar_xp
                for k in range(G.m()):
                    if k!=p:
                        Q_tmp[k] = sum(np.dot(S[k][s[k]],G.Q()[p][k])*ce[s] for s in aux_it if s[p]==bar_xp)
                #try:
                s_p, u_max, _ = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],c_tmp,Q_tmp,G.A()[p],G.b()[p],[],p,False,Best_m[p],True)
                s_p = Convert_to_maximal(s_p,G,p)
                #except Exception:
                #    print("Timed out while checking best response")
                #    return ce, [], S, count, time()-time_aux, False
                if Profits_pure_strategy[p][bar_xp]+10**-6<= u_max: # ce constraint is violated
                    aux = False
                    S_new[p].append(s_p)
                    Numb_stra_S = deepcopy(Numb_stra)
                    Numb_stra[p] = Numb_stra[p]+1
                    U_depend = Utilities_Poymatrix(G.m(),G.Q(),U_depend,S,S_new,Numb_stra_S)
                    U_p, S = IndUtilities(G.m(), G.c(), G.Q(), S, U_p, S_new)
                    S_new = [[] for _ in range(G.m())]
                    list_best.append(p)
                    list_best = list_best[:aux_p]+list_best[aux_p+1:]
                    break
            aux_p = aux_p+1
        if aux:
            final_time = time()-time_aux
            # verify it is a NE
            CE_is_NE = Verify_CE_NE(G.m(),ce,Profits,U_depend,U_p,Numb_stra)
            return ce, Profits, S,count,final_time, CE_is_NE
        count = count +1
    if time()-time_aux>3600:
        print("Time Limit Exceeded")
    else:
        print(" Maximum number of iterations was attained")
    return ce_previous, [], S,count,time()-time_aux,False

## verify if CE is NE (or can be transformed in one)
# we can do this by solving the feasibiliy problem associated with ce
# This might not be enough as there can be multiple NE for the same support
def Verify_CE_NE(m,ce,Profits_ce,U_depend,U_p,Numb_stra):
    it = np.nditer(np.ones(tuple(Numb_stra)),flags=['multi_index'])
    aux_it = [it.multi_index for _ in it]
    A_supp = [set([]) for _ in range(m)]
    for s in aux_it:
        if ce[s]>10**-4:
            for p in range(m):
                A_supp[p] = A_supp[p].union(set([s[p]]))
    A_supp = [tuple(a) for a in A_supp]
    ne, Profits = FeasibilityProblem_Gurobi(m,A_supp, U_depend,U_p,Numb_stra,None,Profits_ce)
    if ne == []:
        return False
    else:
        return True
        # dist = 0
        # for s in aux_it:
        #     ne_to_ce = 1
        #     for p in range(m):
        #         ne_to_ce = ne_to_ce *ne[p][s[p]]
        #     dist = dist + abs(ne_to_ce-ce[s])
        # if dist<=10**-4:
        #     return True
        # else:
        #     return False


# def Verify_CE_NE(m,ce,Profits_ce,U_depend,U_p,Numb_stra):
#     it = np.nditer(np.ones(tuple(Numb_stra)),flags=['multi_index'])
#     aux_it = [it.multi_index for _ in it]
#     # each player must be indifferent among her strategies
#     for p in range(m):
#         List_profits_p = []
#         Max_profit_p = Profits_ce[p]
#         for bar_s1 in range(Numb_stra[p]):
#             Profit_bar_s1 = sum(ce[s]*(U_p[bar_s1]+sum(U_depend[p][k][bar_s1,s[k]] for k in range(m) if k!=p)) for s in aux_it if s[p] == bar_s1)
#             List_profits_p.append(Profit_bar_s1)
#             if Profit_bar_s1 >=10**-3 and Profit_bar_s1>= Profits_ce[p]-10**-4 and  Profit_bar_s1<= Profits_ce[p]+10**-4
#
#         Max_profit_p=max(List_profits_p)





#######################################################################################################################

#######################################################
##        COMPUTE INDIVIDUAL PROFITS                 ##
#######################################################

# INPUT
# m = number of players
# c = linear objective function coefficients for each player (list of vectors)
# S = list of strategies for each player
# U_p = list of individual profits for each player
# S_new = new strategies to be added to S and to compute individual profit

# OUTPUT
# U_p = list of players individual profits
# S = new set of strategies
def IndUtilities(m, c, Q, S, U_p, S_new):
    for p in range(m):
        for s in S_new[p]:
            U_p[p].append(float(np.dot(c[p],s)-0.5*np.dot(s,np.dot(Q[p][p],s))))
            S[p].append(s)
    return U_p,S

#######################################################################################################################

#######################################################
##        POLYMATRIX PART OF THE PROFITS             ##
#######################################################

# INPUT
# m = number of players
# Q = bilinear coefficients in the objective function for each player (list of matrices)
# p = player for which we are fixing the strategy
# U_p = list of individual profits for each player
# U_depend = list of the players' profit
# S = strategies of each player (list)
# s = profile of strategies being fixed
# numb = last strategy fixed
# Numb_stra_S = number of strategies in S[p]
# OUTPUT
# U_depend = matrice of utilities (in fact it is a dictionary)
def Utilities_Poymatrix(m,Q,U_depend,S,S_new,Numb_stra_S):
    for p in range(m):
        for k in range(p+1,m):
            for sp in enumerate(S_new[p]):
                for sk in enumerate(S[k]+S_new[k]):
                    U_depend[p][k][(Numb_stra_S[p]+sp[0],sk[0])] = float(np.dot(sk[1],np.dot(Q[p][k],sp[1])))
                    U_depend[k][p][(sk[0],Numb_stra_S[p]+sp[0])] = float(np.dot(sp[1],np.dot(Q[k][p],sk[1])))
        for k in range(p):
            for sp in enumerate(S_new[p]):
                for sk in enumerate(S[k]):
                    U_depend[p][k][(Numb_stra_S[p]+sp[0],sk[0])] = float(np.dot(sk[1],np.dot(Q[p][k],sp[1])))
                    U_depend[k][p][(sk[0],Numb_stra_S[p]+sp[0])] = float(np.dot(sp[1],np.dot(Q[k][p],sk[1])))
    return U_depend

#######################################################################################################################

#######################################################
##        COMPUTE Nash Equilibrium                   ##
#######################################################

# INPUT
# S = set of strategies for each player (list)
# M = (p, numb, sigma)
# Back = computation continues from previous computed equilibrium sigma (if Back = True)
# U_depend = polymatrix
# U_p = individual profits
# m = number of players
# n_I = number of binary variables for each player (list)
# n_C = number of continuous variables for each player (list)
# Numb_stra = size of S; number of strategies available for each player (list)
# opt_solver = 0 then use CLEP, = 1 then use Gurobi
# Supp_Stra = M_pos[2] strategies to consider in the support(new strategies should not be considered = S_new of M_pos)
# OUTPUT
# ne = a Nash equilibrium with strategy S[p][numb] of player p in the support

from itertools import combinations_with_replacement, combinations, product,chain

def ComputeCE(U_depend,U_p,m,n_I,n_C,Numb_stra,opt_solver,m_ce=None):
    m_ce = grb.Model("Correlated Equilibrium")
    m_ce.setParam("Threads",2)
    m_ce.setParam("OutputFlag",False)
    m_ce.ModelSense = -1 # maximize
    m_ce.update()
    # we use the mapping for sigma(player 1 strategy, player 2 strategy, ...
    it = np.nditer(np.ones(tuple(Numb_stra)),flags=['multi_index'])
    aux_it = [it.multi_index for _ in it]
    sigma = {s: m_ce.addVar(lb=0,vtype="C",obj=sum(U_p[p][s[p]]+sum(U_depend[p][k][s[p],s[k]] for k in range(m) if k!=p) for p in range(m))) for s in aux_it}
    m_ce.update()
    # create constraints
    # sigma is a probability distribution
    m_ce.addConstr(sum(sigma.values())==1)
    m_ce.update()
    # correlated equilibria constraints
    for p in range(m):
        for s1 in range(Numb_stra[p]):
            for s2 in range(Numb_stra[p]):
                if s1!=s2:
                    m_ce.addConstr(sum(sigma[s]*(U_p[p][s1] - U_p[p][s2]+sum(U_depend[p][k][s1,s[k]]-U_depend[p][k][s2,s[k]]  for k in range(m) if k!=p)) for s in aux_it if s1==s[p]) >=0)
                    m_ce.update()
    m_ce.optimize()
    ce = np.zeros(tuple(Numb_stra))
    Profits =[0 for p in range(m)]
    Profits_pure_strategy = [[0 for _ in range(Numb_stra[p])] for p in range(m)]
    if m_ce.status not in [3,4]:
        for s in aux_it:
            ce[s] = sigma[s].x
            for p in range(m):
                Profits[p] = Profits[p]+ce[s]*(U_p[p][s[p]]+sum(U_depend[p][k][s[p],s[k]] for k in range(m) if k!=p))
                Profits_pure_strategy[p][s[p]] = Profits_pure_strategy[p][s[p]] + ce[s]*(U_p[p][s[p]]+sum(U_depend[p][k][s[p],s[k]] for k in range(m) if k!=p))
    return ce, Profits,Profits_pure_strategy


def FeasibilityProblem_Gurobi(m,A_supp, U_depend,U_p,Numb_stra,m_p = None,Profits_ce=[]):
    #print "\n\n Solving Problem with Supports: ", A_supp
    if m_p == None:
        # initiate model
        m_p = grb.Model("FeasibilityProblem")
        m_p.setParam("Threads", 2)
        # no pritting of the output
        m_p.setParam( 'OutputFlag', False )
        # set objective function direction
        m_p.ModelSense = -1 # maximize
        m_p.update()
        # probability variables
        sigma = [{sp:m_p.addVar(lb=0,vtype="C",name="sigma_"+str(p)+"_"+str(sp)) for sp in A_supp[p]} for p in range(m)]
        m_p.update()

        ########################################################################################################
        ############# WHEN FEASIBILITY PROBLEM HAS MORE THAN ONE SOLUTION ######################################
        ###### MAXIMIZE THE NUMBER OF VARIABLES WITH POSITIVE PROBABILITY ######################################
        # aux = [m_p.addVar(obj = 1, lb=0,vtype="C",name="aux_"+str(p)) for p in range(m)] # aux <= sigma_p_sp
        # m_p.update()
        # for p, sp in enumerate(A_supp):
        #     for s in sp:
        #         m_p.addConstr(aux[p] <= sigma[p][s])
        #         m_p.update()
        ########################################################################################################
        ########################################################################################################

        # profit variables
        v = [m_p.addVar(lb=-1*grb.GRB.INFINITY,vtype="C",name="v_"+str(p)) for p in range(m)]
        m_p.update()
        for p in range(m):
            m_p.addConstr(grb.quicksum(sigma[p].values())==1)
            m_p.update()
        for p, S_p in enumerate(Numb_stra):
            for sp in range(S_p):
                if sp in A_supp[p]:
                    m_p.addConstr(U_p[p][sp]+grb.quicksum(sigma[k][sk]*U_depend[p][k][(sp,sk)] for k in range(m) if k != p for sk in A_supp[k]) == v[p])
                    m_p.update()
                    # vp must be close to Profit of ce
                    m_p.addConstr(v[p]>=Profits_ce[p]-10**-4)
                    m_p.update()
                    m_p.addConstr(v[p]<=Profits_ce[p]+10**-4)
                    m_p.update()
                else:
                    m_p.addConstr(U_p[p][sp]+grb.quicksum(sigma[k][sk]*U_depend[p][k][(sp,sk)] for k in range(m) if k != p for sk in A_supp[k]) <= v[p])
                    m_p.update()
        #m_p.write("apagar.lp")
        m_p.optimize()
        ne = []
        Profits = []
        #print "Solution status for Feasibility Problem: ", m_p.status
        if m_p.status not in [3,4]:
            ne = [[0 for _ in range(Numb_stra[p])] for p in range(m)]
            for p, sp in enumerate(Numb_stra):
                for j in range(sp):
                    if j in A_supp[p]:
                        ne[p][j]= sigma[p][j].x
                Profits.append(v[p].x)
        return ne, Profits


if __name__ == "__main__":
    # # create normal form game
    # #2 players
    m = 2
    n_I = [2,2]
    n_C = [0,0]
    n_constr = [2,2]
    c = [np.array([0,0]),np.array([0,0])]
    Q = [[np.zeros((2,2)), np.array([[5,6],[2,1]])],[np.array([[5,6],[2,1]]),np.zeros((2,2))]]
    A = [np.array([[1,1],[-1,-1]]),np.array([[1,1],[-1,-1]])]
    b = [np.array([1,-1]),np.array([1,-1])]
    type="normal form game"
    G = Game("empty")
    G.Create(m,n_I,n_C,n_constr,c,Q,A,b,type)
    # Computer one pure Nash equilibrium
    ce_1, Profits_SGM,S,numb_iter,cpu_time_not_dfs,CE_is_NE_1 = IterativeSG_NOT_DFS(G,50)
    #print(CE_is_NE)

    # let's compute the CE that is not NE
    S = [[[1,0],[0,1]],[[1,0],[0,1]]]
    ce_2, Profits_SGM,S,numb_iter,cpu_time_not_dfs,CE_is_NE_2 = IterativeSG_NOT_DFS(G,50,1,S)

    # knapsack example
    # m=2
    # n=20
    # ins=0
    # filename  ="Instances/Knapsack/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".npy"
    # G = Game("empty")
    # G.Read_Game(filename)
    # max_iter=50
    # ce, Profits_SGM,S,numb_iter,cpu_time_not_dfs, CE_is_NE  = IterativeSG_NOT_DFS(G,max_iter)
