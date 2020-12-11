from Instances import *
from Initial_str import *
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

# m-SGM
def IterativeSG(G,max_iter,opt_solver=1, S=[]):
    r"""Create instances in a standard format.

    Parameters:
    ----------
    G: Game class (see Instances file)
    max_iter: maximum number of sampled games to be solved
    opt_solver: 0 if cplex is used and 1 otherwise (use gurobi); in the paper it is always 1.
    S: intial set of strategies (optinal)
    Returns:
    -------
    ne: list of vectors indicating the Nash equilibrium associated with the strategies in S
    Profits: List of profits for each player under ne
    S: final set of strategies
    count: number of iterations, i.e., sampled games solved
    numb_back: number of backtracking
    cpu_time: computational time
    """
    # create log file
    # f = open("log_ComputeOneNE.txt",'w')
    # f.close()
    # STEP 0 - INITIALIZATION
    # initialize set of strategies
    if S == []:
        S, U_p, Best_m = InitialStrategies(G,opt_solver)
        #S, U_p, Best_m = InitialStrategiesII(G,opt_solver)
    else:
        U_p, S = IndUtilities(G.m(), G.c(), G.Q(), [[] for _ in range(G.m())], [[] for _ in range(G.m())], S)
        Best_m = CreateModels(G.m(), G.n_I(), G.n_C(), G.n_constr(), G.c(), G.Q(), G.A(), G.b())
    S_new = [[] for p in range(G.m())]
    if  [[]] in S:
        print("ERROR: There is a player without feasible strategies")
        return [],[],S,0,0,0
    Numb_stra = [1 for p in range(G.m())]
    U_depend = [[{} for k in range(G.m())] for p in range(G.m())]
    # set mandatory strategy in the support
    Numb_stra_S = [0]*G.m()
    M_pre = [G.m()-1,0,deepcopy(Numb_stra_S),deepcopy(S_new),[[]]*G.m(),deepcopy(U_depend),deepcopy(S_new),[0,0],0] # strategy S[0] must be in the supp of player m
    # set memory
    Memory = [deepcopy(M_pre)]
    Back = False # start computation of equilibria from the start
    numb_back = 0
    # STEP 2 - COMPUTE EQUILIBRIA OF RESTRICTED (SAMPLED) GAME
    # compute Nash equilibrium taking account M and Back
    count = 1
    U_depend = Utilities_Poymatrix(G.m(),G.Q(),U_depend,S_new,S,Numb_stra_S)
    # M_pos = [player p, best response of p, Numb_stra of current node, strategies S of current node, U_p of S, U_depend of S, new best responses from S]
    M_pos = [None,None,deepcopy(Numb_stra),deepcopy(S),deepcopy(U_p),deepcopy(U_depend),deepcopy(S_new),None,0]
    list_best = list(range(G.m()))
    time_aux = time()
    ne = []
    ne_previous = ne[:]
    while True and count <= max_iter and time()-time_aux<=3600:
        # f = open("log_ComputeOneNE.txt",'a')
        # f.write("\n\n Processing node ... "+str(count))
        # f.write("\nComputing equilibria of game of size "+str(Numb_stra)+" .... \n")
        # f.close()
        print("\n\n Processing node ... ", count)
        print("Computing equilibria.... \n")
        ne_previous = ne[:]
        # under unix: to limit time
        #signal.signal(signal.SIGALRM, signal_handler)
        #signal.alarm(3600-int(time()-time_aux))   #  seconds
        try:
            ne, Profits,S0_j = ComputeNE(M_pre,Back,U_depend, U_p,G.m(),G.n_I(), G.n_C(),Numb_stra,opt_solver,M_pos[2],M_pos[-1])
        #except Exception, msg: python 2.7
        except Exception:
            print("Time limit exceeded")
            return ne_previous, [], S,count,numb_back,time()-time_aux
        if ne == []: # fail to compute such equilibrium
            # f = open("log_ComputeOneNE.txt",'a')
            # f.write("\nFail to compute equilibrium\n#############################\n Backtracking ....\n#############################\n")
            # f.close()
            print(" Fail to compute equilibrium")
            print("#############################")
            print("\n Backtracking ....\n")
            print("#############################")
            S_new = deepcopy(M_pre[-3])
            U_p = deepcopy(M_pre[4])
            S = deepcopy(M_pre[3])
            U_depend = deepcopy(M_pre[5])
            Numb_stra_S = deepcopy(M_pre[2])
            U_depend = Utilities_Poymatrix(G.m(),G.Q(),U_depend,S,S_new,Numb_stra_S)
            U_p, S = IndUtilities(G.m(), G.c(), G.Q(), S, U_p, S_new)
            Numb_stra = [Numb_stra_S[p]+len(S_new[p]) for p in range(G.m())] # update number of strategies
            M_pos = deepcopy(M_pre)
            Memory.pop()
            M_pre = deepcopy(Memory[-1])
            Back = True
            numb_back = numb_back+1
            #count = count -1
        else:
            # f = open("log_ComputeOneNE.txt",'a')
            # f.write("\n Equilibrium computed sucessfully\n")
            # f.close()
            print(" Equilibrium computed sucessfully")
            #M_pre[2] = ne
            #Memory[-1] = M
            Back = False # ne found
            aux = True # no player has incentive to deviate
            S_new = [[] for p in range(G.m())] # set of new strategies to be considered
            Profile = [np.array([sum(ne[int(k+sum(Numb_stra[:p]))]*S[p][k][i] for k in range(Numb_stra[p]))  for i in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())]
            aux_p = 0
            while aux and aux_p<G.m(): # FEED BEST RESPONSES WITH NE solution
                p = list_best[aux_p]
                s_p, u_max, _ = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,Best_m[p])
                s_p = Convert_to_maximal(s_p,G,p)
                if Profits[p] +10**-6 <= u_max:
                    aux = False
                    S_new[p].append(s_p)
                    Numb_stra_S = deepcopy(Numb_stra)
                    # update M_pos
                    M_pos[0] = p
                    M_pos[1] = Numb_stra[p]
                    M_pos[-2] = deepcopy(ne)
                    M_pos[-3][p].append(s_p)
                    M_pos[-1] = S0_j
                    Memory.append(deepcopy(M_pos))
                    Numb_stra[p] = Numb_stra[p]+1
                    M_pre = deepcopy(M_pos)
                    U_depend = Utilities_Poymatrix(G.m(),G.Q(),U_depend,S,S_new,Numb_stra_S)
                    U_p, S = IndUtilities(G.m(), G.c(), G.Q(), S, U_p, S_new)
                    S_new = [[] for _ in range(G.m())]
                    M_pos = [None,None,deepcopy(Numb_stra),deepcopy(S),deepcopy(U_p),deepcopy(U_depend),deepcopy(S_new),None,0]
                    print(list_best)
                    list_best.append(p)
                    list_best = list_best[:aux_p]+list_best[aux_p+1:]
                aux_p = aux_p+1
            #count = count+1
            if aux:
                # f = open("log_ComputeOneNE.txt",'a')
                # f.write("\n Total time: "+str(time()-cpu_time))
                # f.close()
                return ne, Profits, S,count,numb_back,time()-time_aux
        count = count +1
    if time()-time_aux>3600:
        print("Time Limit Exceeded")
    else:
        print(" Maximum number of iterations was attained")
    return ne_previous, [], S,count,numb_back,time()-time_aux

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
    ne: list of vectors indicating the Nash equilibrium associated with S
    Profits: List of profits for each player under ne
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
        return [],[],S,0,0,0
    Numb_stra = [len(S[p]) for p in range(G.m())]
    U_depend = [[{} for k in range(G.m())] for p in range(G.m())]
    # set mandatory strategy in the support
    Numb_stra_S = [0]*G.m()
    # STEP 2 - COMPUTE EQUILIBRIA OF RESTRICTED GAME
    # compute Nash equilibrium taking account M and Back
    count = 1
    U_depend = Utilities_Poymatrix(G.m(),G.Q(),U_depend,S_new,S,Numb_stra_S)
    list_best = list(range(G.m()))
    time_aux = time()
    ne = [0 for p in range(G.m()) for _ in range(Numb_stra[p]) ]
    deviator = G.m()-1
    ne = []
    ne_previous = ne[:]
    #return U_depend,U_p,Numb_stra,ne,deviator,S
    while True and count <= max_iter and time()-time_aux<=3600:
        print("\n\n Processing node ... ", count)
        print("Computing equilibria.... \n")
        ne_previous = ne[:]
        #signal.signal(signal.SIGALRM, signal_handler)
        #signal.alarm(3600-int(time()-time_aux))   #  seconds
        try:
            ne, Profits = ComputeNE_NOT_DFS(U_depend,U_p,G.m(),G.n_I(),G.n_C(),Numb_stra,opt_solver,ne,deviator)
            ### modify ###
            #return ne,Profits,S,count,time()-time_aux
        #except Exception, msg: python2.7
        except Exception:
            print("Time limit exceeded")
            return ne_previous, [], S,count,time()-time_aux
        print(" Equilibrium computed sucessfully")
        aux = True # no player has incentive to deviate
        S_new = [[] for p in range(G.m())] # set of new strategies to be considered
        Profile = [np.array([sum(ne[int(k+sum(Numb_stra[:p]))]*S[p][k][i] for k in range(Numb_stra[p]))  for i in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())]
        aux_p = 0
        while aux and aux_p<G.m(): # FEED BEST RESPONSES WITH NE solution
            p = list_best[aux_p]
            s_p, u_max, _ = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,Best_m[p])
            s_p = Convert_to_maximal(s_p,G,p)
            if Profits[p] +10**-6 <= u_max:
                aux = False
                S_new[p].append(s_p)
                Numb_stra_S = deepcopy(Numb_stra)
                Numb_stra[p] = Numb_stra[p]+1
                U_depend = Utilities_Poymatrix(G.m(),G.Q(),U_depend,S,S_new,Numb_stra_S)
                U_p, S = IndUtilities(G.m(), G.c(), G.Q(), S, U_p, S_new)
                S_new = [[] for _ in range(G.m())]
                list_best.append(p)
                list_best = list_best[:aux_p]+list_best[aux_p+1:]
                deviator = p
            aux_p = aux_p+1
        if aux:
            return ne, Profits, S,count,time()-time_aux
        count = count +1
    if time()-time_aux>3600:
        print("Time Limit Exceeded")
    else:
        print(" Maximum number of iterations was attained")
    return ne_previous, [], S,count,time()-time_aux

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

def ComputeNE(M,Back,U_depend,U_p,m, n_I,n_C,Numb_stra,opt_solver,Supp_Stra,start):
    ##### HEURISTIC 6 ####
    size_pre_ne = [sum(1 for j in M[-2][int(sum(M[2][:p])):int(sum(M[2][:p+1]))] if j>10**-5) for p in range(m)]
    S0 = Heuristic6(Supp_Stra,m,size_pre_ne,n_I,n_C)
    ##### HEURISTIC 6 ####
    for S0_j,s0 in enumerate(S0[start:]):
        # now we have the support sizes
        # if Back is true, how can we start s0?
        A_supp = [None for p in range(m)]
        # ATTENTION: THE SETS IN D ARE NOT SORTED
        D = []
        for p in range(m):
            # D[p] is the set of strategies in the support of player p with size s0[p]
            # from S[p] take s0[p] strategies
            if p != M[0]:
                #D.append([candidate for candidate in combinations(range(Supp_Stra[p]),s0[p])])
                # HEURISTIC II
                D.append([candidate for candidate in combinations(HeuristicII(Supp_Stra[p],M[-2][int(sum(M[2][:p])):int(sum(M[2][:p+1]))],-1,M[-3][p]),s0[p])])
            else: # M[1] must be in the support of player M[0]
                #D.append([candidate+(M[1],) for candidate in combinations(range(M[1])+range(M[1]+1,Supp_Stra[p]),s0[p]-1)])
                # HEURISTIC II: give priority to strategies choosen in the previous ne
                D.append([candidate+(M[1],) for candidate in combinations(HeuristicII(Supp_Stra[p],M[-2][int(sum(M[2][:p])):int(sum(M[2][:p+1]))],M[1],M[-3][p]),s0[p]-1)])
        ne, Profits = Recursive_Backtracking(m,A_supp,D,0,U_depend,U_p,opt_solver,Numb_stra)
        if ne != []: # Nash equilibrium found!
            return ne, Profits,start+S0_j
    return [], [],start+S0_j

def ComputeNE_NOT_DFS(U_depend,U_p,m,n_I,n_C,Numb_stra,opt_solver,ne_previous,deviator):
    ##### HEURISTIC 6 ####
    size_pre_ne = [sum(1 for j in ne_previous[int(sum(Numb_stra[:p])):int(sum(Numb_stra[:p+1]))] if j >10**-5) for p in range(deviator)]+[sum(1 for j in ne_previous[int(sum(Numb_stra[:deviator])):int(sum(Numb_stra[:deviator+1])-1)] if j >10**-5)]+[sum(1 for j in ne_previous[int(sum(Numb_stra[:p])):int(sum(Numb_stra[:p+1]))] if j >10**-5) for p in range(deviator+1,m)]
    S0 = Heuristic6(Numb_stra,m,size_pre_ne,n_I,n_C)
    ##### HEURISTIC 6 ####
    Numb_stra_previous = deepcopy(Numb_stra)
    Numb_stra_previous[deviator] = Numb_stra_previous[deviator]-1
    for S0_j,s0 in enumerate(S0):
        # now we have the support sizes
        # if Back is true, how can we start s0?
        A_supp = [None for p in range(m)]
        # ATTENTION: THE SETS IN D ARE NOT SORTED
        D = []
        for p in range(m):
            # D[p] is the set of strategies in the support of player p with size s0[p]
            # from S[p] take s0[p] strategies
            # HEURISTIC II
            if p != deviator:
                D.append([candidate for candidate in combinations(HeuristicII(Numb_stra[p],ne_previous[int(sum(Numb_stra_previous[:p])):int(sum(Numb_stra_previous[:p+1]))],-1,[]),s0[p])])
            else:
                D.append([candidate for candidate in combinations(HeuristicII(Numb_stra[p],ne_previous[int(sum(Numb_stra_previous[:p])):int(sum(Numb_stra_previous[:p+1]))],-1,[1]),s0[p])])
        ne, Profits = Recursive_Backtracking(m,A_supp,D,0,U_depend,U_p,opt_solver,Numb_stra)
        if ne != []: # Nash equilibrium found!
            return ne, Profits
    return [], []

def HeuristicII(Supp_Stra_p,M_ne,M_1,S_new_p):
    order_index_str = chain(range(M_1),range(M_1+1,Supp_Stra_p))
    M_ne_aux = M_ne+[0 for _ in S_new_p]
    if M_ne_aux !=[]:
        return sorted(order_index_str, key = lambda x:M_ne_aux[x])
    else:
        return order_index_str

def Heuristic6(Supp_Stra,m,size_ne,n_I,n_C):
    S0 = []
    # it is n_I[p]+n_C[p]+1 in case the objectives are linear: an optimum is in a facet which has dimension n_I+n_C-1
    # and thus, any point of it can be written as a convex combinatiuon of n_I+n_C extreme points of that facet
    for s0 in product(*[range(1,min(Supp_Stra[p]+1,n_I[p]+n_C[p]+2)) for p in range(m)]):
        S0.append(list(s0))
    if m == 2:
        return sorted(S0,key =lambda x:(abs(x[0]-x[1]),max(abs(size_ne[0] -x[0]),abs(size_ne[1] -x[1])),max(abs(size_ne[0]+1-x[0]),abs(size_ne[1]+1-x[1])),x[0]+x[1]))
    else:
        return sorted(S0,key =lambda x:(max(abs(size_ne[p]-x[p]) for p in range(m)),max(abs(size_ne[p]+1-x[p]) for p in range(m)),sum(x),max(abs(x[i]-x[j]) for i in range(m) for j in range(i,m))))

#######################################################################################################################

#######################################################
##        RECURSIVE BACKTRACKING                     ##
#######################################################

# INPUT
# m = number of players
# A_supp = set of strategies in the support for each player (list)
# D = candidates to be a support for each player (list)
# i = player for whom we are fixing the strategy
# S = set of strategies for each player in the restricted game (list)
# U_depend = polymatrix of utilities
# U_p = indepedent utilities
# opt_solver = 0 then use CPLEX, = 1 then use Gurobi
# Numb_stra = number strategies for each player (list)
# OUPUT
# ne - Nash equilibrium (list)
# Profits - profits for each player in the equilibrium ne (list)

def Recursive_Backtracking(m,A_supp,D,i,U_depend, U_p,opt_solver,Numb_stra):
    if i == m: # this is, we have fixed the support for each player
        # Solve Feasibility Problem
        return FeasibilityProblem(m,A_supp,U_depend,U_p,opt_solver,Numb_stra) # ne, Profits
    else:
        while D[i]!=[]:
            d_i = D[i].pop() # remove d_i from D[i]
            A_supp[i] = d_i
            D_new = RS([[A_supp[p]] for p in range(i+1)]+deepcopy(D[i+1:]), Numb_stra,U_depend, U_p,m)
            if D_new != None:
                ne, Profits = Recursive_Backtracking(m,deepcopy(A_supp),deepcopy(D_new),i+1,U_depend, U_p,opt_solver,Numb_stra)
                if ne !=[]:
                    return ne, Profits
    return [],[]

#######################################################################################################################

#######################################################
##        FEASIBILITY PROBLEM                       ##
#######################################################

# INPUT
# m = number of players
# A_supp  = strategies to which each player associates positive probability (list)
# U_depend = polymatrix of utilities
# U_p = independent utilities
# opt_solver = 0 then use CPLEX, = 1 then use Gurobi
# Numb_stra = number of strategies for each player (list)
# OUTPUT
# ne = Nash equilibrium (list)
# Profits = profit of each player for the equilibrium ne (list)

def FeasibilityProblem(m,A_supp, U_depend,U_p, opt_solver,Numb_stra):
    return FeasibilityProblem_Gurobi(m,A_supp, U_depend,U_p,Numb_stra)

def FeasibilityProblem_Gurobi(m,A_supp, U_depend,U_p,Numb_stra,m_p = None):
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
                else:
                    m_p.addConstr(U_p[p][sp]+grb.quicksum(sigma[k][sk]*U_depend[p][k][(sp,sk)] for k in range(m) if k != p for sk in A_supp[k]) <= v[p])
                    m_p.update()
        #m_p.write("apagar.lp")
        m_p.optimize()
        ne = []
        Profits = []
        #print "Solution status for Feasibility Problem: ", m_p.status
        if m_p.status not in [3,4]:
            for p, sp in enumerate(Numb_stra):
                for j in range(sp):
                    if j in A_supp[p]:
                        ne.append(sigma[p][j].x)
                    else:
                        ne.append(0)
                Profits.append(v[p].x)
        return ne, Profits

#######################################################################################################################

#######################################################
##    REMOVAL OF STRICTLY DOMINATED STRATEGIES       ##
#######################################################

# INPUT
# D = list of strategies to which each player is restricted to play (list; D[p] is a list that contains tuples = sets of strategies)
# Numb_stra = strategies available for each player (list)
# U_depend = polymatrix of utilities
# U_p = indepedent utilities
# OUTPUT
# D_new =  D_new[p] strategies that are not strictly dominated given D[-p]

# REMARK: conditionally dominated
# sp in S[p] is conditionally dominated given a profile of sets of available actions R[-p] contained in S[-p],
# if the following conditions holds:
# there is sp' in Sp[p], forall s[-p] in R[-p]: Profit[p](sp,s[-p]) < Profit[p](sp',s[-p])

def RS(D,Numb_stra,U_depend, U_p,m):
    changed = True
    while changed:
        changed = False
        for p in range(m):
            for a in set(dp for Ap in D[p] for dp in Ap): # for all pure strategies in the possible supports
                for a_prime in chain(range(a),range(a+1,Numb_stra[p])):
                    aux = True # it is conditionally dominated
                    # all possible outcomes
                    Outcomes_minus_p = [set(d_minus_p for A_minus_p in D[k] for d_minus_p in A_minus_p) for k in range(p)]+[set([0])]+[set(d_minus_p for A_minus_p in D[k] for d_minus_p in A_minus_p) for k in range(p+1,m)]
                    for s in product(*Outcomes_minus_p):
                        # if a is conditionally dominated by a_prime given D[-p]:
                        if U_p[p][a]+sum(U_depend[p][k][(a,s[k])] for k in range(m) if k !=p) >= U_p[p][a_prime]+sum(U_depend[p][k][(a_prime,s[k])] for k in range(m) if k !=p):
                            aux = False
                            break
                    if aux: # it is dominated
                        D_new_p = []
                        for Ap in D[p]:
                            if a not in Ap:
                                D_new_p.append(Ap)
                        D[p] = D_new_p
                        changed = True
                        if D[p] == []:
                            #print "### Not solving feasibility problem "
                            return None
    return D

if __name__ == "__main__":
    G = Game("KEG",2,20,26,50,3)
    S_ini = [[[]] for p in range(G.m())]
    Profile_ini = [np.array([1 for k in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())]
    for p in range(G.m()):
        S_ini[p][0], _ ,_ = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile_ini,p,False)
    ne, Profits_mSGM,S,numb_iter,numb_back,cpu_time  = IterativeSG(G,50,1,S_ini)
