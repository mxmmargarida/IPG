#from numpy import *
import numpy as np

class Game(object):
    r"""Create instances in a standard format.

    Parameters:
    ----------
    type: string with Knapsack or LotSizing
    m: number of players (optional),
    n: number of items for Knapsack, number of period for LotSizing (optional), number of vertices KEG
    ins: number associated with the instance (Knapsack only), must be <numb_ins
    numb_ins: number of instances (Knapsack only)
    K: maximum cycle length allowed (KEG only)
    Returns:
    -------
    m: number of players
    n_I: list of number of binary variables for each player i=0,..,m-1
    n_C: list of number of continuous variables for each player i=0,..,m-1
    n_constr: list of number of constraints for each player i=0,..,m-1
    c: list of coeficients for the linear part of the obj function,
        for each player i=0,..,m-1
    Q: list of matrices for the bilinear part of the obj function,
        for each player i=0,..,m-1
    A: list of constraint matrices for each player i=0,..,m-1
    b: list of vectors with the rhs  for the constraints of each player i=0,..,m-1
    """
    def __init__(self, type, m = 2, n = 10, ins = 0, numb_ins = 10, K=0):
        if type == 'Knapsack':
            m, n_I, n_C, n_constr, c, Q, A, b = Knapsack_Game(m,n,ins,numb_ins)
        elif type == 'LotSizing':
            # n serves as number of items or number of periods: n = T
            T = n
            m, n_I, n_C, n_constr, c, Q, A, b, A_market, B, F, H, C, M = LotSizing_Game(m,T)
        elif type == 'KEG':
            m, n_I, n_C, n_constr, c, Q, A, b = Two_KEG_RandomGame(n,ins,K)
        elif type == "empty":
            m = 0
            n_I = []
            n_C = []
            n_constr = []
            c, Q, A, b = [], [], [], []
        else:
            print("Not valid instance")
            raise NameError('Give a proper type to the game')
        self.__m = m
        self.__n_I = n_I
        self.__n_C = n_C
        self.__n_constr = n_constr
        self.__c = c
        self.__Q = Q
        self.__A = A
        self.__b = b
        self.__type = type
        self.__ins = ins
        self.__numb_ins = numb_ins

    # give parameters of a player
    def Player_n_I(self,p):
        if p > self.__m:
            raise NameError('That player does not exist')
        else:
            return self.__n_I[p]
    def Player_n_C(self,p):
        if p > self.__m:
            raise NameError('That player does not exist')
        else:
            return self.__n_C[p]
    def Player_n_constr(self,p):
        if p > self.__m:
            raise NameError('That player does not exist')
        else:
            return self.__n_constr[p]
    def Player_c(self,p):
        if p > self.__m:
            raise NameError('That player does not exist')
        else:
            return self.__c[p]
    def Player_Q(self,p):
        if p > self.__m:
            raise NameError('That player does not exist')
        else:
            return self.__Q[p]
    def Player_A(self,p):
        if p > self.__m:
            raise NameError('That player does not exist')
        else:
            return self.__A[p]
    def Player_b(self,p):
        if p > self.__m:
            raise NameError('That player does not exist')
        else:
            return self.__b[p]

    def Numb_players(self):
        return self.__m

    def type(self):
        return self.__type
    def ins(self):
        return self.__ins
    def numb_ins(self):
        return self.__numb_ins
    def b(self):
        return self.__b
    def A(self):
        return self.__A
    def Q(self):
        return self.__Q
    def c(self):
        return self.__c
    def n_constr(self):
        return self.__n_constr
    def n_C(self):
        return self.__n_C
    def n_I(self):
        return self.__n_I
    def m(self):
        return self.__m

    def Save_Game(self,m=2,n=10,ins=0):
        # save file with instance
        filename ='Instances/'+self.__type+"/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".npy"
        with open(filename,"wb") as f:
            np.save(f,self.__m)
            np.save(f,self.__n_I)
            np.save(f,self.__n_C)
            np.save(f,self.__n_constr)
            np.save(f,self.__c)
            np.save(f,self.__Q)
            if self.__type == "KEG": # then number of restrictions can vary an numpy raises an error if we save A all together
                aux = [self.__A[p] for p in range(self.__m)]
                finalc = np.empty(len(aux),dtype=object)
                finalc[:]=aux
                np.save(f,finalc)
            else:
                np.save(f,self.__A)
            np.save(f,self.__b)
            np.save(f,self.__type)
            np.save(f,self.__ins)
            np.save(f,self.__numb_ins)

    def Read_Game(self, filename):
        if self.__type=="empty":
            with open(filename,"rb") as f:
                self.__m = int(np.load(f))
                self.__n_I = list(np.load(f))
                self.__n_C = list(np.load(f))
                self.__n_constr = list(np.load(f))
                self.__c = list(np.load(f,allow_pickle=True))
                self.__Q = list(np.load(f,allow_pickle=True))
                self.__A = list(np.load(f,allow_pickle=True))
                self.__b = list(np.load(f))
                self.__type = str(np.load(f))
                self.__ins = int(np.load(f))
                self.__numb_ins = int(np.load(f))
        else:
            raise NameError("It is not an empty game")

    # create game manually
    def Create(self,m,n_I,n_C,n_constr,c,Q,A,b,type="empty",ins=1,numb_ins=1):
        self.__m = m
        self.__n_I = n_I
        self.__n_C = n_C
        self.__n_constr = n_constr
        self.__c = c
        self.__Q = Q
        self.__A = A
        self.__b = b
        self.__type = type
        self.__ins = ins
        self.__numb_ins = numb_ins


    # recover info on lot sizing game
    def A_market(self):
        if self.__type=="LotSizing": # recover A_market
            T = self.__n_I[0] # number of periods equal to number of binary variables
            # market size: part of linear objective
            return list(self.__c[0][2*T:3*T])
    def B(self):
        if self.__type=="LotSizing": # recover B
            T = self.__n_I[0] # number of periods equal to number of binary variables
            # market slope: part of quadratic objective
            return [int(-0.5*self.__Q[0][0][2*T+i,2*T+i]) for i in range(T)]
    def F(self):
        if self.__type=="LotSizing": # recover F
            T = self.__n_I[0] # number of periods equal to number of binary variables
            # F : setup costs
            return [list(self.__c[p][:T]) for p in range(self.__m)]
    def H(self):
        if self.__type=="LotSizing": # recover H
            T = self.__n_I[0] # number of periods equal to number of binary variables
            # H : inventory costs
            return [list(self.__c[p][3*T:4*T]) for p in range(self.__m)]
    def C(self):
        if self.__type=="LotSizing": # recover C
            T = self.__n_I[0] # number of periods equal to number of binary variables
            # C : production costs
            return [list(self.__c[p][T:2*T]) for p in range(self.__m)]
    def M(self):
        if self.__type=="LotSizing": # recover M
            T = self.__n_I[0] # number of periods equal to number of binary variables
            # M : production capacity per period
            return [[self.__A[p][2*T+i,i] for i in range(T)] for p in range(self.__m)]

    def __str__(self):
        return self.__type+" game"


################################################
############ KNAPSACK RANDOM GAME ##############
################################################
def Knapsack_Game(m,n, ins, numb_ins):
    # number of integer decision variables
    n_I = [n for _ in range(m)]
    # number of continuous decision variables
    n_C = [0 for _ in range(m)]
    # number of constraints
    n_constr = [1 for _ in range(m)]
    # linear gains c^p
    c = [np.random.randint(-50,50,size=n) for _ in range(m)]
    Q = [[np.diag(np.random.randint(-100,100,size=n)) for _ in range(m)] for _ in range(m)]
    for p in range(m):
        Q[p][p] = np.zeros((n,n))
    # Knapsack Capacity constraint
    A = [np.random.randint(-100,100,size=(1,n)) for _ in range(m)]
    # Knapsack Capacity
    b = [np.array([int(((ins*1.)/(numb_ins+1))*np.sum(A[p]))]) for p in range(m)]
    return m, n_I, n_C, n_constr, c, Q, A, b

#######################################################
##     GENERATE RANDOM LOT SIZING GAME               ##
#######################################################

def LotSizing_Game(m,T):
    ##### Generate problem parameters ########
    # Market Price
    B = [-1*np.random.randint(1,3) for t in range(T)]
    A_market = [np.random.randint(20,30) for t in range(T)]
    # Setup Costs
    F = [[-1*np.random.randint(10,20) for t in range(T)] for p in range(m)]
    # Variable Costs
    C = [[-1*np.random.randint(5,10) for t in range(T)] for p in range(m)]
    # Inventory Holding Costs
    #H = [[-1*randint(5,10) for t in range(T)] for p in range(m)]
    H = [[0 for t in range(T)] for p in range(m)]

    # Production Capacity
    M = [[-1*sum((A_market[t]*1.)/(-1*B[t]) for t in range(T)) for j in range(T)] for p in range(m)]
    # number of integer decision variables
    n_I = [T for p in range(m)]
    # number of continuous decision variables
    n_C = [3*T+1 for p in range(m)]
    # number of constraints
    n_constr = [3*T+2 for p in range(m)]
    # linear gains c^p
    c = [np.array(F[p]+C[p]+A_market+[0]+H[p]) for p in range(m)]
    Q = [[np.diag([0]*T+[0]*T+B+[0]*(T+1)) for k in range(m)] for p in range(m)]
    B_new = [-2*i for i in B]
    for p in range(m):
        Q[p][p] = np.diag([0]*T+[0]*T+B_new+[0]*(T+1))
    A=[]
    b = []
    for p in range(m):
        A_p1 = np.concatenate((np.zeros((T,T)),np.diag([1]*T),np.diag([-1]*T),np.concatenate((np.diag([1]*T),np.zeros((T,1))),axis=1)+np.concatenate((np.zeros((T,1)),np.diag([-1]*T)),axis=1)),axis=1)
        A_p2 = -1 * A_p1
        A_p3 = np.concatenate((np.diag(M[p]),np.diag([1]*T),np.zeros((T,T)),np.zeros((T,T+1))),axis=1)
        aux = np.zeros((2,T+1))
        aux[0,0] = 1
        aux[1,T] = 1
        A_p4 = np.concatenate((np.zeros((2,T)),np.zeros((2,T)),np.zeros((2,T)),aux), axis = 1)
        A_p = np.concatenate((A_p1,A_p2,A_p3,A_p4), axis = 0)
        A.append(A_p)
        b.append(np.zeros((3*T+2,1)))
    return m, n_I, n_C, n_constr, c, Q, A, b,A_market, B, F, H, C, M

#######################################################
##     GENERATE Kidney Exchange GAME from data set  ##
#######################################################
def Two_KEG_RandomGame(size, ins,K):
    r"""Create instances in a standard format.

    Parameters:
    ----------
    size: total number of vertices (should be 10,20,30,40,50,60,70,80,90) MUST BE EVEN
    ins: instance to read (should be between 1 and 51),
    K: maximum length for cycles
    Returns:
    -------
    m: number of players
    n_I: list of number of binary variables for each player i=0,..,m-1
    n_C: list of number of continuous variables for each player i=0,..,m-1
    n_constr: list of number of constraints for each player i=0,..,m-1
    c: list of coeficients for the linear part of the obj function,
        for each player i=0,..,m-1
    Q: list of matrices for the bilinear part of the obj function,
        for each player i=0,..,m-1
    A: list of constraint matrices for each player i=0,..,m-1
    b: list of vectors with the rhs  for the constraints of each player i=0,..,m-1
    """
    if size <=60:
        aux =str(size)+"_"+str(ins).zfill(2)+ ".input/"+str(size)+"-instance-"+ str(ins)+".input"
    else:
        aux =str(size)+"_"+str(ins).zfill(2)+ ".input/"+str(size)+"_" +str(ins).zfill(2)+".input"
    filename = "Instances/KEPDataSet/"+aux
    G, num_V = read_kep(filename)
    cycles_K = get_all_cycles(G,K)
    #the nodes from 0 to num_v/2-1 belong to player A
    cycles_A, cycles_B, cycles_IA = IdentifyCyles(cycles_K,num_V)
    if (cycles_A == [] and cycles_IA == []) or (cycles_B == [] and cycles_IA == []):
        return 500, None, None, None, None, None,None, None
    # number of integer decision variables
    n_I = [len(cycles_A)+len(cycles_IA), len(cycles_B)+len(cycles_IA)]
    # number of continuous decision variables
    n_C = [0,0]
    # number of constraints
    n_constr = [int(num_V/2),int(num_V/2)]
    b = [np.ones(int(num_V/2)), np.ones(int(num_V/2))]
    A = [np.zeros([int(num_V/2), n_I[0]]),np.zeros([int(num_V/2), n_I[1]])]
    Q = [[np.zeros([n_I[0],n_I[0]]),np.zeros([n_I[1],n_I[0]])],[np.zeros([n_I[0],n_I[1]]),np.zeros([n_I[1],n_I[1]])]]
    c = [np.zeros(n_I[0]),np.zeros(n_I[1])]
    for i,cycle_c in enumerate(cycles_A):
        c[0][i] = len(cycle_c)
        for v in cycle_c:
            A[0][v,i] = 1
    for i,cycle_c in enumerate(cycles_B):
        c[1][i] = len(cycle_c)
        for v in cycle_c:
            A[1][v-(int(num_V/2)),i] = 1
    for i,cycle_c in enumerate(cycles_IA):
        wA, wB = 0,0
        for v in cycle_c:
            if v <= num_V/2-1:
                wA = wA+1
                A[0][v,len(cycles_A)+i] = 1
            else:
                wB = wB+1
                A[1][v-(int(num_V/2)),len(cycles_B)+i] = 1
        Q[0][1][len(cycles_B)+i,len(cycles_A)+i] = wA
        Q[1][0][len(cycles_A)+i, len(cycles_B)+i] = wB
    return 2, n_I, n_C, n_constr, c, Q, A, b

# identify each player set of cycles
def IdentifyCyles(cycles_K,num_V):
    cycles_A, cycles_B, cycles_IA = [], [], []
    for c in cycles_K:
        if max(c)<=num_V/2-1:
            cycles_A.append(c)
        elif min(c)>num_V/2-1:
            cycles_B.append(c)
        else:
            cycles_IA.append(c)
    return cycles_A, cycles_B, cycles_IA

### JPP code
# READ INSTANCE
# INPUT
# filename - it is a string
# OUTPUT
# G - it is a incident list; a dictionary
# num_V - number of nodes
def read_kep(filename):
    #read file in the 'standard' kep format
    f = open(filename)
    num_V, num_E = map(int,f.readline().split())
    G = {i:[] for i in range(num_V)}
    for _ in range(num_E):
        v1, v2, w = map(int,f.readline().split()) # ignore arcs' weights
        G[v1].append(v2)
    return G, num_V

# make cycles to start in the node with smallest label
def normalize(cycle):
    cmin = min(cycle)
    while cycle[0] != cmin:
        v = cycle.pop(0)
        cycle.append(v)

# spc = ""
def all_cycles(cycles, path, node, tovisit, adj,K):
    global spc
    # spc += "    "
    # print K, spc, path, "-->", node, "tovisit:", tovisit, "graph:", adj
    for i in adj[node]:
        if i in path:
            j = path.index(i)
            cycle = path[j:]+[node]
            normalize(cycle)
            # print spc, "added cycle", cycle
            cycles.add(tuple(cycle))

        if i in tovisit:
            # print spc, "going through", node, "-", i
            if K-1 > 0:
                all_cycles(cycles,path+[node],i,tovisit-set([i]),adj,K-1)
    # spc = spc[:-4]
    return cycles


def get_all_cycles(adj,K):
    tovisit = set(adj.keys())
    visited = set([])
    cycles = set([])
    for i in tovisit:
        tmpvisit = set(tovisit)
        tmpvisit.remove(i)
        first = i
        all_cycles(cycles,[],first,tmpvisit,adj,K)
    return cycles



if __name__ == "__main__":
    np.random.seed(1)
    m = 2
    n = 5
    ins = 2
    G_KP = Game('Knapsack',m,n,ins)
    G_LS = Game('LotSizing',m,n)
    #G_KP.Save_Game(m,n,ins)
    G = Game("empty")
    #G_KP.Save_Game(m,n,ins)
    #filename = "Instances/"+G_KP.type()+"/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".npy"
    #G.Read_Game(filename)

    m = 2
    n = 30
    ins = 49
    K = 3
    G_KEG = Game('KEG',m,n,ins,50,K)
    G_KEG.Save_Game(m,n,ins)
    G = Game("empty")
    filename = "Instances/"+G_KEG.type()+"/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".npy"
    G.Read_Game(filename)
