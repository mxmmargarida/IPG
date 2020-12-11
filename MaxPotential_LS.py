import gurobipy as grb

#######################################################
##   Potential Function for Lot Sizing game          ##
#######################################################
# Potential Function for the Lot Sizing Game
# INPUT
# a - demand parameter
# b - demand parameter
# F - fixed cots
# H - inventory costs
# C - production costs
# n - number of players
# T - number of periods
# OUTPUT
# cpu_time - computational time to find a Nash equilibrium
# X - Nash equilibrium (production plan; it is always a pure equilibrium)
# Profits - utilities associated with the equilibirum X
def MaxPotentialFunction(a,b,F,H,C,M,n,T):
    # initiate model
    m = grb.Model("Potential Function")
    # no output of progress
    m.setParam( 'OutputFlag', False )
    m.setParam("TimeLimit", 3600)
    m.setParam("Threads", 2)
    # create decision variables
    x = {} # quantity to be produced by player p at time t
    y = {} # set up binary variables for player p at time t
    q = {} # quantity placed in the market by player p at time t
    h = {} # stock variables for player p at time t
    for p in range(1,n+1): h[0,p], h[T,p] = 0,0
    QuadPart = grb.QuadExpr(0)
    for t in range(1,T+1):
        Qt = grb.LinExpr(0)
        for p in range(1,n+1):
            x[t,p] = m.addVar(lb=0,vtype="C", name="x"+str(t)+"_"+str(p))
            y[t,p] = m.addVar(vtype="B",name="y"+str(t)+"_"+str(p))
            q[t,p] = m.addVar(lb=0,vtype="C", name="q"+str(t)+"_"+str(p))
            if t !=T:
                h[t,p] = m.addVar(lb=0,vtype="C", name="h"+str(t)+"_"+str(p))
            m.update()
            m.addConstr(x[t,p] <= -1*M[p-1][t-1]*y[t,p])
            m.update()
            m.addConstr(h[t-1,p]+x[t,p]-q[t,p]-h[t,p]==0)
            m.update()
            QuadPart.add(0.5*b[t-1]*q[t,p]*q[t,p]+H[p-1][t-1]*h[t,p]+C[p-1][t-1]*x[t,p]+F[p-1][t-1]*y[t,p]+a[t-1]*q[t,p])
            Qt.add(q[t,p])
        QuadPart.add(0.5*b[t-1]*Qt*Qt)
    m.setObjective(QuadPart)
    m.update()
    m.ModelSense = -1 # maximize
    m.update()
    m.optimize()
    if m.Status == 9: # time limit exceeded
        return "tl",[],[]
    X = [[] for p in range(n)]
    for t in range(1,T+1):
        for p in range(1,n+1):
            if t == T:
                X[p-1]= X[p-1]+ [y[t,p].x, x[t,p].x, q[t,p].x, 0]
            else:
                X[p-1]= X[p-1]+[y[t,p].x, x[t,p].x, q[t,p].x, h[t,p].x]
    Profits = [sum(F[p][t]*X[p][4*t]+C[p][t]*X[p][4*t+1]+H[p][t]*X[p][4*t+3]+(a[t]+b[t]*sum(X[k][4*t+2] for k in range(n)))*X[p][4*t+2] for t in range(T)) for p in range(n)]
    print('Optimal Potential value = ',m.ObjVal)
    return m.Runtime, X, Profits

if __name__ == "__main__":
    from Compute_NE import *
    m = 2 # 2 players
    n = 10 # 10 periods
    filename  ="Instances/LotSizing/Game_"+str(2)+"_"+str(10)+"_"+str(0)+".npy"
    G = Game("empty")
    G.Read_Game(filename)
    cpu_time_pot,X_pot,Profits_pot= MaxPotentialFunction(G.A_market(),G.B(),G.F(),G.H(),G.C(),G.M(),m,n)
    ne, Profits_SGM,S,numb_iter,cpu_time_not_dfs  = IterativeSG_NOT_DFS(G,50)
