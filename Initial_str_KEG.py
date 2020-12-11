# ONLY FOR KEG: it removes quadratic term in best reaction which impacts the running time of creating the optimization model

from Instances import *
# get optimization software
import gurobipy as grb
import numpy as np

#######################################################
##        COMPUTE INITIAL SET OF STRATEGIES          ##
## Generate initial strategies: monopoly strategies  ##
#######################################################

# INPUT
# G : game class (see Instances.py)

# OUTPUT
# S = list of strategies
# U_p = individial profit for each player and each strategy in S
# Best_m = list of the players best reaction models

def InitialStrategies(G,opt_solver=1):
    ### for each player produce the optimal strategy if she was alone in the game ###
    S = [[[]] for p in range(G.m())] # list of strategies
    U_p = [[[]] for p in range(G.m())] # associated individual profit
    # Profile is the profile of strategies
    Profile = [np.array([0 for k in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())]
    # Best reaction models
    Best_m = []
    for p in range(G.m()):
        try:
            S[p][0], U_p[p][0], Model_p = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False)
        except:
            print("Player ", p+1, " has no feasible solution or the problem is unbounded")
        Best_m.append(Model_p)
    return S, U_p, Best_m

#######################################################################################################################
#### ALTERNATIVE INITIALIZATION: social optimum or potential #######################################################################
#######################################################################################################################
# social = 1 then social optimum, otherwise, potential
def InitialStrategiesII(G,opt_solver=1,social=1):
    if social:
        S, U_p = SocialOptimumGurobi(G.m(), G.n_I(), G.n_C(), G.n_constr(), G.c(), G.Q(), G.A(), G.b())
    else:
        S, U_p = PotentialNEGurobi(G.m(), G.n_I(), G.n_C(), G.n_constr(), G.c(), G.Q(), G.A(), G.b())
    return S, U_p, CreateModels(G.m(), G.n_I(), G.n_C(), G.n_constr(), G.c(), G.Q(), G.A(), G.b())

#######################################################################################################################
#### ALTERNATIVE INITIALIZATION: pure NE for potential part
#### Potential  ###########################################################
#######################################################################################################################
def PotentialNEGurobi(m, n_I, n_C, n_constr, c, Q, A, b):
    m_Pot = grb.Model("PotentialNE")
    m_Pot.setParam( 'OutputFlag', False )
    m_Pot.setParam("Threads", 2)
    # set objective function direction
    m_Pot.ModelSense = -1 # maximize
    m_Pot.update()
    x = [np.array([m_Pot.addVar(vtype="B", name="x_"+str(p)+'_'+str(i)) for i in range(n_I[p])]+[m_Pot.addVar(lb=0, vtype="C", name="x_"+str(p)+'_'+str(i)) for i in range(n_I[p],n_I[p]+n_C[p])]) for p in range(m)]
    m_Pot.update()
    QuadPart = grb.QuadExpr(0)
    for p in range(m):
        for k in range(n_constr[p]):
            m_Pot.addConstr(np.dot(x[p],A[p][k]),grb.GRB.LESS_EQUAL, b[p][k])
            m_Pot.update()
        if n_I[p]+n_C[p] ==1:
            QuadPart = QuadPart+grb.QuadExpr(x[p][0]*c[p][0]-0.5*x[p][0]*Q[p][p]*x[p][0])
        else:
            QuadPart = QuadPart+grb.QuadExpr(grb.quicksum(0.5*np.dot(x[j],np.dot((Q[p][j]+Q[j][p].T),x[p].T)) for j in range(p))+np.dot(x[p],c[p])-0.5*(np.dot(x[p],np.dot(Q[p][p],x[p].T))))
    m_Pot.setObjective(QuadPart)
    m_Pot.update()
    #m_Pot.write("apagar.lp")
    m_Pot.optimize()
    try:
        S = [[[x[p][k].x for k in range(n_I[p]+n_C[p])]] for p in range(m)]
        U_p = [[float(np.dot(c[p],S[p][0])-0.5*np.dot(S[p][0],np.dot(Q[p][p],S[p][0])))] for p in range(m)]
        return S,U_p
    except:
        print("No feasible profile of strategies", m_Pot.status)
    return None

#######################################################################################################################
#### ALTERNATIVE INITIALIZATION: pure NE for potential part ###########################################################
#### Social optimum #######################
#######################################################################################################################

def SocialOptimumGurobi(m, n_I, n_C, n_constr, c, Q, A, b):
    m_SO = grb.Model("SocialOptimum")
    m_SO.setParam( 'OutputFlag', False )
    m_SO.setParam("Threads", 2)
    # set objective function direction
    m_SO.ModelSense = -1 # maximize
    m_SO.update()
    x = [np.array([m_SO.addVar(vtype="B", name="x_"+str(p)+'_'+str(i)) for i in range(n_I[p])]+[m_SO.addVar(lb=0, vtype="C", name="x_"+str(p)+'_'+str(i)) for i in range(n_I[p],n_I[p]+n_C[p])]) for p in range(m)]
    m_SO.update()
    QuadPart = grb.QuadExpr(0)
    for p in range(m):
        for k in range(n_constr[p]):
            m_SO.addConstr(np.dot(x[p],A[p][k]),grb.GRB.LESS_EQUAL, b[p][k])
            m_SO.update()
        if n_I[p]+n_C[p] ==1:
            QuadPart = QuadPart+grb.QuadExpr(x[p][0]*c[p][0]-0.5*x[p][0]*Q[p][p]*x[p][0])
        else:
            QuadPart = QuadPart+grb.QuadExpr(grb.quicksum(np.dot(x[j],np.dot(Q[p][j],x[p].T)) for j in range(m) if j !=p)+np.dot(x[p],c[p])-0.5*(np.dot(x[p],np.dot(Q[p][p],x[p].T))))
    m_SO.setObjective(QuadPart)
    m_SO.update()
    #m_SO.write("apagar.lp")
    m_SO.optimize()
    try:
        S = [[[x[p][k].x for k in range(n_I[p]+n_C[p])]] for p in range(m)]
        U_p = [[float(np.dot(c[p],S[p][0])-0.5*np.dot(S[p][0],np.dot(Q[p][p],S[p][0])))] for p in range(m)]
        return S,U_p
    except:
        print("No feasible profile of strategies", m_SO.status)
    return None

######################################################


def CreateModels(m, n_I, n_C, n_constr, c, Q, A, b):
    Profile = [np.array([0 for k in range(n_I[p]+n_C[p])]) for p in range(m)]
    Best_m = []
    for p in range(m):
        _,_,Model_p = BestReactionGurobi(m,n_I[p],n_C[p],n_constr[p],c[p],Q[p],A[p],b[p],Profile,p,True)
        Best_m.append(Model_p)
    return Best_m

##################################################################################
##############     RESTRICTED STRATEGY METHOD      ###############################
##############                                     ###############################
#############  to compute (EXACT) Nash equilibria  ###############################
##################################################################################

# Compute Best Reaction of player against the strategy 'Profile'
def BestReactionGurobi(m,n_I_p,n_C_p,n_constr_p,c_p,Q_p,A_p,b_p,Profile,p,create, m_p = None,CE_verify=False):
    if CE_verify:
        xk_Qkp = sum(Q_p[k] for k in range(m) if k!=p)
    else:
        xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p) # np.array
    if m_p == None:
        # initiate model
        m_p = grb.Model("MIQPG")
        # no pritting of the output
        m_p.setParam( 'OutputFlag', False )
        m_p.setParam("Threads", 2)
        #m_p.setParam('BarHomogeneous', 1)
        #m_p.setParam('DualReductions',0)
        # set objective function direction
        m_p.ModelSense = -1 # maximize
        m_p.update()
        # binary variables
        x = [] # decision vector
        for i in range(n_I_p):
            #x.append(m.addVar(vtype="B", obj = float(c_p[i]+xk_Qkp[i]), name="x"+str(i)))
            x.append(m_p.addVar(vtype="B", name="x"+str(i)))
            m_p.update()
        for i in range(n_I_p,n_C_p+n_I_p):
            #x.append(m.addVar(lb=0, vtype="C", obj = float(c_p[i]+xk_Qkp[i]), name="x"+str(i)))
            x.append(m_p.addVar(lb=0, vtype="C", name="x"+str(i)))
            m_p.update()
        x = np.array(x)
        # constraints
        for k in range(n_constr_p):
            #m_p.addConstr(np.array(x).np.dot(A_p[k]), grb.GRB.LESS_EQUAL, b_p[k])
            #m_p.addConstr(x.np.dot(A_p[k]), grb.GRB.LESS_EQUAL, b_p[k])
            m_p.addConstr(np.dot(A_p[k],x), grb.GRB.LESS_EQUAL, b_p[k])
            m_p.update()
        if n_I_p+n_C_p ==1:
            #QuadPart = grb.QuadExpr(x[0]*c_p[0]+ xk_Qkp*x[0]-0.5*x[0]*Q_p[p]*x[0])
            QuadPart = grb.QuadExpr(x[0]*c_p[0])
        else:
            #QuadPart = grb.QuadExpr(x.np.dot(c_p)-0.5*(x.np.dot(Q_p[p].np.dot(x.T)))+xk_Qkp.np.dot(x.T))
            #QuadPart = grb.QuadExpr(x.np.dot(c_p)-0.5*(x.np.dot(Q_p[p].np.dot(x.T))))
            QuadPart = grb.QuadExpr(np.dot(c_p,x))
        m_p.setObjective(QuadPart)
        m_p.update()
    if n_I_p+n_C_p ==1 and type(xk_Qkp) is not np.ndarray:
        if CE_verify and m_p!=None:
            # when we use CE, we change objective function in the indepedent part
            #QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
            x_tmp = m_p.getVars()
            QuadPart = grb.QuadExpr(x_tmp[0]*c_p[0])
            m_p.setObjective(QuadPart)
            m_p.update()
            m_p.setObjective(m_p.getObjective()+xk_Qkp*x_tmp)
            m_p.update()
        else:
            m_p.setObjective(m_p.getObjective()+xk_Qkp*m_p.getVars()[0])
            m_p.update()
    else:
        if CE_verify and m_p!=None:
            x_tmp = np.array(m_p.getVars())
            #QuadPart = grb.QuadExpr(np.dot(c_p,m_p.getVars())-0.5*(np.dot(np.dot(m_p.getVars().T,Q_p[p]),m_p.getVars())))
            QuadPart = grb.QuadExpr(np.dot(c_p,x_tmp))
            #QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
            m_p.setObjective(QuadPart)
            m_p.update()
            m_p.setObjective(m_p.getObjective()+np.dot(x_tmp,xk_Qkp))
            m_p.update()
        else:
            m_p.setObjective(m_p.getObjective()+np.dot(m_p.getVars(),xk_Qkp))
            m_p.update()
    #m_p.write("apagar.lp")
    # create is always false
    if create:
        if n_I_p+n_C_p ==1 and type(xk_Qkp) is not np.ndarray:
            if CE_verify:
                # when we use CE, we change objective function in the indepedent part
                #QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
                x_tmp = m_p.getVars()
                QuadPart = grb.QuadExpr(x_tmp[0]*c_p[0])
                # overrite objective
                m_p.setObjective(QuadPart)
                m_p.update()
            m_p.setObjective(m_p.getObjective()-xk_Qkp*m_p.getVars()[0])
        else:
            if CE_verify:
                #QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
                x_tmp = np.array(m_p.getVars())
                QuadPart = grb.QuadExpr(np.dot(c_p,x_tmp))
                # overrite objective
                m_p.setObjective(QuadPart)
                m_p.update()
            m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()))
        m_p.update()
        return None,None,m_p
    if not CE_verify:
        # warm start
        for j,aux_var in enumerate(m_p.getVars()):
            aux_var.start = Profile[p][j]
            m_p.update()
    m_p.optimize()
    try:
        #return [x[i].x for i in range(n_I_p+n_C_p)],m_p.ObjVal, m_p
        sol = [i.x for i in m_p.getVars()]
        value = m_p.ObjVal
        if n_I_p+n_C_p ==1 and type(xk_Qkp) is not np.ndarray:
            # this is important for NE verification
            if CE_verify:
                m_p.setObjective(m_p.getObjective()-xk_Qkp*x_tmp)
            else:
                m_p.setObjective(m_p.getObjective()-xk_Qkp*m_p.getVars()[0])
        else:
            if CE_verify:
                m_p.setObjective(m_p.getObjective()-np.dot(x_tmp,xk_Qkp))
            else:
                # this is important for NE verfication
                m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()))
        m_p.update()
        return sol, value, m_p
    except:
        print("Wow! The best reaction problem has no feasible solution. The status code is: ", m_p.status)

if __name__ == "__main__":
    G = Game("KEG",2,80,2,50,3)

    # Verify best response
    S = [[[]] for p in range(G.m())] # list of strategies
    U_p = [[[]] for p in range(G.m())] # associated individual profit
    # Profile is the profile of strategies
    Profile = [np.array([0 for k in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())]
    # Best reaction models
    Best_m = []
    p=0
    S[p][0], U_p[p][0], Model_p = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False)
    #S, U_p, Best_m = InitialStrategies(G,1)

    # verify initial STRATEGIES
    #S, U_p, Best_m = InitialStrategies(G)

    # verify initial strategies II: uses potential part of the game
    #S_II, U_p_II, Best_m_II = InitialStrategiesII(G,1,1)
