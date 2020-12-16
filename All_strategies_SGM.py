__author__ = 'Carvalho'

from Compute_NE_all_strategies import *

# compute all feasible solutions
def All_Strategies(A,b,n,m):
    S = [[] for _ in range(m)]
    for numb in range(n+1): # number of items selected
        for set_items in combinations(range(n), numb): # select numb items
            aux = [0]*n
            for i in set_items:
                aux[i] = 1
            for p in range(m):
                if sum(A[p][0][i] for i in set_items) <= b[p]:
                    S[p].append(aux)
    return S


if __name__ == "__main__":
    #from RestrictedStrategyMethod import Knapsack_RandomGame
    #for m in [2,3]:
    for m in [2]:
        #for n in [5,7,10]:
        for n in [5]:
            f_knapsack = open("KnapsackGeneralTable_only_PNS.txt",'a')
            f_knapsack.write(str(n)+" & "+str(m)+" & ")
            f_knapsack.close()

            # DFS: use the information about the player with incentive to deviate
            total_time = 0
            total_iter = 0
            total_size_0 = 0
            total_size_1 = 0
            total_size_2 = 0
            numb_pNE = 0
            numb_mNE = 0
            aux = 0

            # all strategies
            time_all = 0
            total_size_0_PNS = 0
            total_size_1_PNS = 0
            total_size_2_PNS = 0
            numb_pNE_PNS = 0
            numb_mNE_PNS = 0
            aux_PNS = 0

            #for ins in range(10):
            for ins in [0]:
                filename  ="Instances/KnapsackSmall/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".npy"
                G = Game("empty")
                G.Read_Game(filename)
                #G = Game("Knapsack",m,n,ins)
                #G.Save_Game(m,n,ins)

                #######################################################################################################
                # DFS: EXECUTE
                max_iter=50
                #ne, Profits,S,numb_iter,numb_back,cpu_time  = IterativeSG(G,max_iter)
                ne, Profits,S,numb_iter,cpu_time  = IterativeSG_NOT_DFS(G,max_iter)

                f_knapsack = open("KnapsackGeneralTable_only_PNS.txt",'a')
                #f_knapsack.write('ins & time & iter & pNE & mNE & size S & numb_back \\\ \n')
                if cpu_time>3600:
                    f_knapsack.write(" & & "+str(ins)+' & '+'tl & '+str(numb_iter)+' & 0 & 0 & '+str([len(S[i]) for i in range(m)])+' & ')
                else:
                    if numb_iter >= max_iter:
                        f_knapsack.write(" & & "+str(ins)+' & '+"%0.2f" % cpu_time+' & '+str(numb_iter)+' & 0 & 0 & '+str([len(S[i]) for i in range(m)])+' & ')
                    else:
                        total_time = total_time + cpu_time
                        total_iter = total_iter + numb_iter
                        total_size_0 = total_size_0 + len(S[0])
                        total_size_1 = total_size_1 + len(S[1])
                        if m>2:
                            total_size_2 = total_size_2+len(S[2])
                        aux = aux +1
                        if sum(1 for s in ne if s>=1-10**-6) <m:
                            f_knapsack.write(" & & "+str(ins)+' & '+"%0.2f" % cpu_time+' & '+str(numb_iter)+' & 0 & 1 & '+str([len(S[i]) for i in range(m)])+' & ')
                            numb_mNE = numb_mNE+1
                        else:
                            f_knapsack.write(" & & "+str(ins)+' & '+"%0.2f" % cpu_time+' & '+str(numb_iter)+' & 1 & 0 & '+str([len(S[i]) for i in range(m)])+' & ')
                            numb_pNE = numb_pNE+1
                f_knapsack.close()
                # DFS: DONE
                #######################################################################################################
                # Generate all possible strategies
                S_all = All_Strategies(G.A(),G.b(),n,m)
                ne_PNS,Profits_PNS,S_PNS,numb_iter_not,cputime_not_dfs = IterativeSG_NOT_DFS(G,max_iter,1,S_all)
                f_knapsack = open("KnapsackGeneralTable_only_PNS.txt",'a')
                if cputime_not_dfs>3600 or Profits_PNS==[]:
                    f_knapsack.write(' & tl & 0  &  0 '+str([len(S_PNS[i]) for i in range(m)])+' & \\\ \n')
                else:
                    total_size_0_PNS = total_size_0_PNS + len(S_PNS[0])
                    total_size_1_PNS = total_size_1_PNS + len(S_PNS[1])
                    if m>2:
                        total_size_2_PNS = total_size_2_PNS+len(S_PNS[2])
                    aux_PNS = aux_PNS+1
                    time_all = time_all+cputime_not_dfs
                    if sum(1 for s in ne if s>=1-10**-6) <m:
                        numb_mNE_PNS = numb_mNE_PNS+1
                        f_knapsack.write(' & '+"%0.2f" % cputime_not_dfs+' & 0 & 1 & '+str([len(S_PNS[i]) for i in range(m)])+' & & \\\ \n')
                    else:
                        numb_pNE_PNS = numb_pNE_PNS+1
                        f_knapsack.write(' & '+"%0.2f" % cputime_not_dfs+' & 1 & 0 & '+str([len(S_PNS[i]) for i in range(m)])+' & & \\\ \n')
                f_knapsack.close()



                #######################################################################################################
            # statistics to DFS
            total_time = (total_time*1.)/aux
            total_iter = (total_iter*1.)/aux
            total_size_0 = (total_size_0*1.)/aux
            total_size_1 = (total_size_1*1.)/aux
            total_size_2 = (total_size_2*1.)/aux
            f_knapsack = open("KnapsackGeneralTable_only_PNS.txt",'a')
            if m >2:
                f_knapsack.write(" & & & "+ "%0.2f" % total_time+' & '+"%0.2f" % total_iter+' & '+"%0.2f" % total_size_0+' & '+"%0.2f" % total_size_1+' & '+"%0.2f" % total_size_2+' & '+str(numb_pNE)+' & '+str(numb_mNE))
            else:
                f_knapsack.write(" & & & "+ "%0.2f" % total_time+' & '+"%0.2f" % total_iter+' & '+"%0.2f" % total_size_0+' & '+"%0.2f" % total_size_1+' &  & '+str(numb_pNE)+' & '+str(numb_mNE))
            f_knapsack.close()

            time_all = (time_all*1.)/aux_PNS
            total_size_0_PNS = (total_size_0_PNS*1.)/aux_PNS
            total_size_1_PNS = (total_size_1_PNS*1.)/aux_PNS
            total_size_2_PNS = (total_size_2_PNS*1.)/aux_PNS
            f_knapsack = open("KnapsackGeneralTable_only_PNS.txt",'a')
            if m>2:
                f_knapsack.write(' & '+ "%0.2f" % time_all+' & '+"%0.2f" % total_size_0_PNS+' & '+"%0.2f" % total_size_1_PNS+' & '+"%0.2f" % total_size_2_PNS+" & "+str(numb_pNE_PNS)+" & "+str(numb_mNE_PNS)+' \\\ \n')
            else:
                f_knapsack.write(' & '+ "%0.2f" % time_all+' & '+"%0.2f" % total_size_0_PNS+' & '+"%0.2f" % total_size_1_PNS+' &  & '+str(numb_pNE_PNS)+" & "+str(numb_mNE_PNS)+' \\\ \n')
            f_knapsack.close()
