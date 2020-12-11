from Compute_NE import *
from KEG_SocialOPT import *

for m in [2]:
    #for n in [20,40,80,100]:
    for n in [20]:
        f_KEG = open("KEG_table.txt",'a')
        f_KEG.write('|V| ='+str(n)+' \n\n')
        f_KEG.write(' & ins & time & iter & pNE & mNE & size S & numb_back & time & iter & pNE & mNE & size S \\\ \n')
        f_KEG.close()

        # m-SGM
        total_time = 0
        total_iter = 0
        total_size_0 = 0
        total_size_1 = 0
        numb_pNE = 0
        numb_mNE = 0
        avg_numb_back = 0
        aux = 0 # counts instances solved

         # SGM
        total_time_not_dfs = 0
        total_iter_not_dfs = 0
        total_size_0_not_dfs = 0
        total_size_1_not_dfs = 0
        numb_pNE_not_dfs = 0
        numb_mNE_not_dfs = 0
        aux_not_dfs = 0


        # utilities in the ne
        U_player_A_NE = 0
        U_player_B_NE = 0
        U_player_A_NE_not_dfs = 0
        U_player_B_NE_not_dfs = 0
        Social_welfare = 0
        Social_welfare_not_dfs = 0
        time_social = 0

        aux_no_game = 0

        #for ins in range(1,51):
        for ins in [33]:
            # read instance
            if ins!=19 or n!=20: # there is a player without strategies
                filename  ="Instances/KEG/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".npy"
                G = Game("empty")
                G.Read_Game(filename)

                # Profits alone
                Profile_alone = [np.array([0 for k in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())]
                value_p = [0,0]
                for p in range(G.m()):
                    _, value_p[p] ,_ = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile_alone,p,False)
                # save "alone" result
                filename ="Instances/KEG/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".txt"
                f_result = open(filename, "a")
                f_result.write("Profit alone = "+str(value_p)+"\n")
                f_result.close()

                # Start with the following set of strategies: anything can be selected
                S_ini = [[[]] for p in range(G.m())]
                Profile_ini = [np.array([1 for k in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())]
                for p in range(G.m()):
                    S_ini[p][0], _ ,_ = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile_ini,p,False)

                # RUN m-SGM
                # restrict to 3600 seconds use command line
                max_iter=50
                ne, Profits_mSGM,S,numb_iter,numb_back,cpu_time  = IterativeSG(G,max_iter,1,S_ini)
                # support size for 2 player game
                size_supp_ne = [sum(1 for s in ne[:len(S[0])] if s>10**-4),sum(1 for s in ne[len(S[0]):len(S[0])+len(S[1])] if s>10**-4)]

                # SAVE TABLE OF RESULTS
                f_KEG = open("KEG_table.txt",'a')
                if cpu_time>3600 or Profits_mSGM==[]:
                    f_KEG.write("&"+str(ins)+' & '+'tl & '+str(numb_iter)+' & 0 & '+str(size_supp_ne) +" & "+str([len(S[i]) for i in range(m)])+' & '+str(numb_back) +' & &')
                elif numb_iter>=max_iter:
                     f_KEG.write("&"+str(ins)+' & '+"%0.2f" % cpu_time+' & '+str(numb_iter)+' & 0 & '+str(size_supp_ne)+" & "+str([len(S[i]) for i in range(m)])+' & '+str(numb_back) +' & &')
                else:
                    total_time = total_time + cpu_time
                    total_iter = total_iter + numb_iter
                    total_size_0 = total_size_0 + len(S[0])
                    total_size_1 = total_size_1 + len(S[1])
                    avg_numb_back = avg_numb_back + numb_back
                    if Profits_mSGM[0]>10**-6: # otherwise the decrease is zero
                        U_player_A_NE = U_player_A_NE + (1 - (value_p[0]/Profits_mSGM[0])) # decrease by acting alone
                    if Profits_mSGM[1]>10**-6:
                        U_player_B_NE = U_player_B_NE + (1- (value_p[1]/Profits_mSGM[1]))
                    # SOCIAL OPTIMUM
                    Social_OPT,_, cput_time_social =  SocialOptimum(n, ins,3)
                    Social_welfare = Social_welfare + (sum(Profits_mSGM)/Social_OPT) # percentage of social welfare decrease in the game
                    time_social = time_social + cput_time_social

                    aux = aux +1
                    if sum(1 for s in ne if s>=1-10**-6) <m: #mixed NE
                        f_KEG.write("&"+str(ins)+' & '+"%0.2f" % cpu_time+' & '+str(numb_iter)+' & 0 & '+str(size_supp_ne) +" & "+str([len(S[i]) for i in range(m)])+' & '+str(numb_back) +' & &')
                        numb_mNE = numb_mNE+1
                    else:
                        f_KEG.write("&"+str(ins)+' & '+"%0.2f" % cpu_time+' & '+str(numb_iter)+' & 1 & 0 &'+str([len(S[i]) for i in range(m)])+' & '+str(numb_back) +' & & ')
                        numb_pNE = numb_pNE+1
                f_KEG.close()

                # SAVE INSTANCE SOLUTIONS
                filename ="Instances/KEG/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".txt"
                f_result = open(filename, "a")
                f_result.write("Social_OPT="+str(Social_OPT)+" \n")
                f_result.write("cpu_time_social="+str(cput_time_social)+"\n")
                f_result.write("m-SGM results \n")
                f_result.write("ne="+str(ne)+"\n")
                f_result.write("Profits="+str(Profits_mSGM)+"\n")
                f_result.write("S="+str(S)+"\n")
                f_result.write("numb_iter="+str(numb_iter)+"\n")
                f_result.write("numb_back="+str(numb_back)+"\n")
                f_result.write("cpu_time="+str(cpu_time)+"\n")
                f_result.write("max_iter="+str(max_iter)+"\n\n")
                f_result.close()

                # RUN SGM
                # restrict to 3600 seconds use command line
                max_iter=50
                ne, Profits_SGM,S,numb_iter,cpu_time_not_dfs  = IterativeSG_NOT_DFS(G,max_iter,1,S_ini)
                # support size for 2 player game
                size_supp_ne = [sum(1 for s in ne[:len(S[0])] if s>10**-4),sum(1 for s in ne[len(S[0]):len(S[0])+len(S[1])] if s>10**-4)]

                # SAVE TABLE OF RESULTS
                f_KEG = open("KEG_table.txt",'a')
                if cpu_time_not_dfs>3600 or Profits_SGM==[]:
                    f_KEG.write('tl & '+str(numb_iter)+' & 0 & '+str(size_supp_ne) +" & "+str([len(S[i]) for i in range(m)])+'& \n ')
                elif numb_iter>=max_iter:
                     f_KEG.write("%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & 0 & '+str(size_supp_ne)+" & "+str([len(S[i]) for i in range(m)])+'& \n ')
                else:
                    total_time_not_dfs= total_time_not_dfs+ cpu_time_not_dfs
                    total_iter_not_dfs= total_iter_not_dfs+ numb_iter
                    total_size_0_not_dfs= total_size_0_not_dfs+ len(S[0])
                    total_size_1_not_dfs= total_size_1_not_dfs+ len(S[1])
                    if Profits_SGM[0]>10**-6:
                        U_player_A_NE_not_dfs = U_player_A_NE_not_dfs + (1-(value_p[0]/Profits_SGM[0])) # decrease by acting alone
                    if Profits_SGM[1]>10**-6:
                        U_player_B_NE_not_dfs = U_player_B_NE_not_dfs + (1-(value_p[1]/Profits_SGM[1]))
                    # SOCIAL OPTIMUM
                    Social_welfare_not_dfs = Social_welfare_not_dfs + (sum(Profits_SGM)/Social_OPT) # percentage of social welfare decrease in the game
                    aux_not_dfs= aux_not_dfs+1
                    if sum(1 for s in ne if s>=1-10**-6) <m: #mixed NE
                        f_KEG.write("%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & 0 & '+str(size_supp_ne) +" & "+str([len(S[i]) for i in range(m)])+'& \n ')
                        numb_mNE_not_dfs = numb_mNE_not_dfs+1
                    else:
                        f_KEG.write("%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & 1 & 0 &'+str([len(S[i]) for i in range(m)])+'& \n ')
                        numb_pNE_not_dfs= numb_pNE_not_dfs+1
                f_KEG.close()

                # SAVE INSTANCE SOLUTIONS
                filename ="Instances/KEG/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".txt"
                f_result = open(filename, "a")
                f_result.write("SGM results \n")
                f_result.write("ne="+str(ne)+"\n")
                f_result.write("Profits="+str(Profits_SGM)+"\n")
                f_result.write("S="+str(S)+"\n")
                f_result.write("numb_iter="+str(numb_iter)+"\n")
                f_result.write("numb_back="+str(numb_back)+"\n")
                f_result.write("cpu_time="+str(cpu_time_not_dfs)+"\n")
                f_result.write("max_iter="+str(max_iter)+"\n\n")
                f_result.close()
            else:
                aux_no_game = aux_no_game +1

        # statistics utilities
        U_player_A_NE = (U_player_A_NE*1.)/aux# increase on profit due to the game
        U_player_B_NE = (U_player_B_NE*1.)/aux
        # SOCIAL OPTIMUM
        Social_welfare = (Social_welfare*1.)/aux # percentage of social welfare decrease in the game
        U_player_A_NE_not_dfs = (U_player_A_NE_not_dfs*1.)/aux_not_dfs# increase on profit due to the game
        U_player_B_NE_not_dfs = (U_player_B_NE_not_dfs*1.)/aux_not_dfs
        # SOCIAL OPTIMUM
        Social_welfare_not_dfs = (Social_welfare_not_dfs*1.)/aux_not_dfs # percentage of social welfare decrease in the game

        # statistics to DFS
        total_time = (total_time*1.)/aux
        total_iter = (total_iter*1.)/aux
        total_size_0 = (total_size_0*1.)/aux
        total_size_1 = (total_size_1*1.)/aux
        avg_numb_back = (avg_numb_back*1.)/aux
        f_KEG = open("KEG_table.txt",'a')
        f_KEG.write(str(n) +' & '+"%0.2f" % time_social+' & '+"%0.2f" % Social_welfare+' & '+"%0.2f" % U_player_A_NE+' & '+"%0.2f" % U_player_B_NE+" & "+ "%0.2f" % total_time+' & '+"%0.2f" % total_iter+' & '+"%0.2f" % total_size_0+' & '+"%0.2f" % total_size_1+' &  '+"%0.2f" % avg_numb_back+" & "+str(numb_pNE)+" & "+str((numb_mNE+numb_pNE)/(50-aux_no_game))+" & ")
        f_KEG.close()
        # statistics to NOT_DFS
        total_time_not_dfs= (total_time_not_dfs*1.)/aux_not_dfs
        total_iter_not_dfs= (total_iter_not_dfs*1.)/aux_not_dfs
        total_size_0_not_dfs= (total_size_0_not_dfs*1.)/aux_not_dfs
        total_size_1_not_dfs= (total_size_1_not_dfs*1.)/aux_not_dfs
        f_KEG = open("KEG_table.txt",'a')
        f_KEG.write("%0.2f" % Social_welfare_not_dfs+' & '+"%0.2f" % U_player_A_NE_not_dfs+' & '+"%0.2f" % U_player_B_NE_not_dfs+' & '+"%0.2f" % total_time_not_dfs+' & '+"%0.2f" % total_iter_not_dfs+" & "+"%0.2f" % total_size_0_not_dfs+' & '+"%0.2f" % total_size_1_not_dfs+' &  '+str(numb_pNE_not_dfs)+' &  '+str((numb_mNE_not_dfs+numb_pNE_not_dfs)/(50-aux_no_game))+' \\\ \n')
        f_KEG.close()
