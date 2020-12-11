from Compute_CE_maximal import *
from os import path
from KEG_SocialOPT import *

f_KEG = open("KEG_CE_table.txt",'a')
f_KEG.write(' Results for maximal CE \\\ \n')
f_KEG.close()

for m in [2]:
    #for n in [20,40,80,100]:
    for n in [20]:
        f_KEG = open("KEG_CE_table.txt",'a')
        f_KEG.write('|V| ='+str(n)+' \n\n')
        f_KEG.write(' & ins & time & iter & NE? & supp CE & size S & Decrease on Social Welfare by choosing NE instead & \\\ \n')
        f_KEG.close()

        # SGM
        total_time_not_dfs = 0
        total_iter_not_dfs = 0
        total_size_0_not_dfs = 0
        total_size_1_not_dfs = 0
        total_supp_ce = 0
        numb_ce_ne_not_dfs = 0 # number of ce that are ne
        aux_not_dfs = 0
        decrease_ce = 0
        aux_decrease = 0
        aux_no_game = 0
        #for ins in range(1,51):
        for ins in [1]:
            # read instance
            if ins!=19 or n!=20: # there is a player without strategies
                filename  ="Instances/KEG/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".npy"
                G = Game("empty")
                G.Read_Game(filename)


                # Start with the following set of strategies: anything can be selected
                S_ini = [[[]] for p in range(G.m())]
                Profile_ini = [np.array([1 for k in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())]
                for p in range(G.m()):
                    S_ini[p][0], _ ,_ = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile_ini,p,False)

                # OR USE ini of CE
                # filename ="Instances/KEG/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".txt"
                # if path.exists(filename):
                #     aux_stop =False
                #     with open(filename, "r") as f:
                #         for f_tmp in f:
                #             if f_tmp[:2]=="S=" and not aux_stop:
                #                 S_aux = eval(f_tmp.split("=")[1])
                #                 for p in range(G.m()):
                #                     S_ini[p][0]=S_aux[p][0]
                #                 aux_stop = True


                # RUN SGM
                # restrict to 3600 seconds use command line
                max_iter=1000
                ce, Profits_SGM,S,numb_iter,cpu_time_not_dfs, CE_is_NE  = IterativeSG_NOT_DFS(G,max_iter,1,S_ini)
                # support size for 2 player game
                Numb_stra = [len(S[p]) for p in range(m)]
                size_supp_ce = np.count_nonzero(ce > 10**-4)

                # SAVE TABLE OF RESULTS
                f_KEG = open("KEG_CE_table.txt",'a')
                if cpu_time_not_dfs>3600 or Profits_SGM==[]:
                    f_KEG.write(" & "+ str(ins)+ " & "+'tl & '+str(numb_iter)+' &  & '+str(size_supp_ce) +" & "+str(Numb_stra)+' &  & \\\ \n')
                elif numb_iter>=max_iter:
                     f_KEG.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' &  & '+str(size_supp_ce)+" & "+str(Numb_stra)+' & & \\\ \n')
                else:
                    total_time_not_dfs= total_time_not_dfs+ cpu_time_not_dfs
                    total_iter_not_dfs= total_iter_not_dfs+ numb_iter
                    total_size_0_not_dfs= total_size_0_not_dfs+ Numb_stra[0]
                    total_size_1_not_dfs= total_size_1_not_dfs+ Numb_stra[1]
                    total_supp_ce = total_supp_ce + size_supp_ce
                    aux_not_dfs= aux_not_dfs+1

                    # compute social increase
                    filename ="Instances/KEG/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".txt"
                    if path.exists(filename):
                        aux_stop =False
                        with open(filename, "r") as f:
                            for f_tmp in f:
                                if f_tmp[:7]=="Profits" and not aux_stop:
                                    Profits_ne = eval(f_tmp.split("=")[1])
                                    aux_stop = True
                        if sum(Profits_SGM)>0 and Profits_ne!=[]:
                            decrease = (1-(sum(Profits_ne)/sum(Profits_SGM)))
                            decrease_ce = decrease_ce+decrease
                            aux_decrease = aux_decrease +1

                            if CE_is_NE:
                                f_KEG.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & YES & '+str(size_supp_ce) +" & "+str(Numb_stra)+" & "+ "%0.2f" % decrease+' & \\\ \n')
                                numb_ce_ne_not_dfs = numb_ce_ne_not_dfs+1
                            else:
                                f_KEG.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & NO & '+str(size_supp_ce) +" & "+str(Numb_stra)+" & "+"%0.2f" % decrease+' & \\\ \n')
                        else:
                            decrease = "-" # we do not have a NE to compare
                            if CE_is_NE:
                                f_KEG.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & YES & '+str(size_supp_ce) +" & "+str(Numb_stra)+" & "+ decrease+' & \\\ \n')
                                numb_ce_ne_not_dfs = numb_ce_ne_not_dfs+1
                            else:
                                f_KEG.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & NO & '+str(size_supp_ce) +" & "+str(Numb_stra)+" & "+decrease+' & \\\ \n')
                    else:
                        decrease = "-" # we do not have a NE to compare
                        if CE_is_NE:
                            f_KEG.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & YES & '+str(size_supp_ce) +" & "+str(Numb_stra)+" & "+ decrease+' & \\\ \n')
                            numb_ce_ne_not_dfs = numb_ce_ne_not_dfs+1
                        else:
                            f_KEG.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & NO & '+str(size_supp_ce) +" & "+str(Numb_stra)+" & "+decrease+' & \\\ \n')

                f_KEG.close()

                # SAVE INSTANCE SOLUTIONS
                filename ="Instances/KEG/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".txt"
                f_result = open(filename, "a")
                f_result.write("\n Second run with good S_ini \n SGM results CE \n")
                f_result.write("ce="+str(ce)+"\n")
                f_result.write("Profits="+str(Profits_SGM)+"\n")
                f_result.write("size_supp_ce="+str(size_supp_ce)+"\n")
                f_result.write("numb_iter="+str(numb_iter)+"\n")
                f_result.write("cpu_time="+str(cpu_time_not_dfs)+"\n")
                f_result.write("NE?="+str(CE_is_NE)+"\n")
                f_result.write("max_iter="+str(max_iter)+"\n\n")
                f_result.close()
            else:
                aux_no_game = aux_no_game +1

        # statistics utilities
        total_time_not_dfs = (total_time_not_dfs*1.)/aux_not_dfs
        total_iter_not_dfs = (total_iter_not_dfs*1.)/aux_not_dfs
        total_size_0_not_dfs = (total_size_0_not_dfs*1.)/aux_not_dfs
        total_size_1_not_dfs = (total_size_1_not_dfs*1.)/aux_not_dfs
        total_supp_ce = (total_supp_ce*1.)/aux_not_dfs
        decrease_ce = (decrease_ce*1.)/aux_decrease
        f_KEG = open("KEG_CE_table.txt",'a')
        f_KEG.write(' & avg  & '+ "%0.2f" % total_time_not_dfs+' & '+"%0.2f" % total_iter_not_dfs+' & '+str(numb_ce_ne_not_dfs)+' & '+str(total_supp_ce)+" & "+"%0.2f" % total_size_0_not_dfs+' & '+"%0.2f" % total_size_1_not_dfs+' & '+"%0.2f" % decrease_ce+' \\\ \n')
        f_KEG.close()
