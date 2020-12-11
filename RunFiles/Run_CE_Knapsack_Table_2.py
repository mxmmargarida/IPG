from Compute_CE import *
from os import path

for m in [2]:
    for n in [20,40,80,100]:
    #for n in [20]:
        f_knapsack = open("CE_KP_table_2player.txt",'a')
        f_knapsack.write('m ='+str(m)+' players, n = '+str(n)+' items \n\n')
        f_knapsack.write(' & ins & time & iter & NE? & supp CE & size S & Decrease on Social Welfare by choosing NE instead & \\\ \n')
        f_knapsack.close()

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

        for ins in range(10):
        #for ins in [0]:
            # read instance
            filename  ="Instances/Knapsack/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".npy"
            G = Game("empty")
            G.Read_Game(filename)

            # RUN SGM
            # restrict to 3600 seconds use command line
            max_iter=1000
            ce, Profits_SGM,S,numb_iter,cpu_time_not_dfs, CE_is_NE  = IterativeSG_NOT_DFS(G,max_iter)
            # support size for 2 player game
            Numb_stra = [len(S[p]) for p in range(m)]
            it = np.nditer(np.ones(tuple(Numb_stra)),flags=['multi_index'])
            aux_it = [it.multi_index for _ in it]
            #size_supp_ce = sum([1 for s in aux_it if ce[s]>10**-4])
            size_supp_ce = np.count_nonzero(ce > 10**-4)

            # SAVE TABLE OF RESULTS
            f_knapsack = open("CE_KP_table_2player.txt",'a')
            if cpu_time_not_dfs>3600 or Profits_SGM==[]:
                f_knapsack.write(" & "+ str(ins)+ " & "+'tl & '+str(numb_iter)+' &  & '+str(size_supp_ce) +" & "+str(Numb_stra)+' &  & \\\ \n')
            elif numb_iter>=max_iter:
                 f_knapsack.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' &  & '+str(size_supp_ce)+" & "+str(Numb_stra)+' & & \\\ \n')
            else:
                total_time_not_dfs= total_time_not_dfs+ cpu_time_not_dfs
                total_iter_not_dfs= total_iter_not_dfs+ numb_iter
                total_size_0_not_dfs= total_size_0_not_dfs+ Numb_stra[0]
                total_size_1_not_dfs= total_size_1_not_dfs+ Numb_stra[1]
                total_supp_ce = total_supp_ce + size_supp_ce
                aux_not_dfs= aux_not_dfs+1


                # compute social increase
                filename ="Instances/Knapsack/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".txt"
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
                            f_knapsack.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & YES & '+str(size_supp_ce) +" & "+str(Numb_stra)+" & "+ "%0.2f" % decrease+' & \\\ \n')
                            numb_ce_ne_not_dfs = numb_ce_ne_not_dfs+1
                        else:
                            f_knapsack.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & NO & '+str(size_supp_ce) +" & "+str(Numb_stra)+" & "+"%0.2f" % decrease+' & \\\ \n')
                    else:
                        decrease = "-" # we do not have a NE to compare
                        if CE_is_NE:
                            f_knapsack.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & YES & '+str(size_supp_ce) +" & "+str(Numb_stra)+" & "+ decrease+' & \\\ \n')
                            numb_ce_ne_not_dfs = numb_ce_ne_not_dfs+1
                        else:
                            f_knapsack.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & NO & '+str(size_supp_ce) +" & "+str(Numb_stra)+" & "+decrease+' & \\\ \n')
                else:
                    decrease = "-" # we do not have a NE to compare
                    if CE_is_NE:
                        f_knapsack.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & YES & '+str(size_supp_ce) +" & "+str(Numb_stra)+" & "+ decrease+' & \\\ \n')
                        numb_ce_ne_not_dfs = numb_ce_ne_not_dfs+1
                    else:
                        f_knapsack.write(" & "+ str(ins)+ " & "+"%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & NO & '+str(size_supp_ce) +" & "+str(Numb_stra)+" & "+decrease+' & \\\ \n')

            f_knapsack.close()

            # SAVE INSTANCE SOLUTIONS
            filename ="Instances/Knapsack/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".txt"
            f_result = open(filename, "a")
            f_result.write("SGM results CE \n")
            f_result.write("ce="+str(ce)+"\n")
            f_result.write("Profits="+str(Profits_SGM)+"\n")
            f_result.write("S="+str(S)+"\n")
            f_result.write("numb_iter="+str(numb_iter)+"\n")
            f_result.write("cpu_time="+str(cpu_time_not_dfs)+"\n")
            f_result.write("NE?="+str(CE_is_NE)+"\n")
            f_result.write("max_iter="+str(max_iter)+"\n\n")
            f_result.close()

        # statistics
        total_time_not_dfs = (total_time_not_dfs*1.)/aux_not_dfs
        total_iter_not_dfs = (total_iter_not_dfs*1.)/aux_not_dfs
        total_size_0_not_dfs = (total_size_0_not_dfs*1.)/aux_not_dfs
        total_size_1_not_dfs = (total_size_1_not_dfs*1.)/aux_not_dfs
        total_supp_ce = (total_supp_ce*1.)/aux_not_dfs
        decrease_ce = (decrease_ce*1.)/aux_decrease
        f_knapsack = open("CE_KP_table_2player.txt",'a')
        f_knapsack.write(' & avg  & '+ "%0.2f" % total_time_not_dfs+' & '+"%0.2f" % total_iter_not_dfs+' & '+str(numb_ce_ne_not_dfs)+' & '+str(total_supp_ce)+" & "+"%0.2f" % total_size_0_not_dfs+' & '+"%0.2f" % total_size_1_not_dfs+' & '+"%0.2f" % decrease_ce+' \\\ \n')
        f_knapsack.close()
