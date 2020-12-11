from Compute_NE import *
from MaxPotential_LS import MaxPotentialFunction

for m in [2]:
    for n in [10,20,50,100]:
    #for n in [10]:
        f_LotSizing = open("LS_table.txt",'a')
        f_LotSizing.write('m ='+str(m)+' players, n = '+str(n)+' items \n\n')
        f_LotSizing.write(' & ins & time & iter & pNE & mNE & size S & numb_back & time & iter & pNE & mNE & size S \\\ \n')
        f_LotSizing.close()

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

        # POTENTIAL FUNCTION MAXIMIZATION
        total_time_pot = 0
        aux_pot = 0

        for ins in range(10):
        #for ins in [0]:
            # read instance
            filename  ="Instances/LotSizing/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".npy"
            G = Game("empty")
            G.Read_Game(filename)

            # RUN m-SGM
            # restrict to 3600 seconds use command line
            max_iter=50
            ne, Profits_mSGM,S,numb_iter,numb_back,cpu_time  = IterativeSG(G,max_iter)
            # support size for 2 player game
            size_supp_ne = [sum(1 for s in ne[:len(S[0])] if s>10**-4),sum(1 for s in ne[len(S[0]):len(S[0])+len(S[1])] if s>10**-4)]

            # SAVE TABLE OF RESULTS
            f_LotSizing = open("LS_table.txt",'a')
            if cpu_time>3600 or Profits_mSGM==[]:
                f_LotSizing.write("&"+str(ins)+' & '+'tl & '+str(numb_iter)+' & 0 & '+str(size_supp_ne) +" & "+str([len(S[i]) for i in range(m)])+' & '+str(numb_back) +' & &')
            elif numb_iter>=max_iter:
                 f_LotSizing.write("&"+str(ins)+' & '+"%0.2f" % cpu_time+' & '+str(numb_iter)+' & 0 & '+str(size_supp_ne)+" & "+str([len(S[i]) for i in range(m)])+' & '+str(numb_back) +' & &')
            else:
                total_time = total_time + cpu_time
                total_iter = total_iter + numb_iter
                total_size_0 = total_size_0 + len(S[0])
                total_size_1 = total_size_1 + len(S[1])
                avg_numb_back = avg_numb_back + numb_back
                aux = aux +1
                if sum(1 for s in ne if s>=1-10**-6) <m: #mixed NE
                    f_LotSizing.write("&"+str(ins)+' & '+"%0.2f" % cpu_time+' & '+str(numb_iter)+' & 0 & '+str(size_supp_ne) +" & "+str([len(S[i]) for i in range(m)])+' & '+str(numb_back) +' & &')
                    numb_mNE = numb_mNE+1
                else:
                    f_LotSizing.write("&"+str(ins)+' & '+"%0.2f" % cpu_time+' & '+str(numb_iter)+' & 1 & 0 &'+str([len(S[i]) for i in range(m)])+' & '+str(numb_back) +' & & ')
                    numb_pNE = numb_pNE+1
            f_LotSizing.close()

            # SAVE INSTANCE SOLUTIONS
            filename ="Instances/LotSizing/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".txt"
            f_result = open(filename, "a")
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
            ne, Profits_SGM,S,numb_iter,cpu_time_not_dfs  = IterativeSG_NOT_DFS(G,max_iter)
            # support size for 2 player game
            size_supp_ne = [sum(1 for s in ne[:len(S[0])] if s>10**-4),sum(1 for s in ne[len(S[0]):len(S[0])+len(S[1])] if s>10**-4)]

            # SAVE TABLE OF RESULTS
            f_LotSizing = open("LS_table.txt",'a')
            if cpu_time_not_dfs>3600 or Profits_SGM==[]:
                f_LotSizing.write('tl & '+str(numb_iter)+' & 0 & '+str(size_supp_ne) +" & "+str([len(S[i]) for i in range(m)])+'&  ')
            elif numb_iter>=max_iter:
                 f_LotSizing.write("%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & 0 & '+str(size_supp_ne)+" & "+str([len(S[i]) for i in range(m)])+'& ')
            else:
                total_time_not_dfs= total_time_not_dfs+ cpu_time_not_dfs
                total_iter_not_dfs= total_iter_not_dfs+ numb_iter
                total_size_0_not_dfs= total_size_0_not_dfs+ len(S[0])
                total_size_1_not_dfs= total_size_1_not_dfs+ len(S[1])
                aux_not_dfs= aux_not_dfs+1
                if sum(1 for s in ne if s>=1-10**-6) <m: #mixed NE
                    f_LotSizing.write("%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & 0 & '+str(size_supp_ne) +" & "+str([len(S[i]) for i in range(m)])+'& ')
                    numb_mNE_not_dfs = numb_mNE_not_dfs+1
                else:
                    f_LotSizing.write("%0.2f" % cpu_time_not_dfs+' & '+str(numb_iter)+' & 1 & 0 &'+str([len(S[i]) for i in range(m)])+'&  ')
                    numb_pNE_not_dfs= numb_pNE_not_dfs+1
            f_LotSizing.close()

            # SAVE INSTANCE SOLUTIONS
            filename ="Instances/LotSizing/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".txt"
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

            # Potential function maximization
            cpu_time_pot,X_pot,Profits_pot= MaxPotentialFunction(G.A_market(),G.B(),G.F(),G.H(),G.C(),G.M(),m,n)
            # SAVE INSTANCE SOLUTIONS
            filename ="Instances/LotSizing/Game_"+str(m)+"_"+str(n)+"_"+str(ins)+".txt"
            f_result = open(filename, "a")
            f_result.write("Potential results \n")
            f_result.write("X_pot="+str(X_pot)+"\n")
            f_result.write("Profits_pot="+str(Profits_pot)+"\n")
            f_result.write("cpu_time_pot="+str(cpu_time_pot)+"\n")
            f_result.close()
            # SAVE TABLE OF RESULTS
            if cpu_time_pot!="tl":
                total_time_pot = total_time_pot + cpu_time_pot
                aux_pot = aux_pot+1
                f_LotSizing = open("LS_table.txt",'a')
                f_LotSizing.write(' & '+"%0.2f" % cpu_time_pot+'  & \\\ \n')
                f_LotSizing.close()
            else:
                f_LotSizing = open("LS_table.txt",'a')
                f_LotSizing.write(' & tl  & \\\ \n')
                f_LotSizing.close()

        # statistics to DFS
        total_time = (total_time*1.)/aux
        total_iter = (total_iter*1.)/aux
        total_size_0 = (total_size_0*1.)/aux
        total_size_1 = (total_size_1*1.)/aux
        avg_numb_back = (avg_numb_back*1.)/aux
        f_LotSizing = open("LS_table.txt",'a')
        f_LotSizing.write(' & '+ str(n) +' & '+ "%0.2f" % total_time+' & '+"%0.2f" % total_iter+' & '+"%0.2f" % total_size_0+' & '+"%0.2f" % total_size_1+' & & '+"%0.2f" % avg_numb_back+" & "+str(numb_pNE)+" & ")
        f_LotSizing.close()
        # statistics to NOT_DFS
        total_time_not_dfs= (total_time_not_dfs*1.)/aux_not_dfs
        total_iter_not_dfs= (total_iter_not_dfs*1.)/aux_not_dfs
        total_size_0_not_dfs= (total_size_0_not_dfs*1.)/aux_not_dfs
        total_size_1_not_dfs= (total_size_1_not_dfs*1.)/aux_not_dfs
        f_LotSizing = open("LS_table.txt",'a')
        f_LotSizing.write("%0.2f" % total_time_not_dfs+' & '+"%0.2f" % total_iter_not_dfs+" & "+"%0.2f" % total_size_0_not_dfs+' & '+"%0.2f" % total_size_1_not_dfs+' & & '+str(numb_pNE_not_dfs)+' & ')
        f_LotSizing.close()

        # statistics to POTENTIAL FUNCTION MAXIMIZATION
        total_time_pot = (total_time_pot*1.)/aux_pot
        f_LotSizing = open("LS_table.txt",'a')
        f_LotSizing.write("%0.2f" % total_time_pot+' & '+str(aux_pot)+'\\\ \n')
        f_LotSizing.close()
