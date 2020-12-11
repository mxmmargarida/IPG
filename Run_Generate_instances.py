# Generate problem instances
from Compute_NE import *

# np.random.seed(1)
#
# for m in [2]:
#     for n in [20,40,80,100]:
#         for ins in range(10):
#             G_KP = Game("Knapsack",m,n,ins)
#             G_KP.Save_Game(m,n,ins)
#
# for m in [3]:
#     for n in [10,20,40]:
#         for ins in range(10):
#             G_KP = Game("Knapsack",m,n,ins)
#             G_KP.Save_Game(m,n,ins)
#
# for m in [2,3]:
#     for T in [10,20,50,100]:
#         for ins in range(10):
#             G_LS = Game("LotSizing",m,T,ins)
#             G_LS.Save_Game(m,T,ins)

for m in [2]:
    for n in [100]:
        for ins in range(1,51):
            # max cycles of length 3
            G_KEG = Game("KEG",m,n,ins,50,3)
            if G_KEG.m() != 500:
                G_KEG.Save_Game(m,n,ins)
