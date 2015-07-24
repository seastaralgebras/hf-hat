load hf-hat.sage

"""
Sigma_2_3_4=HeegaardDiagram([
        [0,1,10,9],
        [9,8,6,7,4,0],
        [5,6,8,7,0,4],
        [4,3,10,5],
        [1,0,7,6],
        [2,1,6,5],
        [3,4,7,8],
        [2,3,8,9],
        [1,2,9,10,3,2,5,10]
        ], 1)

Sigma_2_3_7=HeegaardDiagram([
        [3,4,7,8,10,11],
        [8,5,13,14,9,10],
        [5,6,16,13],
        [14,15,3,2,12,9],
        [1,2,1,0,5,8],
        [0,4,9,12,6,5],
        [12,11,15,14,7,6],
        [11,10,16,15],
        [14,13,0,1,8,7],
        [2,1,2,3,11,12],
        [15,16,6,7,4,0,13,16,10,9,4,3]
        ], 1, [3,2,1,0,4,11,10,9,12,7,6,5,8,15,14,13,16])


Sigma_2_3_7_m=HeegaardDiagram([
        [1,0,15,16],
        [2,1,16,17],
        [3,2,17,18],
        [4,3,18,10,5,6],
        [5,4,6,7,10,18],
        [6,5,18,17],
        [7,6,17,16],
        [8,7,16,15],
        [9,8,15,14],
        [1,2,13,14],
        [2,3,12,13],
        [3,4,11,12],
        [4,5,10,11],
        [7,8,11,10],
        [8,9,12,11],
        [9,0,13,12],
        [0,1,14,15,0,9,14,13]
        ], 1, [0,9,8,7,6,5,4,3,2,1,18,17,16,15,14,13,12,11,10])

trefoil=HeegaardDiagram([[0,1,2,0,2,1,0,2],[1,0],[1,2]],2)
L_3_1=branched_double(trefoil,1)
"""

T_3_7_nice=HeegaardDiagram([
        [[1,2,32,0],[0,1]]+
        [[a,a+1,15-a,14-a] for a in range(3,7)]+[[7,8],]+
        [[a,a+1,a-11,a-12] for a in range(13,31)]+
        [[a,a-1,(54-a)%33,53-a] for a in range(21,27)]+[[27,26],]+
        [[2,3,12,13,1,0,20,19,31,32],]
        ], 2)
#Sigma_2_3_7_again=branched_double(T_3_7_nice)

