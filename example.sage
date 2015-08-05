#This is an example file for the hf-hat program.

#Copyright (C) 2015 Sucharit Sarkar.
#Contact: sucharit@math.princeton.edu

#Copying and distribution of this file, with or without modification,
#are permitted in any medium without royalty provided the copyright
#notice and this notice are preserved.

##############

#First you have to load the main program.
load hf-hat.sage

#Declare a Heegaard diagram as HeegaardDiagram(boundary_intersections,num_pointed_regions)
#num_pointed_regions is the number of pointed regions in the Heegaard diagram (through which domains won't pass).
#Then number the regions (both pointed and unpointed) starting at 0, so that all the pointed regions are numbered last. Also number the intersections between alpha and beta circles starting at 0 in some arbitrary order.
#boundary_intersections is a list whose i-th element is a list of intersection points that lie on the boundary of the i-th region (in the same order as per the boundary orientation of that region).
S234=HeegaardDiagram([
        [0,1,10,9],
        [9,8,6,7,4,0],
        [5,6,8,7,0,4],
        [4,3,10,5],
        [1,0,7,6],
        [2,1,6,5],
        [3,4,7,8],
        [2,3,8,9],
        [1,2,9,10,3,2,5,10]
        ], 1)#Heegaard diagram for Sigma(2,3,4).

#The following is a list of variables in the HeegaardDiagram class file.
S234.boundary_intersections#List whose i-th element is a list of intersection points that lie on the boundary of the i-th region (in the same order as per the boundary orientation of that region).
S234.regions#the regions.
S234.regions_un#the unpointed regions.
S234.intersections#the intersection points between alpha and beta circles.
S234.is_nice#True or False, depending on whether the Heegaard diagram is nice or not.        
S234.euler_measures_2#List whose i-th element is twice the Euler measures of the i-th region.
S234.euler_measures_2_un#Vector whose i-th coordinate is the twice the Euler measure of the i-th unpointed region.

        #boundary_mat[i][j] is the coefficient of the boundary of (boundary of the i-th region restricted to alpha circles) at the j-th intersection point  (unpointed data stored as a matrix).  
S234.boundary_mat#Double list so that boundary_intersections[R][0::2].count(p)+S234.boundary_intersections[R][1::2].count(p) for p in S234.intersections] for R in S234.regions]
        S234.boundary_mat_un=matrix(ZZ,len(S234.regions_un),len(S234.intersections),S234.boundary_mat[:len(S234.regions_un)])
        #point_measures_4[i][j] is the point measure of i-th region at the j-th point, times 4  (unpointed data stored as a matrix).  
        S234.point_measures_4=[[S234.boundary_intersections[R].count(p) for p in S234.intersections] for R in S234.regions]
        S234.point_measures_4_un=matrix(ZZ,len(S234.regions_un),len(S234.intersections),S234.point_measures_4[:len(S234.regions_un)])


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

T_3_7_nice=HeegaardDiagram(
    [[1,2,32,0],[0,1]]+
    [[a,a+1,14-a,15-a] for a in range(3,7)]+[[7,8],]+
    [[a,a+1,a-11,a-12] for a in range(13,31)]+
    [[a,a-1,(54-a)%33,53-a] for a in range(21,27)]+[[27,26],]+
    [[2,3,12,13,1,0,20,19,31,32],]
    , 2)
S237=branched_double(T_3_7_nice,1)

