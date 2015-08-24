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
#boundary_intersections is a list whose i-th element is a list of intersection points that lie on the boundary of the i-th region (in the same order as per the boundary orientation of that region), so that the first two points are on some alpha circle.
#(The current implementation requires that each alpha circle contains at least 3 intersections, and each beta circle contains at least 3 intersections.)
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
##############
S234.boundary_intersections#List whose i-th element is a list of intersection points that lie on the boundary of the i-th region (in the same order as per the boundary orientation of that region), so that the first two points are on some alpha circle.
S234.regions#The regions.
S234.regions_un#The unpointed regions.
S234.intersections#The intersection points between alpha and beta circles.
S234.is_nice#True or False, depending on whether the Heegaard diagram is nice or not.        
S234.euler_measures_2#List whose i-th element is twice the Euler measures of the i-th region.
S234.euler_measures_2_un#Vector whose i-th coordinate is the twice the Euler measure of the i-th unpointed region.
S234.boundary_mat#Double list so that boundary_mat[i][j] is the coefficient of the boundary of (boundary of the i-th region restricted to alpha circles) at the j-th intersection point.
S234.boundary_mat_un#matrix with boundary_mat_un[i][j] the coefficient of the boundary of (boundary of the i-th unpointed region restricted to alpha circles) at the j-th intersection point.
S234.point_measures_4#Double list so that point_measures_4[i][j] is 4 times the point measure of the i-th region at the j-th point.
S234.point_measures_4_un#matrix with point_measures_4_un[i][j] being 4 times the point measure of i-th unpointed region at the j-th point.
S234.alphas#The alpha circles.
S234.betas#The beta circles.
S234.intersections_on_alphas#List whose i-th element is the ordered list of intersections on alpha_i (according to the orientation of alpha_i)
S234.intersections_on_betas#Ditto
S234.intersection_incidence#List whose i-th element is [a,b,n] so that the i-th intersection point lies in alpha_a and beta_b, and has alpha.beta intersection number n.
S234.intersection_matrix#Matrix whose (a,b) entry is the intersection number between alpha_a and beta_b.
S234.region_graph_alpha#Graph whose vertices are regions, and edges are alpha arcs separating two regions. Edges labeled by (a,p), where a is the alpha circle containing that alpha arc, and p is the first (according to the orientation of a) intersection point on that alpha arc.
S234.region_graph_beta#Ditto
S234.generators#The CF-hat generators.
S234.generator_reps#List whose i-th element represents the i-th generator as a tuple of intersections, whose i-th point is an intersection point on alpha_i.
#The following variables that will not be initialized initially. Some functions will initialize them, such as generate_domains().
##############
S234.SpinC_structures#The SpinC structures (ordered arbitrarily).
S234.SpinC#List whose i-th element is the SpinC structure of the i-th generator.
S234.genSpinC#List whose i-th element is a list of generators living in the i-th SpinC structure.
S234.domains_stored#Double list so that domains_stored[i][j] stores (m,D) where D is a domain (as a vector over the unpointed regions) from the i-th generator to the j-th generator, and m is its Maslov index. If there are multiple domains (happens when H^1 not 0), will store only one (the program isn't really fit for that case). If there are none (happens when H_1 not 0, i.e., have multiple SpinC structures), stores None.
S234.abs_gr#List whose i-th element is the absolute gradings of the i-th generator; can be reset manually by set_abs_gr().
S234.gr_min#Minimum absolute grading.
S234.gr_max#Ditto
S234.grading_spread#Dictionary with grading_spread[gr] being a list of generators with abs_gr=gr
#The following variables will also not be initialized initially. Some functions will initialize them, such as generate_complexes().
##############
S234.chain_complex#This is the CF-hat chain complex. It is a list whose i-th element is the complex in the i-th SpinC structure. (Currently, only works if every positive index 1 domain contributes.)

#The following is a list of functions in the HeegaardDiagram class file.
#We are running the following functions with i=10,j=0,gr=0.
[i,j,gr]=[10,0,0]
##############
S234.generate_domains()#Initializes a bunch of variables such as domains_stored and the SpinC information.
S234.can_contribute(i,j)#Returns True if the domain from i-th generator to j-th generator has Maslov index 1 and is positive; False otherwise.
S234.does_contribute(i,j)#Returns True if the domain from i-th generator to j-th generator does contribute to the CF-hat differential; False signifies no knowledge.
S234.domain_type(i,j)#Analyses the domain from i-th generator to j-th generator (only if can_contribute(i,j)). Returns (Euler characteristic, number of boundary components).
S234.set_abs_gr(i,gr)#Sets abs_gr of the i-th generator to value gr.
try:
    S234.generate_complexes()#Initializes the variable chain_complex.
    S234.compute_homology()#Computes HF-hat. Returns a list, whose i-th element is the homology in the i-th SpinC structure, which is a dictionary whose value at gr is the dim of HF-hat in grading gr.
except:
    pass
##############
S234#Prints details about the Heegaard diagram.
S234.print_differentials()#Prints what the program knows about the differentials in CF-hat.
S234.pretty_print()#Prints the same in a pretty format. Returns a graphics. See below for further documentation.


#Sometimes there is a Z/2-action on the Heegaard diagram. The program can work with a single such action.
#Declare as HeegaardDiagram(boundary_intersections,num_pointed_regions,image_of_intersections)
#image_of_intersections is a list whose i-th element is the Z/2-action image of the i-th intersection point.
S237=HeegaardDiagram([
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
        ], 1, [3,2,1,0,4,11,10,9,12,7,6,5,8,15,14,13,16])#Heegaard diagram for Sigma(2,3,7) with a Z/2-action.

#The following is a list of extra variables and functions for the Z/2-action.
##############
S237.there_is_action#True if there is a Z/2-action; False otherwise.
S237.image_of_intersections#List whose i-th element is the Z/2-action image of the i-th intersection point.
S237.image_of_generators#List whose i-th element is the Z/2-action image of the i-th generator.
S237.image_of_SpinC_structures#List whose i-th element is the Z/2-action image of the i-th SpinC_structure.
#The following variables and functions will not be initialized initially. Some functions will initialize them, such as generate_complexes().
##############
S237.mapping_cone_complex#This is a list whose i-th element is the mapping cone complex of Id plus the Z/2-action in the i-th SpinC structure if that SpinC structure is Z/2 equivariant; otherwise the i-th element is None.
try:
    S237.compute_mc_homology()#Computes homology of the mapping cone complex of Id plus the action. Returns a list whose i-th element is the homology of the mapping cone complex in SpinC structure i, as a dictionary. (Returns None is the i-th SpinC structure is not Z/2-equivariant.)
except:
    pass

#The program can also compute double branched covers along knots, returning a Heegaard diagram of the double branched cover along with the Z/2-action.
#(Currently the program requires the knot to be Z/2-nullhomologous.)
#Run the function as branched_double(H,num_pointed_regions)
#H is the starting Heegaard diagram with exactly two pointed regions.
#num_pointed_regions is the number of basepoints on the output Heegaard diagram. If one, puts it on the preimage of the last region of H; if two, puts them on the preimages of both the pointed regions of H.
trefoil=HeegaardDiagram([[0,1,2,0,2,1,0,2],[1,0],[1,2]],2)
L31=branched_double(trefoil,1)


##############
#The pretty_print function tries to display the CF-hat chain complex
#in human-readable format. The red arrows indicate that the program
#knows that the differential contributes; blue signifies no knowledge.
#Accepts optional parameters output_file, interchanges, cancellations,
#relevant_SpinC. If output_file specified, the program stores the
#graphics in that output file (default is None, when the program
#displays it). The program only shows the complex in SpinC structure
#which are in the list relevant_SpinC (default is show all SpinC
#structures).
S237.pretty_print(output_file="S237_0.pdf")
#You may make the diagram prettier via interchanges, which is a list
#of tuples of generators to be swapped (the interchanges will be done
#in the order provided).
S237.pretty_print(output_file="S237_1.pdf",interchanges=[(12,9),(13,1),(16,5),(14,0),(17,4),(15,6),(6,18),(6,0),(15,3),(15,14)])
#You may also simplify the complex by cancellations, which is a list
#of tuples of generators where the differential coefficient is known
#to be 1, and it performs the cancellations in the order
#provided. (Cancellations can produce new zig-zag arrows, which are
#shown in green.) 
S237.pretty_print(output_file="S237_2.pdf",interchanges=[(12,9),(13,1),(16,5),(14,0),(17,4),(15,6),(6,18),(6,0),(15,3),(15,14)],cancellations=[(6,4),(18,23),(0,11),(15,17),(3,21),(14,10),(5,1),(16,13)])
#The program tries to print Z/2-equivariantly if
#Z/2-action given, and relevant_SpinC is fixed under the
#Z/2-action. So tweak and cancel equivariantly in that case if you
#want the Z/2 symmetry to remain.
        
