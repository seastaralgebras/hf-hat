## Copyright (C) 2015 Sucharit Sarkar.
## Contact: sucharit@math.princeton.edu

## This file is part of hf-hat.

## hf-hat is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## hf-hat is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with hf-hat; see COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

#######BASIC CALCULATIONS#######
genus=len(intersections)
#for now, same as the number of alpha-circles.
all_intersections=[(i,j) for i in range(genus) for j in intersections[i]]
#The j-th intersection point lies in alpha_a and beta_b where (a,b)=all_intersections[j].
regions=len(boundaries)
#number of regions
for i in range(regions):
    n=len(boundaries[i])
    if n%2==1:
        raise Exception('Odd number of points on boundary.')
    for j in range(n/2):
        [a,b]=boundaries[i][2*j:2*j+2]#these should be consecutive points on some alpha-circle.
        c=boundaries[i][2*j-1]#should be on the same beta-circle as a.
        al=all_intersections[a][0]
        if al!=all_intersections[b][0]:
            raise Exception('Not on the same alpha-circle.')
        if all_intersections[a][1]!=all_intersections[c][1]:
            raise Exception('Not on the same beta-circle')
        if abs(a-b) not in [1,len(intersections[al])-1]:
            raise Exception('Not consecutive points.')
euler_measures_2=vector(ZZ,regions,[2-len(ll)/2 for ll in boundaries])
#Euler measures of the regions, times 2.            
boundary_mat=matrix(ZZ,regions,len(all_intersections),[[-boundaries[i][0::2].count(j)+boundaries[i][1::2].count(j) for j in range(len(all_intersections))] for i in range(regions)])
#boundary_mat[i][j] is the coefficient of the boundary of (boundary of the i-th region restricted to alpha circles) at the j-th intersection point.
point_measures_4=matrix(ZZ,regions,len(all_intersections),[[boundaries[i].count(j) for j in range(len(all_intersections))] for i in range(regions)])
#point_measures[i][j] is the point measure of i-th region at the j-th point, times 4.

generators=[]
#hf-generators will be represented as a tuple of g intersection points; the i-th point will be on alpha_i (i.e., the tuple is sorted).
permutations=Arrangements(range(genus),genus).list()
for permu in permutations:
    possibilities=[[j for j in range(len(all_intersections)) if all_intersections[j]==(i,permu[i])] for i in range(genus)]
    generators+=list(cartesian_product_iterator(possibilities))
        
def find_domain(initial,final):
    #finds a domain D from the initial generator to the final generator; returns (Maslov_index(D),D)
    target_vect=vector(ZZ,len(all_intersections))
    target_vect_abs=vector(ZZ,len(all_intersections))
    for j in range(len(all_intersections)):
        target_vect[j]=final.count(j)-initial.count(j)
        target_vect_abs[j]=final.count(j)+initial.count(j)
    
    answer=boundary_mat.solve_left(target_vect)#returns ValueError on failure
    answer=vector(ZZ,regions,answer)#forcing answer to be over integers, returns TypeError on failure
    maslov=2*answer.dot_product(euler_measures_2)+answer*point_measures_4*target_vect_abs
    if maslov%4!=0:
        raise Exception('Maslov index not an integer!')
    maslov=maslov/4

    return (maslov,answer)

SpinC=[0,]+(len(generators)-1)*[None,]
#SpinC[i] will be the SpinC structure of the i-th generator (SpinC structures numbered arbitrarily starting at 0, currently set to be the SpinC structure of the 0-th generator).

domains_stored=dict()
#domains_stored[(i,j)] stores (m,D) where D is a domain from the i-th generator to the j-th generator, and m is its Maslov index. If there are multiple domains (happens when H^1 not 0), will store only one (the program isn't really fit for that case). If there are none (happens when H_1 not 0, i.e., have multiple SpinC structures), stores None.
for initial_ind,initial in enumerate(generators):
    for final_ind,final in enumerate(generators):
        S=SpinC[final_ind]
        if final_ind<initial_ind:
            if domains_stored[(final_ind,initial_ind)]==None:
                domains_stored[(initial_ind,final_ind)]=None
            else:
                (m,D)=domains_stored[(final_ind,initial_ind)]
                domains_stored[(initial_ind,final_ind)]=(-m,-D)
        elif final_ind==initial_ind:
            domains_stored[(initial_ind,final_ind)]=(0,vector(ZZ,regions))
        else:
            if S==None:
                try:
                    domains_stored[(initial_ind,final_ind)]=find_domain(initial,final)
                    SpinC[final_ind]=SpinC[initial_ind]
                except (ValueError, TypeError):#Will return an error if cannot find a domain
                    domains_stored[(initial_ind,final_ind)]=None
                    if final_ind==initial_ind+1:
                        SpinC[final_ind]=max(SpinC)+1
            elif S==SpinC[initial_ind]:
                first=next(g for g in range(len(generators)) if SpinC[g]==S)
                (m1,D1)=domains_stored[(first,initial_ind)]
                (m2,D2)=domains_stored[(first,final_ind)]
                domains_stored[(initial_ind,final_ind)]=(m2-m1,D2-D1)
            else:
                domains_stored[(initial_ind,final_ind)]=None

genSpinC=[[g for g in range(len(generators)) if SpinC[g]==S] for S in range(max(SpinC)+1)]
#genSpinC[i] is a list of generators indices living in the i-th SpinC structure.

def can_contribute(i,j):
    #checks if the domain from i-th generator to j-th generator has Maslov index 1 and is positive.
    try:
        (m,D)=domains_stored[(i,j)]
        if m==1 and [a for a in D if a<0]==[]:
            return True
        else:
            return False
    except TypeError:#when no domains between initial and final.
        return False

def domain_type(i,j):
    #Analyses the domain from i-th generator to j-th generator (only if can_contribute(i,j)). Returns (Euler characteristic, number of boundary components)
    if not can_contribute(i,j):
        raise Exception("No positive Maslov index 1 domain")
    (m,D)=domains_stored[(i,j)]
    e2=D.dot_product(euler_measures_2)#twice the Euler measure of D
    iota=m-e2#the intersection number with fat diagonal

    alpha_coefficients=[sum([(D*point_measures_4)[foo] for foo in range(len(all_intersections)) if all_intersections[foo][0]==a]) for a in range(genus)]
    trivial=alpha_coefficients.count(0)#number of trivial disks
    euler=genus-trivial-iota#Euler char of Lipshitz surface
    
    bdy_comps=genus*[None,]
    while None in bdy_comps:
        first_unknown=bdy_comps.index(None)
        try:
            bdy_comps[first_unknown]=max(bdy_comps)+1
        except TypeError:
            bdy_comps[first_unknown]=0
            
        current_circle=first_unknown
        
        reached_end_of_loop=False
        while not reached_end_of_loop:
            next_circle=next(foo for foo in range(genus) if all_intersections[generators[j][current_circle]][1]==all_intersections[generators[i][foo]][1])
            bdy_comps[next_circle]=bdy_comps[current_circle]
            current_circle=next_circle
            if current_circle==first_unknown:
                reached_end_of_loop=True
                
    return (euler,max(bdy_comps)+1-trivial)

def does_contribute(i,j):
    #returns True if somehow we know that the domain from the i-th generator to the j-th generator does contribute. (False significies no knowledge.)

    (m,D)=domains_stored[(i,j)]
    (e,b)=domain_type(i,j)

    if (e,b)==(1,1):
        return True

    return False

if not abs_gradings_manual:
    (chosen_generators,their_abs_gradings)=([gS[0] for gS in genSpinC],len(genSpinC)*[0,])
else:
    chosen_generators=[generators.index(g) for g in chosen_generators]#convert generators to their indices.
    if sorted([SpinC[g] for g in chosen_generators])!=range(len(genSpinC)):
        raise Exception("Chosen generators for absolute gradings do not contain exactly one from each SpinC structure.")
#Absolute gradings; we arbitrarily set it to 0 on one generator in each SpinC grading. If you know absolute gradings from other sources, add them in the data file.
abs_gr=len(generators)*[None,]
for other_ind,other in enumerate(generators):
    chosen_generator=next(g for g in genSpinC[SpinC[other_ind]] if g in chosen_generators)
    its_abs_grading=their_abs_gradings[chosen_generators.index(chosen_generator)]
    abs_gr[other_ind]=domains_stored[(other_ind,chosen_generator)][0]+its_abs_grading
(gr_min,gr_max)=[min(abs_gr),max(abs_gr)]
#the min and max gradings
grading_spread=dict([(gr,[g for g in range(len(generators)) if abs_gr[g]==gr]) for gr in Set(abs_gr)])
#generators in each grading

if there_is_action:
    image_of_generators=[tuple(sorted([image_of_intersections[p] for p in g])) for g in generators]
    image_of_generator_indices=[generators.index(g) for g in image_of_generators]
else:
    image_of_generators=None

#######OUTPUT#######

def print_action():
    if there_is_action:
        for i in range(len(generators)):
            print repr(i)+" under the action maps to "+repr(image_of_generator_indices[i])
    else:
        print "No action specified."
    return None

def print_differentials():
    for ind_S,S in enumerate(genSpinC):
        print "In SpinC structure "+repr(ind_S)
        for i in S:
            print repr(i)+" in grading "+repr(abs_gr[i])+" under the differential (possibly) maps to "+repr([j for j in S if can_contribute(i,j)])+" with domain types (Euler char, num of boundaries): "+repr([domain_type(i,j) for j in S if can_contribute(i,j)])
    return None


def pretty_print(output_file=None,interchanges=[],cancellations=[],SpinC_structures=range(len(genSpinC))):
    #Tries to print the chain complex in human-readable format. Only prints the complex whose SpinC structure is in the list SpinC_structures (default is all). User tweaking (via interchanges, with a list of tuples of generators to be swapped) recommended; the interchanges will be done in the order provided. Cancellations is a list of tuples of generators where the differential coefficient is known to be 1, and it performs the cancellations in the order provided. Tries to print Z/2-equivariantly if Z/2-action given. So tweak and cancel equivariantly in that case if you want the Z/2 symmetry to remain. If output_file=None, shows the graphics, else stores it.

    replacements=len(generators)*[None,]#the string to be printed instead of generator; empty string kills the generator
    for g in range(len(generators)):
        if SpinC[g] in SpinC_structures:
            replacements[g]=repr(g)
        else:
            replacements[g]=''
    arrows=dict()#the differential arrows in the chain complex
    for g in range(len(generators)):
        for h in range(len(generators)):
            if can_contribute(g,h) and SpinC[g] in SpinC_structures and SpinC[h] in SpinC_structures:
                if does_contribute(g,h):
                    arrows[(g,h)]='red'#Coefficient 1
                else:
                    arrows[(g,h)]='blue'#Coefficient unknown

    #Now the cancellations
    for (a,b) in cancellations:
        if abs_gr[a]!=abs_gr[b]+1:
            raise Exception("Cancelling arrow not in correct grading")
        if arrows[(a,b)]!='red':
            raise Exception("Not clear if the cancelling arrow is coefficient 1")
        
        maps_to_a=[foo for foo in range(len(generators)) if (foo,a) in arrows]
        maps_to_b=[foo for foo in range(len(generators)) if ((foo,b) in arrows and foo!=a)]
        maps_from_a=[foo for foo in range(len(generators)) if ((a,foo) in arrows and foo!=b)]
        maps_from_b=[foo for foo in range(len(generators)) if (b,foo) in arrows]

        for g in maps_to_b:
            if arrows[(g,b)]=='red':
                replacements[g]+='+'+replacements[a]
            else:
                replacements[g]+='+c('+replacements[a]+')'#c denotes some coefficients. All c's are not equal.
            for h in maps_from_a:
                if arrows[(g,b)]=='red' and arrows[(a,h)]=='red':
                    if (g,h) in arrows:
                        if arrows[(g,h)]=='red':
                            del arrows[(g,h)]
                        else:
                            arrows[(g,h)]='green'
                    else:
                        arrows[(g,h)]='red'
                else:
                    arrows[(g,h)]='green'#Coefficient unknown, comes from change of basis.

        for g in maps_to_a:
            del arrows[(g,a)]
        for g in maps_to_b:
            del arrows[(g,b)]
        for h in maps_from_a:
            del arrows[(a,h)]
        for h in maps_from_b:
            del arrows[(b,h)]
            
        replacements[a]=''
        replacements[b]=''
        del arrows[(a,b)]
    
    width=max([len(grading_spread[gr]) for gr in grading_spread])+1

    x_coordinates=len(generators)*[0,]
    #Initial assignment
    for gr in grading_spread:
        ll=[g for g in grading_spread[gr] if SpinC[g] in SpinC_structures]
        if there_is_action:
            fixed=[g for g in ll if image_of_generator_indices[g]==g]
            non_fixed=[g for g in ll if  image_of_generator_indices[g]!=g]
            if fixed==[]:
                distance=width/(len(ll)+1)
            else:
                distance=width/(len(non_fixed)+2)

            for g in fixed:
                x_coordinates[g]=0#there could be a collision here.
            non_fixed_reps=[]
            for g in non_fixed:
                if image_of_generator_indices[g] not in non_fixed_reps:
                    non_fixed_reps.append(g)
            for ind_g,g in enumerate(non_fixed_reps):
                x_coordinates[g]=(ind_g+1)*distance
                x_coordinates[image_of_generator_indices[g]]=-(ind_g+1)*distance

        else:
            distance=width/(len(ll)+1)

            for ind_g,g in enumerate(ll):
                x_coordinates[g]=N(-width/2+(ind_g+1)*distance)

    #Now the interchanges
    for (a,b) in interchanges:
        if abs_gr[a]!=abs_gr[b]:
            raise Exception("Trying to interchange stuff in different gradings")
        temp=x_coordinates[a]
        x_coordinates[a]=x_coordinates[b]
        x_coordinates[b]=temp


    G=Graphics()
    for gr in grading_spread:
        G+=line([(-width/2,gr),(width/2,gr)],alpha=0.1)
        if abs_gradings_manual:
            G+=text(repr(gr),(width/2+1,gr),color="blue")
        else:
            G+=text(repr(gr)+'+C',(width/2+1,gr),color="blue")
            

    for g in range(len(generators)):
        if replacements[g]!='':
            G+=text(replacements[g],(x_coordinates[g],abs_gr[g]),color="black")
            maps_from_g=[h for h in range(len(generators)) if (g,h) in arrows]
            for h in maps_from_g:
                G+=line([(x_coordinates[g],abs_gr[g]),(x_coordinates[h],abs_gr[h])],alpha=0.8,color=arrows[(g,h)])

        
    G.axes(False)

    if output_file==None:
        G.show()
    else:
        G.save(output_file)
