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


class HeegaardDiagram():

    def __init__(self,boundary_intersections,num_pointed_regions,image_of_intersections=None):
        """
        Number all the regions (complement of alpha,beta), starting at
        0, so that all the unpointed regions come first. Number all
        the intersection points between alpha,beta, again starting at
        0 (the numbering is arbitrary). 
        
        The boundary_intersections is a list whose i-th element is a list of
        intersection points that appear on the boundary of the i-th
        region, according to the boundary orientation, so that the
        part of the boundary joining the first two points lie on an
        alpha circle. (We are implicitly assuming that all regions are
        planar with one boundary (which can intersect itself).)

        The num_pointed_regions is the number of pointed regions
        (always the last few). 

        Sometimes, there is a Z/2-action. The program can work with a
        single such Z/2-action. If there is a Z/2-action, input
        image_of_intersections as list whose i-th element is the
        Z/2-action image of the i-th intersection point. (The default
        is None.)
        """

        #Basic initialization. (Will convert boundary_intersections to its lexicographically smallest form.)
        self.boundary_intersections=[]
        for intersections_list in boundary_intersections:
            smallest_indices=[2*ind_i for ind_i,i in enumerate(intersections_list[0::2]) if i==min(intersections_list[0::2])]
            next_intersections=[intersections_list[ind_i+1] for ind_i in smallest_indices]
            smallest_index=next(ind_i for ind_i in smallest_indices if intersections_list[ind_i+1]==min(next_intersections))
            self.boundary_intersections.append(intersections_list[smallest_index:]+intersections_list[:smallest_index])
        if image_of_intersections!=None:
            self.there_is_action=True
            self.image_of_intersections=image_of_intersections
        else:
            self.there_is_action=False

        self.regions=range(len(self.boundary_intersections))#the regions
        self.regions_un=range(len(self.boundary_intersections)-num_pointed_regions)#the unpointed regions
        self.intersections=range(1+max(flatten(self.boundary_intersections)))#the intersection points

        #Euler measures of the regions, times 2 (unpointed data stored as a matrix).            
        self.euler_measures_2=[2-len(self.boundary_intersections[R])/2 for R in self.regions]
        self.euler_measures_2_un=vector(ZZ,len(self.regions_un),self.euler_measures_2[:len(self.regions_un)])
        if min(self.euler_measures_2_un)>-1:
            self.is_nice=True#is a nice diagram
        else:
            self.is_nice=False
        #boundary_mat[i][j] is the coefficient of the boundary of (boundary of the i-th region restricted to alpha circles) at the j-th intersection point  (unpointed data stored as a matrix).  
        self.boundary_mat=[[-self.boundary_intersections[R][0::2].count(p)+self.boundary_intersections[R][1::2].count(p) for p in self.intersections] for R in self.regions]
        self.boundary_mat_un=matrix(ZZ,len(self.regions_un),len(self.intersections),self.boundary_mat[:len(self.regions_un)])
        #point_measures_4[i][j] is the point measure of i-th region at the j-th point, times 4  (unpointed data stored as a matrix).  
        self.point_measures_4=[[self.boundary_intersections[R].count(p) for p in self.intersections] for R in self.regions]
        self.point_measures_4_un=matrix(ZZ,len(self.regions_un),len(self.intersections),self.point_measures_4[:len(self.regions_un)])

        #intersections_on_alphas[i] is the ordered list of intersections on alpha_i (according to the orientation of alpha_i). Similarly, for beta. 
        self.intersections_on_alphas=[]
        while len(flatten(self.intersections_on_alphas))<len(self.intersections):
            start_p=next(p for p in self.intersections if p not in flatten(self.intersections_on_alphas))
            new_circle=[]
            curr_p=start_p
            while curr_p!=start_p or new_circle==[]:
                new_circle.append(curr_p)
                found_next_point=False
                regions_with_curr_point=[(R,ind_p) for R in self.regions for ind_p in range(len(self.boundary_intersections[R])) if curr_p == self.boundary_intersections[R][ind_p]]
                for (R,ind_p) in regions_with_curr_point:
                    if not found_next_point:
                        if ind_p%2==0 and self.boundary_intersections[R][ind_p+1] not in new_circle:
                            found_next_point=True
                            curr_p=self.boundary_intersections[R][ind_p+1]
                        elif ind_p%2==1 and self.boundary_intersections[R][ind_p-1] not in new_circle:
                            found_next_point=True
                            curr_p=self.boundary_intersections[R][ind_p-1]
                if not found_next_point:#must have completed the cycle
                    curr_p=start_p
            self.intersections_on_alphas.append(new_circle)
        self.intersections_on_betas=[]
        while len(flatten(self.intersections_on_betas))<len(self.intersections):
            start_p=next(p for p in self.intersections if p not in flatten(self.intersections_on_betas))
            new_circle=[]
            curr_p=start_p
            while curr_p!=start_p or new_circle==[]:
                new_circle.append(curr_p)
                found_next_point=False
                regions_with_curr_point=[(R,ind_p) for R in self.regions for ind_p in range(len(self.boundary_intersections[R])) if curr_p == self.boundary_intersections[R][ind_p]]
                for (R,ind_p) in regions_with_curr_point:
                    if not found_next_point:
                        if ind_p%2==0 and self.boundary_intersections[R][ind_p-1] not in new_circle:
                            found_next_point=True
                            curr_p=self.boundary_intersections[R][ind_p-1]
                        elif ind_p%2==1 and self.boundary_intersections[R][ind_p+1-len(self.boundary_intersections[R])] not in new_circle:
                            found_next_point=True
                            curr_p=self.boundary_intersections[R][ind_p+1-len(self.boundary_intersections[R])]
                if not found_next_point:#must have completed the cycle
                    curr_p=start_p
            self.intersections_on_betas.append(new_circle)

        self.alphas=range(len(self.intersections_on_alphas))#the alpha circles
        self.betas=range(len(self.intersections_on_betas))#the beta circles
        if sorted(flatten(self.intersections_on_alphas))!=self.intersections:
            raise Exception("Alpha circles don't contain all the intersections")
        if sorted(flatten(self.intersections_on_betas))!=self.intersections:
            raise Exception("Beta circles don't contain all the intersections")

        self.intersection_incidence=[[next(a for a in self.alphas if p in self.intersections_on_alphas[a]),next(b for b in self.betas if p in self.intersections_on_betas[b])] for p in self.intersections]#the i-th intersection point lies in alpha_a and beta_b, and has alpha.beta intersection number n, where [a,b,n]=intersection_incidence[i]
        for p in self.intersections:
            try:
                [a,b]=self.intersection_incidence[p]

                if len(self.intersections_on_alphas[a])>2 and len(self.intersections_on_betas[b])>2:#if any of a or b has 2 or fewer intersections, their orientations are not yet well-defined.
                    [a_pind,b_pind]=[(self.intersections_on_alphas[a]).index(p),(self.intersections_on_betas[b]).index(p)]
                    (R,ind_p)=next((R,ind_p) for R in self.regions for ind_p in range(len(self.boundary_intersections[R])/2) if self.boundary_intersections[R][2*ind_p]==p)
                    prev_b=self.boundary_intersections[R][2*ind_p-1]
                    next_a=self.boundary_intersections[R][2*ind_p+1]
                    intersection_number=-1#if a were oriented from p to next_a, and b oriented from prev_b to p, then the intersection.
                    if self.intersections_on_alphas[a][a_pind-1]==next_a:#a is actually oriented oppositely
                        intersection_number=-intersection_number
                    if self.intersections_on_betas[b][b_pind-1]!=prev_b:#b is actually oriented oppositely
                        intersection_number=-intersection_number
                    (self.intersection_incidence[p]).append(intersection_number)
                    
                elif len(self.intersections_on_alphas[a])==1 or len(self.intersections_on_betas[b])==1:
                    (self.intersection_incidence[p]).append(1)
                else:
                    print "WARNING: The current implementation of intersection numbers doesn't work if some alpha or beta circle has exactly 2 intersections."
                    raise Exception("Couldn't orient the alpha and beta circles")#Comment out exception depending on how we feel.
            except ValueError:#intersection number at p already determined
                pass
        self.intersection_matrix=matrix(ZZ,len(self.alphas),len(self.betas))#the (a,b) entry is the intersection number between alpha circle a and beta circle b
        for a in self.alphas:
            for b in self.betas:
                self.intersection_matrix[a,b]=sum([n for [i,j,n] in self.intersection_incidence if (a==i and b==j)])

        #region_graph_alpha is a graph whose vertices are regions, and edges are alpha arcs separating two regions. Edges labeled by which alpha circle, and the alpha arc (as the first (according to the orientation of that alpha circle) intersection point on that alpha arc). Ditto for beta circles
        self.region_graph_alpha=Graph(len(self.regions),multiedges=True,loops=True)
        for a in self.alphas:
            for ind_p,p in enumerate(self.intersections_on_alphas[a]):
                q=self.intersections_on_alphas[a][ind_p+1-len(self.intersections_on_alphas[a])]
                regions_with_pq=[R for R in self.regions for ind_foo,foo in enumerate(self.boundary_intersections[R][0::2]) if sorted([foo,self.boundary_intersections[R][2*ind_foo+1]])==sorted([p,q])]
                if len(regions_with_pq)!=2:
                    raise Exception("Each alpha arc has two adjacent regions")
                self.region_graph_alpha.add_edge(regions_with_pq+[(a,p),])
        self.region_graph_beta=Graph(len(self.regions),multiedges=True,loops=True)
        for b in self.betas:
            for ind_p,p in enumerate(self.intersections_on_betas[b]):
                q=self.intersections_on_betas[b][ind_p+1-len(self.intersections_on_betas[b])]
                regions_with_pq=[R for R in self.regions for ind_foo,foo in enumerate(self.boundary_intersections[R][0::2]) if sorted([foo,self.boundary_intersections[R][2*ind_foo-1]])==sorted([p,q])]
                if len(regions_with_pq)!=2:
                    raise Exception("Each beta arc has two adjacent regions")
                self.region_graph_beta.add_edge(regions_with_pq+[(b,p),])

        #Now some more error checking. (Definitely not a complete collection.) If diagram is wrong, some error could also be raised by other parts of the program.
        if vector(ZZ,len(self.regions),(len(self.regions))*[1,])*matrix(ZZ,len(self.regions),len(self.intersections),self.point_measures_4)!=vector(ZZ,len(self.intersections),(len(self.intersections))*[4,]):
            raise Exception("Every intersection should have four corners.")
        if len(self.alphas)!=len(self.betas):
            raise Exception("Require same number of alpha and beta circles.")
        if self.there_is_action:
            for p in self.intersections:
                if self.image_of_intersections[self.image_of_intersections[p]]!=p:
                    raise Exception("Not an involution.")
            for R in self.regions:
                image_of_R=[self.image_of_intersections[p] for p in self.boundary_intersections[R]]
                smallest_indices=[2*ind_i for ind_i,i in enumerate(image_of_R[0::2]) if i==min(image_of_R[0::2])]
                next_intersections=[image_of_R[ind_i+1] for ind_i in smallest_indices]
                smallest_index=next(ind_i for ind_i in smallest_indices if image_of_R[ind_i+1]==min(next_intersections))
                image_of_R=image_of_R[smallest_index:]+image_of_R[:smallest_index]
                if image_of_R not in self.boundary_intersections:
                    raise Exception("The Z/2-action is not a valid action.")

        #generators is the list of hf-generators, and generator_reps are their representatives. Each generator is represented as a tuple of intersections; the i-th point will be on alpha_i.
        self.generator_reps=[]
        permutations=Arrangements(self.alphas,len(self.alphas)).list()
        for permu in permutations:
            possibilities=[[p for p in self.intersections_on_alphas[a] if self.intersection_incidence[p][1]==permu[a]] for a in self.alphas]
            self.generator_reps+=list(cartesian_product_iterator(possibilities))
        self.generators=range(len(self.generator_reps))#the generator names

        if self.there_is_action:
            self.image_of_generators=[]#images of the generators under the Z/2-action
            for g in self.generators:
                image=[self.image_of_intersections[self.generator_reps[g][a]] for a in self.alphas]
                image=tuple([next(p for p in image if self.intersection_incidence[p][0]==a) for a in self.alphas])
                self.image_of_generators.append(self.generator_reps.index(image))


        #Some more variables that will not be initialized initially, since we need to solve matrix equations.
        self.SpinC="Not yet initialized"#SpinC[i] will be the SpinC structure of the i-th generator (numbered arbitrarily starting at 0).
        self.SpinC_structures="Not yet initialized"#the SpinC structures
        self.genSpinC="Not yet initialized"#genSpinC[i] is a list of generators indices living in the i-th SpinC structure.
        #domains_stored[i][j] stores (m,D) where D is a domain (as a vector over the unpointed regions) from the i-th generator to the j-th generator, and m is its Maslov index. If there are multiple domains (happens when H^1 not 0), will store only one (the program isn't really fit for that case). If there are none (happens when H_1 not 0, i.e., have multiple SpinC structures), stores None.
        if self.there_is_action:
            self.image_of_SpinC_structures="Not yet initialized"#if_there_is_action, image of SpinC_structures under the action.
        self.domains_stored="Not yet initialized"#domains_stored will be a double list, domains_stored[i][j] will store (m,D) or None; where D is some domain from i to j and m is its Maslov index.
        self.abs_gr="Not yet initialized"#The absolute gradings of the generators; can be reset manually by set_abs_gr.
        self.gr_min="Not yet initialized"
        self.gr_max="Not yet initialized"#the minimum and maximum abs grading
        self.grading_spread="Not yet initialized"#grading_spread[gr] is a list of generators with abs_gr=gr
        self.fully_initialized=False#will be set to true once the above variables are computed.

        #The chain complexes will also not be initialized initially.
        self.chain_complex="Not yet initialized"#This is the HF-hat chain complex. It will be a list whose i-th element is the complex in the i-th SpinC structure.
        if self.there_is_action:
            self.mapping_cone_complex="Not yet initialized"#This is the mapping cone complex of the map Id plus the action. The i-th element is the complex in the i-th SpinC structure if that SpinC structure is Z/2 equivariant, else None.
        self.chain_complexes_initialized=False#will be set to true once the above two variables are initialized.
        

    def __repr__(self):
        
        return "Heegaard diagram with "+repr(len(self.alphas))+" alpha and beta circles intersecting each other in "+repr(len(self.intersections))+" intersection points. The intersection points on the alpha circles are "+repr(self.intersections_on_alphas)+" and the intersection points on the beta circles are "+repr(self.intersections_on_betas)+". There are "+repr(len(self.regions_un))+" unpointed regions, and the intersection points appearing on their boundaries are "+repr(self.boundary_intersections[:len(self.regions_un)])+"."

    def find_domain(self,initial,final):
        #finds a domain D from the initial generator to the final generator; returns (Maslov_index(D),D). Called by generate_domain.
        target_vect=vector(ZZ,len(self.intersections))
        target_vect_abs=vector(ZZ,len(self.intersections))
        for p in self.intersections:
            target_vect[p]=(self.generator_reps[final]).count(p)-(self.generator_reps[initial]).count(p)
            target_vect_abs[p]=(self.generator_reps[final]).count(p)+(self.generator_reps[initial]).count(p)

        answer=self.boundary_mat_un.solve_left(target_vect)#returns ValueError on failure
        answer=vector(ZZ,len(self.regions_un),answer)#forcing answer to be over integers, returns TypeError on failure
        maslov=2*answer.dot_product(self.euler_measures_2_un)+answer*self.point_measures_4_un*target_vect_abs
        if maslov%4!=0:
            raise Exception('Maslov index not an integer!')
        maslov=maslov/4

        return (maslov,answer)

    def generate_domains(self):
        #Fills up domains_stored, and consequently, all the SpinC information.

        if self.fully_initialized:
            return True

        self.domains_stored=[]
        self.SpinC=[0,]+(len(self.generators)-1)*[None,]
        self.abs_gr=(len(self.generators))*[None,]
        self.grading_spread=dict()

        for initial in self.generators:
            self.domains_stored.append([])
            for final in self.generators:
                S=self.SpinC[final]
                if final<initial:
                    if self.domains_stored[final][initial]!=None:
                        (m,D)=self.domains_stored[final][initial]
                        (self.domains_stored[initial]).append((-m,-D))
                    else:
                        (self.domains_stored[initial]).append(None)
                elif final==initial:
                    (self.domains_stored[initial]).append((0,vector(ZZ,len(self.regions_un))))
                else:
                    if S==None:
                        try:
                            (m,D)=self.find_domain(initial,final)
                            (self.domains_stored[initial]).append((m,D))
                            self.SpinC[final]=self.SpinC[initial]
                        except (ValueError, TypeError):#Will return an error if cannot find a domain
                            if final==initial+1:
                                self.SpinC[final]=max(self.SpinC)+1
                            (self.domains_stored[initial]).append(None)
                    elif S==self.SpinC[initial]:
                        first=next(g for g in self.generators if self.SpinC[g]==S)
                        (m1,D1)=self.domains_stored[first][initial]
                        (m2,D2)=self.domains_stored[first][final]
                        (self.domains_stored[initial]).append((m2-m1,D2-D1))
                    else:
                        (self.domains_stored[initial]).append(None)

        self.SpinC_structures=range(max(self.SpinC)+1)
        self.genSpinC=[[g for g in self.generators if self.SpinC[g]==S] for S in self.SpinC_structures]
        if self.there_is_action:
            self.image_of_SpinC_structures=[self.SpinC[self.image_of_generators[self.genSpinC[S][0]]] for S in self.SpinC_structures]

        for S in self.SpinC_structures:
            base_gen=self.genSpinC[S][0]
            self.abs_gr[base_gen]=0
            if self.there_is_action:
                if self.image_of_SpinC_structures[S]<S:
                    self.abs_gr[base_gen]=self.abs_gr[self.image_of_generators[base_gen]]
                elif self.image_of_SpinC_structures[S]==S:#Error check.
                    if self.domains_stored[self.image_of_generators[base_gen]][base_gen][0]!=0:
                        raise Exception("The action does not respect gradings.")
            for g in self.genSpinC[S]:
                self.abs_gr[g]=self.abs_gr[base_gen]+self.domains_stored[g][base_gen][0]
        self.gr_min=min(self.abs_gr)
        self.gr_max=max(self.abs_gr)
        for gr in Set(self.abs_gr):
            self.grading_spread[gr]=[g for g in self.generators if self.abs_gr[g]==gr]

        self.fully_initialized=True

    def set_abs_gr(self,base_gen,gr):
        #Sets abs_gr of generator base_gen to value gr.
        if not self.fully_initialized:
            self.generate_domains()

        change=gr-self.abs_gr[base_gen]
        S=self.SpinC[base_gen]
        if self.there_is_action:
            T=self.image_of_SpinC_structures[S]
        if self.there_is_action and T!=S:
            generators_to_change=self.genSpinC[S]+self.genSpinC[T]
        else:
            generators_to_change=self.genSpinC[S]
        for g in generators_to_change:
            self.abs_gr[g]+=change

        self.gr_min=min(self.abs_gr)
        self.gr_max=max(self.abs_gr)
        self.grading_spread=dict()
        for gr in Set(self.abs_gr):
            self.grading_spread[gr]=[g for g in self.generators if self.abs_gr[g]==gr]
        

        if self.chain_complexes_initialized:
            self.chain_complex[S]=dict([(gr+change,self.chain_complex[S][gr]) for gr in self.chain_complex[S]])
            if self.there_is_action:
                if T!=S:
                    self.chain_complex[T]=dict([(gr+change,self.chain_complex[T][gr]) for gr in self.chain_complex[T]])
                else:
                    self.mapping_cone_complex[S]=dict([(gr+change,self.mapping_cone_complex[S][gr]) for gr in self.mapping_cone_complex[S]])
                    


    def can_contribute(self,i,j):
        #checks if the domain from i-th generator to j-th generator has Maslov index 1 and is positive.
        if not self.fully_initialized:
            self.generate_domains()
        try:
            (m,D)=self.domains_stored[i][j]
            if m==1 and [a for a in D if a<0]==[]:
                return True
            else:
                return False
        except TypeError:#when no domains between initial and final.
            return False

    def domain_type(self,i,j):
        #Analyses the domain from i-th generator to j-th generator (only if can_contribute(i,j)). Returns (Euler characteristic, number of boundary components)
        if not self.can_contribute(i,j):
            raise Exception("No positive Maslov index 1 domain")
        if self.is_nice:
            return (1,1)

        (m,D)=self.domains_stored[i][j]
        e2=D.dot_product(self.euler_measures_2_un)#twice the Euler measure of D
        iota=m-e2#the intersection number with fat diagonal

        alpha_coefficients=[sum([(D*self.point_measures_4_un)[p] for p in self.intersections_on_alphas[a]]) for a in self.alphas]
        num_trivial=alpha_coefficients.count(0)#number of trivial disks
        euler=len(self.alphas)-num_trivial-iota#Euler char of Lipshitz surface

        bdy_comps=len(self.alphas)*[None,]
        while None in bdy_comps:
            first_alpha=bdy_comps.index(None)
            try:
                bdy_comps[first_alpha]=max(bdy_comps)+1
            except TypeError:
                bdy_comps[first_alpha]=0

            current_alpha=first_alpha

            reached_end_of_loop=False
            while not reached_end_of_loop:
                next_alpha=next(a for a in self.alphas if self.intersection_incidence[self.generator_reps[j][current_alpha]][1]==self.intersection_incidence[self.generator_reps[i][a]][1])
                bdy_comps[next_alpha]=bdy_comps[current_alpha]
                current_alpha=next_alpha
                if current_alpha==first_alpha:
                    reached_end_of_loop=True

        return (euler,max(bdy_comps)+1-num_trivial)

    def does_contribute(self,i,j):
        #returns True if somehow we know that the domain from the i-th generator to the j-th generator does contribute. (False signifies no knowledge.)
        (e,b)=self.domain_type(i,j)
        (m,D)=self.domains_stored[i][j]

        if (e,b)==(1,1):
            return True
        #Add more, if possible.

        return False


    def generate_complexes(self):
        if self.chain_complexes_initialized:
            return True

        if not self.fully_initialized:
            self.generate_domains()

        for g in self.generators:
            for h in self.generators:
                if self.can_contribute(g,h) and (not self.does_contribute(g,h)):
                    raise Exception("Not enough data to generate complexes")

        self.chain_complex=[]
        if self.there_is_action:
            self.mapping_cone_complex=[]

        for S in self.SpinC_structures:
            chain_complex=dict()#This will be chain complex in SpinC structure S. chain_complex[gr] will be the matrix from grading gr to gr-1.
            gradings_in_S=sorted(list(Set([self.abs_gr[g] for g in self.genSpinC[S]])))#the gradings that can appear in this SpinC structure.
            for gr in gradings_in_S:
                if (gr-1) in gradings_in_S:
                    generators_higher=[g for g in self.genSpinC[S] if self.abs_gr[g]==gr]
                    generators_lower=[g for g in self.genSpinC[S] if self.abs_gr[g]==gr-1]

                    boundary_matrix=matrix(GF(2),len(generators_higher),len(generators_lower))
                    for ind_g,g in enumerate(generators_higher):
                        for ind_h,h in enumerate(generators_lower):
                            if self.can_contribute(g,h):
                                boundary_matrix[ind_g,ind_h]=1
                    chain_complex[gr]=boundary_matrix
            (self.chain_complex).append(chain_complex)

        if self.there_is_action:
            for S in self.SpinC_structures:
                if self.image_of_SpinC_structures[S]!=S:#Not a fixed SpinC structure. Won't generate the mapping cone complex.
                    (self.mapping_cone_complex).append(None)
                else:#Will add the mapping cone complex of Id + the action.
                    mapping_cone_complex=dict()
                    gradings_in_S=sorted(list(Set([self.abs_gr[g] for g in self.genSpinC[S]]+[self.abs_gr[g]+1 for g in self.genSpinC[S]])))#the gradings that can appear in this SpinC structure.
                    for gr in gradings_in_S:
                        if (gr-1) in gradings_in_S:
                            generators_higher=[g for g in self.genSpinC[S] if self.abs_gr[g]==gr]+[g for g in self.genSpinC[S] if self.abs_gr[g]==gr-1]
                            generators_lower=[g for g in self.genSpinC[S] if self.abs_gr[g]==gr-1]+[g for g in self.genSpinC[S] if self.abs_gr[g]==gr-2]

                            mc_boundary_matrix=matrix(GF(2),len(generators_higher),len(generators_lower))
                            for ind_g,g in enumerate(generators_higher):
                                for ind_h,h in enumerate(generators_lower):
                                    if self.can_contribute(g,h):
                                        mc_boundary_matrix[ind_g,ind_h]=1
                                    if self.image_of_generators[g]!=g and (g==h or self.image_of_generators[g]==h):
                                        mc_boundary_matrix[ind_g,ind_h]=1
                            mapping_cone_complex[gr]=mc_boundary_matrix
                    (self.mapping_cone_complex).append(mapping_cone_complex)

        self.chain_complexes_initialized=True


    def compute_homology(self):
        #Computes homology of the chain complex
        if not self.chain_complexes_initialized:
            self.generate_complexes()

        homology=[]

        for S in self.SpinC_structures:
            homology_in_S=dict()
            ranks_in_S=dict()
            gradings_in_S=sorted(list(Set([self.abs_gr[g] for g in self.genSpinC[S]])))
            for gr in self.chain_complex[S]:
                ranks_in_S[gr]=(self.chain_complex[S][gr]).rank()
            for gr in gradings_in_S:
                homology_in_S[gr]=len([g for g in self.genSpinC[S] if self.abs_gr[g]==gr])
                if (gr-1) in gradings_in_S:
                    homology_in_S[gr]-=ranks_in_S[gr]
                if (gr+1) in gradings_in_S:
                    homology_in_S[gr]-=ranks_in_S[gr+1]
            homology.append(dict([(gr,homology_in_S[gr]) for gr in homology_in_S if homology_in_S[gr]!=0]))

        return homology

    def compute_mc_homology(self):
        #Computes homology of the mapping cone complex of Id + the action
        if not self.there_is_action:
            raise Exception("No map to have a mapping cone complex")

        if not self.chain_complexes_initialized:
            self.generate_complexes()

        mc_homology=[]

        for S in self.SpinC_structures:
            if self.mapping_cone_complex[S]==None:
                mc_homology.append(None)
            else:
                mc_homology_in_S=dict()
                ranks_in_S=dict()
                gradings_in_S=sorted(list(Set([self.abs_gr[g] for g in self.genSpinC[S]]+[self.abs_gr[g]+1 for g in self.genSpinC[S]])))
                for gr in self.mapping_cone_complex[S]:
                    ranks_in_S[gr]=(self.mapping_cone_complex[S][gr]).rank()
                for gr in gradings_in_S:
                    mc_homology_in_S[gr]=len([g for g in self.genSpinC[S] if self.abs_gr[g]==gr]+[g for g in self.genSpinC[S] if self.abs_gr[g]==gr-1])
                    if (gr-1) in gradings_in_S:
                        mc_homology_in_S[gr]-=ranks_in_S[gr]
                    if (gr+1) in gradings_in_S:
                        mc_homology_in_S[gr]-=ranks_in_S[gr+1]
                mc_homology.append(dict([(gr,mc_homology_in_S[gr]) for gr in mc_homology_in_S if mc_homology_in_S[gr]!=0]))

        return mc_homology

    def print_action(self):
        if self.there_is_action:
            for g in self.generators:
                print repr(g)+" under the action maps to "+repr(self.image_of_generators[g])
        else:
            print "No action specified."
        return None

    def print_differentials(self):
        if not self.fully_initialized:
            self.generate_domains()
        for S in self.SpinC_structures:
            print "In SpinC structure "+repr(S)
            for i in self.genSpinC[S]:
                print repr(i)+" in grading "+repr(self.abs_gr[i])+" under the differential (possibly) maps to "+repr([j for j in self.genSpinC[S] if self.can_contribute(i,j)])+" with domain types (Euler char, num of boundaries): "+repr([self.domain_type(i,j) for j in self.genSpinC[S] if self.can_contribute(i,j)])
        return None


    def pretty_print(self,output_file=None,interchanges=[],cancellations=[],relevant_SpinC='all'):
        """
        Tries to print the chain complex in human-readable
        format. Only prints the complex whose SpinC structure is in
        the list relevant_SpinC (default is 'all'). User tweaking (via
        interchanges, with a list of tuples of generators to be
        swapped) recommended; the interchanges will be done in the
        order provided. Cancellations is a list of tuples of
        generators where the differential coefficient is known to be
        1, and it performs the cancellations in the order
        provided. 

        Tries to print Z/2-equivariantly if Z/2-action given, and
        relevant_SpinC is fixed under the Z/2-action. So tweak and
        cancel equivariantly in that case if you want the Z/2 symmetry
        to remain. If output_file=None, shows the graphics, else
        stores it.
        """
        
        print_equivariantly=self.there_is_action
        
        if not self.fully_initialized:
            self.generate_domains()
        if relevant_SpinC=='all':
            relevant_SpinC=self.SpinC
            
        if self.there_is_action:
            for S in relevant_SpinC:
                if self.image_of_SpinC_structures[S] not in relevant_SpinC:
                    print_equivariantly=False

        replacements=len(self.generators)*[None,]#the string to be printed instead of generator; empty string kills the generator
        for g in self.generators:
            if self.SpinC[g] in relevant_SpinC:
                replacements[g]=repr(g)
            else:
                replacements[g]=''
        arrows=dict()#the differential arrows in the chain complex
        for g in self.generators:
            for h in self.generators:
                if self.can_contribute(g,h) and self.SpinC[g] in relevant_SpinC and self.SpinC[h] in relevant_SpinC:
                    if self.does_contribute(g,h):
                        arrows[(g,h)]='red'#Coefficient 1
                    else:
                        arrows[(g,h)]='blue'#Coefficient unknown

        #Now the cancellations
        for (a,b) in cancellations:
            if self.abs_gr[a]!=self.abs_gr[b]+1:
                raise Exception("Cancelling arrow not in correct grading")
            if arrows[(a,b)]!='red':
                print "WARNING: Trying to cancel an arrow where the program doesn't know if the coefficient is 1"
                raise Exception("Not clear if the cancelling arrow is coefficient 1")#Comment or uncomment this line, depending on how we feel.

            maps_to_a=[foo for foo in self.generators if (foo,a) in arrows]
            maps_to_b=[foo for foo in self.generators if ((foo,b) in arrows and foo!=a)]
            maps_from_a=[foo for foo in self.generators if ((a,foo) in arrows and foo!=b)]
            maps_from_b=[foo for foo in self.generators if (b,foo) in arrows]

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
                        arrows[(g,h)]='green'#Coefficient still unknown, comes from change of basis.

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

        width=max([len(self.grading_spread[gr]) for gr in self.grading_spread])+1
        x_coordinates=len(self.generators)*[None,]
        #Initial assignment
        for gr in self.grading_spread:
            ll=[g for g in self.grading_spread[gr] if self.SpinC[g] in relevant_SpinC]
            if print_equivariantly:
                fixed=[g for g in ll if self.image_of_generators[g]==g]
                non_fixed=[g for g in ll if  self.image_of_generators[g]!=g]
                distance=width/(len(non_fixed)+2)

                for g in fixed:
                    x_coordinates[g]=0#there could be a collision here.
                non_fixed_reps=[]
                for g in non_fixed:
                    if self.image_of_generators[g] not in non_fixed_reps:
                        non_fixed_reps.append(g)
                for ind_g,g in enumerate(non_fixed_reps):
                    x_coordinates[g]=(ind_g+1)*distance
                    x_coordinates[self.image_of_generators[g]]=-(ind_g+1)*distance

            else:
                distance=width/(len(ll)+1)

                for ind_g,g in enumerate(ll):
                    x_coordinates[g]=-width/2+(ind_g+1)*distance

        #Now the interchanges
        for (a,b) in interchanges:
            if self.abs_gr[a]!=self.abs_gr[b]:
                raise Exception("Trying to interchange stuff in different gradings")
            temp=x_coordinates[a]
            x_coordinates[a]=x_coordinates[b]
            x_coordinates[b]=temp


        G=Graphics()
        for gr in self.grading_spread:
            G+=line([(-width/2,gr),(width/2,gr)],alpha=0.1)#alpha is transparency here. nothing to do with us.
            G+=text(repr(gr)+'+C',(width/2+1,gr),color="blue")


        for g in self.generators:
            if replacements[g]!='':
                G+=text(replacements[g],(x_coordinates[g],self.abs_gr[g]),color="black")
                maps_from_g=[h for h in self.generators if (g,h) in arrows]
                for h in maps_from_g:
                    G+=line([(x_coordinates[g],self.abs_gr[g]),(x_coordinates[h],self.abs_gr[h])],alpha=0.8,color=arrows[(g,h)])

        G.axes(False)
        if output_file==None:
            G.show()
        else:
            G.save(output_file)



def branched_double(H,num_pointed_regions):
    """
    Inputs a HeegaardDiagram H with 2 basepoints (always placed on the
    last two regions). Returns a Heegaard diagram branched along those
    two points. If num_pointed_regions=1, places a single basepoint on
    the preimage of the last region; if num_pointed_regions=2, places
    two basepoints on the preimages of the last two regions. There is
    a natural Z/2 action on the branched double cover, which is also
    returned. (If H came with a Z/2 action on its own, that is lost.)
    """

    if len(H.regions)-len(H.regions_un)!=2:
        raise Exception("Need two pointed regions")
    if num_pointed_regions not in [1,2]:
        raise Exception("Can only return 1 or 2 basepoints")

    #First we find a path from connecting the two basepoints in the complement of beta circles
    shortest_path=H.region_graph_alpha.shortest_path(H.regions[-2],H.regions[-1])
    cut_edge=[next(e[2] for e in H.region_graph_alpha.edges(labels=True) if sorted(list(e[:2]))==sorted(shortest_path[steps:steps+2])) for steps in range(len(shortest_path)-1)]#a cut_edge between the two marked regions, passing only through alpha circles. (The cut_edge is a list of (alpha_circle,arc_on_that_alpha_circle).)  We want the double branch cover of the complement of the cut_edge to be a trivial cover. So we need to ensure that the cut edge intersects each alpha circle an even number of times.

    cut_edge_vector=vector(GF(2),len(H.alphas),[len([e for e in cut_edge if e[0]==a]) for a in H.alphas])
    incidence_matrix=matrix(GF(2),len(H.alphas),len(H.betas),H.intersection_matrix)
    try:
        extra_betas=incidence_matrix.solve_right(cut_edge_vector)#Need to add these beta circles to cut edge.
    except ValueError:
        raise Exception("Need the knot to be null-homologous mod 2, for the program to compute double branched cover.")

    for b in H.betas:
        if (extra_betas[b])%2==1:#We need to add a parallel copy of this beta circle to cut edge
            for p in H.intersections_on_betas[b]:
                [a,n]=H.intersection_incidence[p][0::2]
                if n==1:
                    cut_edge.append((a,p))
                else:
                    ind_p=(H.intersections_on_alphas[a]).index(p)
                    cut_edge.append((a,H.intersections_on_alphas[a][ind_p-1]))

    #Now we are ready to construct the new Heegaard diagram. The complement of the cut_edge has 2 lifts: the 0 and the 1 lift. So most objects also have 2 lifts. 
    new_intersections=[(p,0) for p in H.intersections]+[(p,1) for p in H.intersections]#(p,i) lies in the i lift.
    new_regions=[(R,0) for R in H.regions_un]+[(R,1) for R in H.regions_un]+[(R,-1) for R in H.regions[-2:]]#(R,i) has the first intersection point in the i lift; when i=-1, it is a lift of a pointed region.

    new_image_of_intersections=[p+len(H.intersections) for p in H.intersections]+[p for p in H.intersections]#The Z-2 action is easy to write down.

    new_boundary_intersections=[]
    for (R,i) in new_regions:
        new_boundary_intersections.append([])
        curr_sheet=abs(i)
        for ind_p in range(0,len(H.boundary_intersections[R]),2):
            [curr_p,next_p]=H.boundary_intersections[R][ind_p:ind_p+2]
            curr_alpha=H.intersection_incidence[curr_p][0]
            (new_boundary_intersections[-1]).append(new_intersections.index((curr_p,curr_sheet)))

            #Now check if we move sheets as we go from curr_p to next_p along this alpha arc
            curr_p_ind=(H.intersections_on_alphas[curr_alpha]).index(curr_p)
            if H.intersections_on_alphas[curr_alpha][curr_p_ind-1]==next_p:
                first_p=next_p#Accordingly to orientation of the alpha arc, next_p appears before curr_p, so is first_p.
            else:
                first_p=curr_p
            for foo in cut_edge:
                if (curr_alpha,first_p) == foo:
                    curr_sheet=1-curr_sheet

            (new_boundary_intersections[-1]).append(new_intersections.index((next_p,curr_sheet)))
            
        if i==-1:
            new_boundary_intersections[-1]+=[new_image_of_intersections[p] for p in new_boundary_intersections[-1]]

    return HeegaardDiagram(new_boundary_intersections,num_pointed_regions,new_image_of_intersections)

