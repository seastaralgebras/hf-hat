intersections=[[0,1,0,1,0,1,0,1,0,1],[1,0,1,0,1,0,1,0,1]]
#Number the alpha and beta circles 0 to g-1. Number all the intersection points starting at 0, so that the points on alpha_i are before points on alpha_j for all i<j, and the intersections on each alpha-circle are numbered cyclically (orientation and starting point immaterial); the j-th point on alpha_i lies in beta_{intersections[i][j]}.

boundaries=[
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
    ]
#Number the regions (complement of alpha and beta), currently except the basepoint, starting at 0; boundaries[i] is a list of intersection points appearing on the boundary of the i-th region, cyclically ordered as per the orientation of the boundary. Furthermore, the arc joining the first two points in the list lies on some alpha-circle (we will partially check for this in the main program). Unfortunately, this means, from now on, we are assuming all regions are planar with one boundary (which can intersect itself).

abs_gradings_manual=False
#(chosen_generators,their_abs_gradings)=
#Absolute gradings. If you want to add abs_gradings manually, set abs_gradings_manual to True, uncomment the above line and add it after the equality sign (give a list of chosen generators (each generator is an ordered tuple of intersection points), one in each SpinC structure, and their absolute gradings). If not, set abs_gradings_manual to False and comment the above line. 

there_is_action=True
image_of_intersections=[0,9,8,7,6,5,4,3,2,1,18,17,16,15,14,13,12,11,10]
#Sometimes there is an action (of Z/2). If so, set there_is_action to True, uncomment the above line and give the images of the intersection points under this action as a list after the equality sign. If not, set there_is_action to False and comment the above line. 


#Main program
load hf-hat.sage


#Printing
#pretty_print()
#Use pretty_print function. You may input SpinC_structures as a list of SpinC structures to concentrate on (default is everything). If output_file specified, stores it; otherwise, displays it. You may provide interchanges as a list of interchanges of generators to make printing prettier. You may provide cancellations as a list of cancellations of generators to simplify the complex. (If there_is_action, try to do interchanges and cancellations equivariantly.)

#prettier version of the whole complex
pretty_print(output_file='hf-hat-output-for-Sigma_2_3_7_other-0.pdf',interchanges=[(19,10),(0,4),(16,13),(25,44),(36,33),(27,42),(15,10),(19,14),(24,38),(5,31),(17,0),(12,4),(25,38),(31,44),(38,13),(16,31),(27,40),(42,29),(3,41),(1,28),(20,3),(9,1),(3,18),(1,11)])
#after the first round of cancellations
pretty_print(output_file='hf-hat-output-for-Sigma_2_3_7_other-1.pdf',interchanges=[(19,10),(0,4),(16,13),(25,44),(36,33),(27,42),(15,10),(19,14),(24,38),(5,31),(17,0),(12,4),(25,38),(31,44),(38,13),(16,31),(27,40),(42,29),(3,41),(1,28),(20,3),(9,1),(3,18),(1,11)],cancellations=[(23,39),(6,30),(10,33),(19,36),(21,41),(8,28),(37,20),(32,9),(40,18),(29,11),(43,22),(26,7),(27,3),(42,1)])
#after the second round of cancellations
pretty_print(output_file='hf-hat-output-for-Sigma_2_3_7_other-2.pdf',interchanges=[(19,10),(0,4),(16,13),(25,44),(36,33),(27,42),(15,10),(19,14),(24,38),(5,31),(17,0),(12,4),(25,38),(31,44),(38,13),(16,31),(27,40),(42,29),(3,41),(1,28),(20,3),(9,1),(3,18),(1,11)],cancellations=[(23,39),(6,30),(10,33),(19,36),(21,41),(8,28),(37,20),(32,9),(40,18),(29,11),(43,22),(26,7),(27,3),(42,1),(24,35),(5,34),(38,0),(31,4)])

