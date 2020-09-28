import kingman
import csv
import random

"""----------------------------------------------------------------
   ----------------------------------------------------------------
   -----         Simulation of one Kingman coalescent         -----
   ----------------------------------------------------------------
   ----------------------------------------------------------------"""

#---- Simulation of a Kingman coalescent: here between 4 individuals.
test = kingman.simulate(4, random_seed=None)

#---- Measure of statistics on the tree: 
#        - time before the MRCA (most recent common ancestor)
#        - number of cherries (coalescence between two individuals that exist at t=0)
#        - length of the external branches (time before the first coalescence of individuals that exist at t=0)
tMRCA = max(test[1])
print("tMRCA=",tMRCA)

cherries = 0
extBranch=0

for indID in range(len(test[0])):
    for ind2ID in range(indID+1, len(test[0])):
        if test[0][indID]==test[0][ind2ID]:
            if (test[1][indID]==0) & (test[1][ind2ID]==0):
                cherries+=1
                extBranch += test[1][test[0][indID]]
                extBranch += test[1][test[0][ind2ID]]
            elif (test[1][indID]==0):
                extBranch += test[1][test[0][indID]]
            elif (test[1][ind2ID]==0):
                extBranch += test[1][test[0][ind2ID]]


print("cherries=",cherries)
print("extBranch",extBranch)



"""----------------------------------------------------------------
   ----------------------------------------------------------------
   -----       Simulation of several Kingman coalescent       -----
   ----------------------------------------------------------------
   ----------------------------------------------------------------"""


#---- open a file to save the statistics
file = open('kingman_stat','a')

for repet in range(5000):
    #---- simulation of a Kingman with 1000 sampled individuals
    test = kingman.simulate(1000, random_seed=None)
    
    #---- Measure of statistics on the tree:tMRCA, Cherries, extBranch 
    tMRCA = max(test[1])
 
    cherries = 0
    extBranch=0
 
    for indID in range(len(test[0])):
        for ind2ID in range(indID+1, len(test[0])):
            if test[0][indID]==test[0][ind2ID]:
                if (test[1][indID]==0) & (test[1][ind2ID]==0):
                    cherries+=1
                    extBranch += test[1][test[0][indID]]
                    extBranch += test[1][test[0][ind2ID]]
                elif (test[1][indID]==0):
                    extBranch += test[1][test[0][indID]]
                elif (test[1][ind2ID]==0):
                    extBranch += test[1][test[0][ind2ID]]
    
    #---- Save the statistics        
    datawriter=csv.writer(file)
    datawriter.writerow([cherries,extBranch,tMRCA])


file.close()



"""----------------------------------------------------------------
   ----------------------------------------------------------------
   -----       Simulation of several Kingman coalescent       -----
   -----         choice of the coalescence coefficient        -----
   ----------------------------------------------------------------
   ----------------------------------------------------------------"""

                    
#---- open a file to save the statistics
file_Coeff = open('kingman_stat_coeff','a')
               
for repet in range(5000):
    
    #---- simulate a Kingman Coalescent with 1000 individuals, with a rate of coalescence = 0.007 * nb of couples
    #    the code is extracted from the package Kingman.
    sample_size=1000
    random.seed()
    time = [0 for j in range(2 * sample_size)]
    parent = [0 for j in range(2 * sample_size)]
    time[0] = -1
    parent[0] = -1
    ancestors = list(range(1, sample_size + 1))
    t = 0
    next_node = sample_size + 1
    for n in range(sample_size, 1, -1):
        t += random.expovariate(0.00765* n * (n - 1))
        for _ in range(2):
            child = random.choice(ancestors)
            parent[child] = next_node
            ancestors.remove(child)
        ancestors.append(next_node)
        time[next_node] = t
        next_node += 1
    
    test = (parent, time)

    
    #---- Measure of statistics on the tree:tMRCA, Cherries, extBranch 
    tMRCA = max(test[1])

    cherries = 0
    extBranch=0

    for indID in range(len(test[0])):
        for ind2ID in range(indID+1, len(test[0])):
            if test[0][indID]==test[0][ind2ID]:
                if (test[1][indID]==0) & (test[1][ind2ID]==0):
                    cherries+=1
                    extBranch += test[1][test[0][indID]]
                    extBranch += test[1][test[0][ind2ID]]
                elif (test[1][indID]==0):
                    extBranch += test[1][test[0][indID]]
                elif (test[1][ind2ID]==0):
                    extBranch += test[1][test[0][ind2ID]]

    datawriter=csv.writer(file_Coeff)
    datawriter.writerow([cherries,extBranch,tMRCA])


file_Coeff.close()