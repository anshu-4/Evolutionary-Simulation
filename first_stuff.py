######
import numpy as np
import matplotlib.pyplot as plt

# Initialise ground 
# Ground is initialised here by fisrst making a 256-long column and repeating it 255 more times vertically
a = np.arange(0,256).reshape((256,1))[::-1]      
ground = np.concatenate([ a for i in range(256) ],axis=1)     # Shape 256x256
del a

# The np array, ground, is now made and is left untampered for the rest of the code.

def orggen():
    one_organism = (np.random.rand(255,2)>=0.5)*1
    one_organism = one_organism.sum(axis=1) # Sum of both alleles - probable outputs are 0, 1, 2
    one_organism = np.where(one_organism==2, 1, one_organism) # Replace 2's with 1
    colour_organism = one_organism.sum()*3 # Get sum of loci and multiply by 3
    return (one_organism, colour_organism) # Return all loci, final colour

organisms_all = [ orggen() for i in range(16) ]
colour_all = np.array([i[1] for i in organisms_all]).reshape((4,4))

def calcprob(colour_all): # Takes in 256x256 ,matrix of rgb sums and returns 0 wherever death elsewhere 1
    org_diff = np.absolute(colour_all - ground)
    org_diff = (0.85/765)*org_diff
    org_diff = org_diff + 0.05
    # y = (0.85/765)*x + 0.05 is the equation that yields probability of death (min = 0.05, max = 0.9) with the input being the difference in ground colour (rgb sum) and the organism colour (rgb sum)
    org_diff = (np.random.rand(4,4) <= org_diff)*1
    return org_diff

prob_death = calcprob(colour_all)
donecells = np.ones((4,4))
flag=True

def positiongen(chosen_one):
    return (chosen_one[0]-1)*4 + (chosen_one[1]-1)
    # This function is used to convert say (1,2) into 6 using the following 256x256 matrix 
    # 0 1 2 3
    # 4 5 6 7
    # 8 9 10 11
    # 12 13 14 15

while flag: # The flag becomes false when donecells becomes a null matrix
    if donecells == np.zeros((4,4)):
        flag==False
        break
    chosen_one = (np.random.randint(0,4), np.random.randint(0,4))
    while prob_death(chosen_one) == 0 and donecells(chosen_one) == 0: # If that cell is empty then change the chosen one
        chosen_one = (np.random.randint(0,4), np.random.randint(0,4))
    # Now the chosen cell in the ground is non-empty for sure
    position = positiongen(chosen_one)
    neighbour = gen_randneighbour(position)
    if prob_cell(neighbour) == 0:
        colour_all(neighbour) = colour_all(chosen_one)
    else:
        colour_all()
        
    if (chosen_one[1]-1)<0:
        neighbours.append(prob_death[chosen_one[0], chosen_one[1]-1])
    if:
        neighbours.append(prob_death[chosen_one[0], chosen_one[1]+1])
    if (chosen_one[0]-1)<0:
        neighbours.append(prob_death[chosen_one[0]-1, chosen_one[1]])
    if:
        neighbours.append(prob_death[chosen_one[0]+1, chosen_one[1]])
    
