######
import numpy as np
import matplotlib.pyplot as plt
import time
import seaborn as sns

# Initialise ground 
# Ground is initialised here by fisrst making a 256-long column and repeating it 255 more times vertically
a = np.arange(0,256).reshape((256,1))[::-1]      
ground = np.concatenate([ a for i in range(256) ],axis=1)     # Shape 256x256. Each cell has a number from 0-255 which denotes r of rgb. Later on, repeat it two more times to get (r,g,b).
del a

# The np array, ground, is now made and is left untampered for the rest of the code.

# Now we make box1 and box2 - which stand for alleles 1 and 2 of all loci of qll organisms. So their dimensions must be 256x256x255 each
# boxes = box1 + box2 holds the phenotype. Since CC = 1+1 = 2, it is converted to 1.
st = time.time()
numm=256
box1 = np.concatenate([ np.concatenate([(np.random.rand(numm-1).reshape(1,1,numm-1)>=0.5)*1 for i in range(numm)],axis=1) for i in range(numm) ], axis=0)
box2 = np.concatenate([ np.concatenate([(np.random.rand(numm-1).reshape(1,1,numm-1)>=0.5)*1 for i in range(numm)],axis=1) for i in range(numm) ], axis=0)

boxes = box1+box2
# Summing boxes along the 3d-row (axis = 2) gives a 256x256 array with numbers from 0-255. Each number represents organism's colour (r of rgb)
boxes_repl = (np.where(boxes==2, 1, boxes)).sum(axis=2)

print((time.time()-st))   # Get the time taken for this subroutine
del st, numm


# calcprob takes in the allele cuboids and returns a 256x256 array with 1s and 0s where 1 indicates alive, 0 means dead.

def calcprob(box1, box2): # Takes in two 256x256x255 matrices and returns two 256x256x255 matrices
    boxes = box1 + box2
    boxes_repl = (np.where(boxes==2, 1, boxes)).sum(axis=2)*3
    org_diff = np.absolute(boxes_repl - ground)   # Absolute value needed since it's the magnitude of diff that matters not the sign
    
    # y = (0.85/765)*x + 0.05 is the equation that yields probability of death (min = 0.05, max = 0.9) with the input being the difference in ground colour (rgb sum) and the organism colour (rgb sum)
    org_diff = (0.85/765)*org_diff 
    org_diff = org_diff + 0.05
    
    # Now, org_diff is a 256x256 matrix with probabilities of death. It needs to be converted to 1s and 0s
    org_diff = (np.random.rand(256,256) <= org_diff)*1  # After this step, live = 0, dead = 1. WE need the reverse. 
    # so we do:
    org_diff = 1-org_diff # Thus, now, all live = 1, all dead = 0
    return org_diff

# Declare a 256x256 ones matrix from which elements become zero as and when they are done reproducing.
donecells = np.ones((256,256))

# Now define functions for updating boxes during one round of replication
def update_onetimeinonegen(boxa,boxb,donecells):
    i, j = np.nonzero(donecells)
    coord = np.random.choice(len(i),2)
    coord = (i[coord[0]], j[coord[1]]) # Get a random coord which corresponds to 1 (or undone yet element) in donecells
    #if (org_diff[coord[0], coord[1]] == np.zeros((255,))).all() == True: # Check if it is an empty space in the ground
    org_diff = calcprob(boxa, boxb)
    #orga = np.repeat(org_diff[:,:,np.newaxis], 255, axis=2)
    for acount in range(256):
        for bcount in range(256):
            if org_diff[acount,bcount] == 0:
                box1[acount,bcount,:] = np.zeros((255,))
                box2[acount,bcount,:] = np.zeros((255,))
    #box1 = np.where(orga==0, 0, boxa)
    #box2 = np.where(orga==0, 0, boxb)
    if (org_diff[coord[0], coord[1]] == 0): # Check if it is an empty space in the ground
        neighbours_of_the_dead = get_neighbours(coord)  
        neighbours_of_the_dead = [i for i in neighbours_of_the_dead if org_diff[i[0],i[1]] == 1 ] # Select only those neighbours which are not zero themselves
        if neighbours_of_the_dead != []:
            #print(neighbours_of_the_dead)
            neighbour = neighbours_of_the_dead[np.random.choice(len(neighbours_of_the_dead))] # Now neighbour chosen, we need to update coord with whatever is there in the neighbour spot
            box1[coord[0],coord[1],:] = box1[neighbour[0], neighbour[1], :] # Whatever is there in box1 in the neighbour's spot, put it our coord's spot
            box2[coord[0],coord[1],:] = box2[neighbour[0], neighbour[1], :] # Repeat the same for box2
            org_diff[coord[0], coord[1]] = 1
        # If neighbours_of_the dead is empty that is no live neighbours, ignore the coord cell and let it remain dead for the generation.
    else: # So in this case, org_diff[coord[0], coord[1]] has to be 1 since it is not 0
        neighbours_of_the_alive = get_neighbours(coord)
        neighbours_of_the_alive = [i for i in neighbours_of_the_alive if org_diff[i[0],i[1]] == 1 ] # Select only those neighbours which are not zero themselves
        if neighbours_of_the_alive != []:
            neighbour = neighbours_of_the_alive[np.random.choice(len(neighbours_of_the_alive))]
            progeny = gen_progeny(box1,box2,coord,neighbour) # List of two progenies basically two 255x2 arrays
            box1[coord[0],coord[1],:] = progeny[0][:,0] # update all first alleles of coord position with progeny 1
            box2[coord[0],coord[1],:] = progeny[0][:,1] # update all second alleles of coord position with progeny 1
            box1[neighbour[0],neighbour[1],:] = progeny[1][:,0] # update all first alleles of coord position with progeny 1
            box2[neighbour[0],neighbour[1],:] = progeny[1][:,1] # update all second alleles of coord position with progeny 1
        # If neighbours_of_the dead is empty that is no live neighbours, ignore the coord cell and let it remain dead for the generation.
    donecells[coord[0],coord[1]] = 0
    return box1, box2, donecells

def gen_progeny(box1,box2,coord,neighbour):
    a1oc = box1[coord[0],coord[1],:].reshape(255,1)
    a2oc = box2[coord[0],coord[1],:].reshape(255,1)
    oc = np.concatenate((a1oc,a2oc),axis=1)
    a1on = box1[neighbour[0],neighbour[1],:].reshape(255,1)
    a2on = box2[neighbour[0],neighbour[1],:].reshape(255,1)
    on = np.concatenate((a1on,a2on),axis=1)
    gam1 = gen_gamete(oc)
    gam2 = gen_gamete(on)
    progeny1 = np.concatenate((gam1,gam2), axis=1)
    gam3 = gen_gamete(oc)
    gam4 = gen_gamete(on)
    progeny2 = np.concatenate((gam3,gam4), axis=1)
    return [progeny1, progeny2]

def gen_gamete(c): # Needs one 255x2 arrays to work. In each of the 255 rows, it selects which allele goes into the gamete.
    what_to_choose_per_locus = (np.random.rand(255) >= 0.5)*1 # In each locus should I choose the 0th element (allele 1) or the first element (allele 2)?
    gamete = np.array([c[i][what_to_choose_per_locus[i]] for i in range(255)])
    return gamete.reshape((255,1))

def get_neighbours(coord):
    if (not coord[0] in (0,255)) and (not coord[1] in (0,255)): # If the cell is not on the border
        neighbours = [ (coord[0]-1, coord[1]), (coord[0], coord[1]+1), (coord[0]+1, coord[1]), (coord[0], coord[1]-1)]
    elif coord[0] == 0:
        if coord[1] == 0: # topleft corner
            neighbours = [ (0,1), (1,0)]
        elif coord[1] == 255: # topright corner
            neighbours = [ (0,254), (1,255)]
        else: # middle of top border
            neighbours = [ (coord[0]+1, coord[1]), (coord[0], coord[1]-1), (coord[0], coord[1]+1)]
    elif coord[1] == 0:
        if coord[0] == 255: # bottomleft corner
            neighbours = [ (254,0), (255,1)]
        else: # middle of left border
            neighbours = [ (coord[0], coord[1]+1), (coord[0]+1, coord[1]), (coord[0]-1, coord[1])]
    elif coord[1] == 255:
        if coord[0] == 255: # bottomright corner
            neighbours = [ (255,254), (254,255)]
        else: # middle of right border
            neighbours = [ (coord[0]-1, coord[1]), (coord[0], coord[1]-1), (coord[0]+1, coord[1])]
    else: # coord[0] == 255 and coord[1] is not 0 or 255 which means middle of the down border
        neighbours = [ (coord[0]-1, coord[1]), (coord[0], coord[1]+1), (coord[0], coord[1]-1)]
    return neighbours
    
# Start reproduction - ours is gen 0 right now
generations = 3
i = 0
for i in range(generations):
    startpaint=time.time()
    boxes = box1 + box2
    boxes_repl = (np.where(boxes==2, 1, boxes)).sum(axis=2)*3
    
    fig, ax = plt.subplots(figsize=(20,20))
    sns.heatmap(ground,cmap='Greys_r',linewidths=0.1, linecolor='yellow')
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    #for j in range(256):
    #    for k in range(256):
    #        plt.scatter(j+0.5,k+0.5,c=(boxes_repl[j,k]/(255*3), boxes_repl[j,k]/(255*3), boxes_repl[j,k]/(255*3)),s=0.5) #boxes_repl[j,:]
    for j in range(256):
        colors=np.concatenate((boxes_repl[j,:].reshape(256,1)/(255*3), boxes_repl[j,:].reshape(256,1)/(255*3), boxes_repl[j,:].reshape(256,1)/(255*3)), axis=1)
        plt.scatter(np.ones(256)*(j+0.5),np.arange(256)+0.5,c=colors,s=0.5) #boxes_repl[j,:]
    plt.savefig('testfig_'+str(i)+'.pdf')
    plt.close()
    print(time.time()-startpaint)
    del startpaint, j
    
    # Update boxes
    start = time.time()
    countt = 0
    while 1 in donecells:
        box1, box2, donecells = update_onetimeinonegen(box1, box2, donecells)
        countt += 1
        print(countt)
        if countt % 1000 == 0:
            print(time.time()-start)
            print(countt)
            start = time.time()
    print(countt)
    del star, countt
    
    # Now that one gen is over make donecells ones again for next round
    donecells = np.ones((256,256))
    
 
 ######  #    #  #####
 #       ##   #  #    #
 #####   # #  #  #    #
 #       #  # #  #    #
 #       #   ##  #    #
 ######  #    #  #####
      

'''def calcprob(colour_all): # Takes in 256x256 ,matrix of rgb sums and returns 0 wherever death elsewhere 1
    org_diff = np.absolute(colour_all - ground)
    org_diff = (0.85/765)*org_diff
    org_diff = org_diff + 0.05
    # y = (0.85/765)*x + 0.05 is the equation that yields probability of death (min = 0.05, max = 0.9) with the input being the difference in ground colour (rgb sum) and the organism colour (rgb sum)
    org_diff = (np.random.rand(4,4) <= org_diff)*1
    return org_diff'''

'''boxes_repl = calcprob(boxes_repl)

prob_death = calcprob(colour_all)
donecells = np.ones((4,4))
flag=True'''

'''def positiongen(chosen_one):
    return (chosen_one[0]-1)*4 + (chosen_one[1]-1)
    # This function is used to convert say (1,2) into 6 using the following 256x256 matrix 
    # 0 1 2 3
    # 4 5 6 7
    # 8 9 10 11
    # 12 13 14 15'''



            
            '''# Now we need a list of neighbours of the dead cell
            try:  # +
                # This tries thinking four neighbour exist. If four don't exist an error is raised so see except
                neighbours_of_the_dead = [ org_diff[(coord[0]-1, coord[1])], org_diff[(coord[0], coord[1]+1)], org_diff[(coord[0]+1, coord[1])], org_diff[(coord[0], coord[1]-1)]]
            except:
                try: # |-
                    # This tries thinking three neighbours exist. If three don't exist, an error is raised so see except
                    neighbours_of_the_dead = [ org_diff[(coord[0], coord[1]+1)], org_diff[(coord[0]+1, coord[1])], org_diff[(coord[0], coord[1]-1)]]
                except: 
                    try: # T
                        neighbours_of_the_dead = [ org_diff[(coord[0]-1, coord[1])], org_diff[(coord[0], coord[1]-1)], org_diff[(coord[0]+1, coord[1])]]
                    except:
                        try: # -|
                            neighbours_of_the_dead = [ org_diff[(coord[0]-1, coord[1])], org_diff[(coord[0], coord[1]+1)], org_diff[(coord[0], coord[1]-1)]]
                        except:
                            try: # inverted T
                                neighbours_of_the_dead = [ org_diff[(coord[0]-1, coord[1])], org_diff[(coord[0], coord[1]+1)], org_diff[(coord[0]+1, coord[1])]]
                            except: # Now we start two neighbout possibilities
                                try: # topleft corner
                                    neighbours_of_the_dead = [ org_diff[(coord[0], coord[1]-1)], org_diff[(coord[0]+1, coord[1])]]
                                except:
                                    try: # topright corner
                                        neighbours_of_the_dead = [ org_diff[(coord[0]-1, coord[1])], org_diff[(coord[0], coord[1]-1)]]
                                    except:
                                        try: # bottomright corner
                                            neighbours_of_the_dead = [ org_diff[(coord[0]-1, coord[1])], org_diff[(coord[0], coord[1]+1)]]
                                        except:
                                            try: # bottomleft corner
                                                neighbours_of_the_dead = [ org_diff[(coord[0], coord[1]+1)], org_diff[(coord[0]+1, coord[1])]]
                                            except:
                                                pass'''
                                                
'''while flag: # The flag becomes false when donecells becomes a null matrix
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
        neighbours.append(prob_death[chosen_one[0]+1, chosen_one[1]])'''
        
'''dfA = pd.DataFrame([[1,0],[1,1],[0,0]], columns=list("XY"), index=list("ABC"))
    dfB = pd.DataFrame([[0,1],[1,1],[0,1]], columns=list("XY"), index=list("ABC"))
    
    assert dfA.shape == dfB.shape
    
    x = np.arange(0,len(dfA.columns))
    y = np.arange(0,len(dfB.index))
    X,Y=np.meshgrid(x,y)
    
    fig, ax = plt.subplots(figsize=(20,20))
    #ax.invert_yaxis()
    ax.imshow(ground, aspect="auto", cmap="Greys_r")
    
    cond = dfB.values == 1
    ax.scatter(X[cond], Y[cond], c="crimson", s=100)
    
    ax.set_xticks(x)
    ax.set_yticks(y)
    ax.set_xticklabels(dfA.columns)
    ax.set_yticklabels(dfA.index)
    plt.show()'''
    
    '''N = 10
    M = 11
    #ylabels = ["".join(np.random.choice(list("PQRSTUVXYZ"), size=7)) for _ in range(N)]
    #xlabels = ["".join(np.random.choice(list("ABCDE"), size=3)) for _ in range(M)]
    
    x, y = np.meshgrid(, ground[])
    s = np.random.randint(0, 180, size=(N,M))
    c = np.random.rand(N, M)-0.5
    
    fig, ax = plt.subplots()
    
    R = s/s.max()/2
    circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
    col = PatchCollection(circles, array=c.flatten(), cmap="RdYlGn")
    ax.add_collection(col)
    
    ax.set(xticks=np.arange(M), yticks=np.arange(N),
           xticklabels=xlabels, yticklabels=ylabels)
    ax.set_xticks(np.arange(M+1)-0.5, minor=True)
    ax.set_yticks(np.arange(N+1)-0.5, minor=True)
    ax.grid(which='minor')
    
    fig.colorbar(col)
    plt.show()
'''
    
'''org_diff = np.where(org_diff==1,2,org_diff)
    org_diff = np.where(org_diff==1,0,org_diff)
    org_diff = np.where(org_diff==2,1,org_diff)
'''

    
# Now one must change box1 and box2 with regards to dead cells since.
# To compare using np.where we need a 3-D array so we make this by copying org_diff 254 more times.

#def orggen():
#    one_organism = (np.random.rand(255,2)>=0.5)*1 # So this generates 255 pairs of alleles (1,0)
#    one_organism = one_organism.sum(axis=1) # Sum of both alleles - probable outputs are 0, 1, 2
#    one_organism = np.where(one_organism==2, 1, one_organism) # Replace 2's with 1 Since this is perfect dominance scenario so Aa = AA = aA and are indistinguishable
#    colour_organism = one_organism.sum()*3 # Get sum of loci and multiply by 3 # The sum of loci can be anywhere from 0 to 255. Hence color can be from (0,0,0) to (255,255,255)
#    return one_organism # Return all loci, final colour
    #return (one_organism, colour_organism) # Return all loci, final colour

#organisms_all = [ orggen() for i in range(256*256) ]
#organisms_all = [ orggen() for i in range(256*256) ]
#colour_all = np.array([i for i in organisms_all]).reshape((256,256))
