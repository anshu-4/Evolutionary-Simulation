#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import random
import numpy as np
nrows=50
ncols=50
fig, ax = plt.subplots(1)
bg=np.zeros((nrows,ncols,3))#background - 3d matrix. x-y values+rbg number
o=np.zeros((nrows,ncols,1))#organism - 3d matrix. x-y +rbg
for i in range(nrows):
    for j in range(ncols):
        o[i][j][0]=random.randrange(0,255)#choosing the phenotype of gen_0 randomly
        c = plt.Circle((i+0.5, j+0.5), 0.25, color=(o[i][j][0]/255.0,o[i][j][0]/255.0,o[i][j][0]/255.0))
        ax.add_artist(c)
        for k in range(3):
           bg[i][j][k]=i/255.0 #adding colors to the ground 
plt.imshow(bg,aspect='equal',extent=(0,ncols,0,nrows))
plt.xlim(0,nrows)
plt.ylim(0,nrows)
#plt.xticks(range(ncols+1))
#plt.yticks(range(ncols+1))
ax.set_xticks(range(ncols+1),minor='True')
ax.set_yticks(range(nrows+1),minor='True')

#plt.grid(linestyle='-',axis="both")
plt.show()


# # DEATHS

# In[ ]:




