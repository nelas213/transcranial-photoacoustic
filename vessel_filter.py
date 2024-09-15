

import numpy as np
import nrrd 
import matplotlib.pyplot as plt
from skimage.filters import frangi
from scipy.ndimage import gaussian_filter

input_file=input_file
output_file=output_file

#%%
head = nrrd.read("input_file/7T_MRT.nrrd")[0].astype(np.float64) #load MRI

#vessel_filtered = frangi(head)

#%%
head/= np.max(head)

def smooth(array):
    array=np.where(array>=.5,.5,array)
    return array
def plot(array,axis=0,label=None):
    if axis==0:
        for i in range(array.shape[axis]):
            plt.imshow(array[i,:,:])
            plt.title(label)
            plt.show()
    if axis==1:
        for i in range(array.shape[axis]):
            plt.imshow(array[:,i,:])
            plt.title(label)
            plt.show()
    if axis==2:
        for i in range(array.shape[axis]):
            plt.imshow(array[:,:,i])
            plt.title(label)
            plt.show()
    
    
newhead = smooth(head)

blurred = gaussian_filter(newhead, (2,1,1))
blurred/= np.max(blurred)

newhead = newhead*(1-blurred)
#newhead = newhead*np.where(newhead > 0, 1, 0)

#%%
def threshold(array,t):
    return np.where(array<t,0,array)



threshead = threshold(head , .050)[:150,:150,110:]

#%%
alphs,bets = [.1,.05,.05,1,.01,.001],[10,15,20,5,20,30]

for alph,bet in zip(alphs,bets):
    vessel_filt = frangi(threshead,alpha = alph,beta = bet)
    plot(vessel_filt,axis=1,label=f"α={alph},β={bet}")

#%%
vessel_filt = frangi(threshead,alpha = .1, beta = 40)
vessel_filt = np.where(vessel_filt>.01,1,0)
plt.imshow(np.max(vessel_filt,axis=1))



#%%
#newnewhead = newhead[:150,:150,110:]
nrrd.write("output_file/vessels.nrrd",vessel_filt)
