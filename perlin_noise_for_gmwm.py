

from perlin_noise import PerlinNoise
import matplotlib.pyplot as plt
import numpy as np
import random
import time
from scipy.ndimage import gaussian_filter, binary_dilation
import nrrd
"""
how to handle:
PIC creates a 3d-perlin noise map. at 155x155x155 takes about 130s. 
to increase details, increase octaves.

mask cuts, where there isnt supposed to be any matter, slope is the "strength" of the cut
"""

start = time.time()

noise = PerlinNoise(octaves=2)


# head,_ = nrrd.read("C:/Users/carlo/OneDrive/Kopf Segmentiererei/new segmentation/neues Phantom/mini_phantom.nrrd")
# head = head[:,:,::-1]


head,_ = nrrd.read("C:/Users/tpas/Desktop/Dateidatei/Simple Segmentation/vessel tree/larger_mini_with_vessel.nrrd")


shape = head.shape[0]
xpix,ypix,zpix = head.shape

I = 0
while 2 not in head[:,:,I]:
    I+=1

#I = 200
PIC = np.array([[[noise([j/xpix, i/ypix, k/zpix]) for j in range(zpix-I+5)] for i in range(ypix)] for k in range(xpix)])+1

PIC = np.concatenate((np.zeros((shape,shape,head.shape[2]-PIC.shape[2])),PIC),axis=2)

#%%




# better mask, uses the phantom.
mask = ~np.where(head==2,1,0).astype(bool)
for _ in range(10):
    mask = binary_dilation(mask)
mask =(~mask).astype(float)
mask = np.where(mask==0,0,gaussian_filter(mask,sigma=20,truncate=100)**2)





def create_mask(slopelength:int=10,position:int=0):
    mask = np.concatenate((np.zeros((shape,shape,shape//2-slopelength//2+position)),np.ones((shape,shape,slopelength)) * np.linspace(0,1,slopelength)),axis=2)
    return np.concatenate((mask, np.ones((shape,shape,shape-mask.shape[2]))),axis=2)

def digitalize(pic):
    pic = np.where(pic>1.,pic,0)
    pic = np.where(pic>1.2,2,pic)
    return pic.astype(dtype=np.int16)


#mask = create_mask(80,30)

pic =PIC*mask
pic = digitalize(pic)

end = time.time()

print(end-start, "s")

#%%
Phantom = np.where(pic==1,6,head)
Phantom = np.where(pic==2,7,Phantom)

nrrd.write("C:/Users/tpas/Desktop/Dateidatei/Simple Segmentation/vessel tree/larger_mini_phantom_grown.nrrd",Phantom)
#%%
for i in range(shape):
    #plt.figure(figsize=(15,15))
    plt.imshow(Phantom[i,:,:]), plt.show()
