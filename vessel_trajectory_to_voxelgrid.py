
import numpy as np
import os
import nrrd
import matplotlib.pyplot as plt

input_file=input_file
output_file=output_file


path = 'path'#path to generated vessel from vessel_tree_generator

TEST = [np.load(path+element) for element in os.listdir(path)]

testtree=TEST[0]
tree_two = TEST[1]
#%%


head,_ = nrrd.read("input_file/larger_mini_csf+skull.nrrd")
head = head[:,:,::-1]
head = np.where(head==1,2,head)

I = head.shape[2]
while 2 in head[:,:,I-1]:
    I-=1
    
    
I=140

#%%




def scale_trajectory(tree):
    x_min,y_min,z_min = np.min(tree[:,:,0]),np.min(tree[:,:,1]),np.min(tree[:,:,2])
    tree[:,:,0]-=x_min
    tree[:,:,1]-=y_min
    tree[:,:,2]-=z_min
    return tree

testtree = scale_trajectory(testtree)
testree = testtree[:,:,::-1]


tree_two = scale_trajectory(tree_two)
def mark_voxels(grid, point, radius, voxel_spacing, zvoxel_spacing):
    size = grid.shape[0]
    zsize = grid.shape[2]
    x, y, z = point
    xi, yi, zi = int(x / voxel_spacing), int(y / voxel_spacing), int(z / zvoxel_spacing)
    radius_in_voxels = int(radius / voxel_spacing)
    
    for i in range(max(0, xi-radius_in_voxels), min(size, xi+radius_in_voxels+1)):
        for j in range(max(0, yi-radius_in_voxels), min(size, yi+radius_in_voxels+1)):
            for k in range(max(0, zi-radius_in_voxels), min(zsize, zi+radius_in_voxels+1)):
                if np.sqrt((xi-i)**2 + (yi-j)**2 + (zi-k)**2) <= radius_in_voxels:
                    grid[i, j, k] = 1

def trajectory_to_voxels(grid, trajectory, radii, voxel_spacing, zvoxel_spacing):
    for point, radius in zip(trajectory, radii):
        mark_voxels(grid, point, radius, voxel_spacing,zvoxel_spacing)
    return grid

size = head.shape[0]  # here the gridsize
zsize = head.shape[2]-I
voxel_spacing = np.max(testtree[:,:,:3])*1e3/size
zvoxel_spacing = np.max(testtree[:,:,:3])*1e3/zsize

# voxel_spacing = .125 # here the spacing

grid = np.zeros((size, size, zsize), dtype=np.int32)
trajectory = testtree[0][:,:3]*1e3
radii = testtree[0][:,3]*2e3
grid = trajectory_to_voxels(grid, trajectory, radii, voxel_spacing, zvoxel_spacing)

trajectory = tree_two[0][:,:3]*1e3
radii = tree_two[0][:,3]*2e3
grid = trajectory_to_voxels(grid, trajectory, radii, voxel_spacing, zvoxel_spacing)


def do_branches(grid,tree):
    branches = tree[1:]
    for branch in branches:
        trajectory = branch[:,:3]*1e3
        radii = branch[:,3]*1.5e3
        grid = trajectory_to_voxels(grid, trajectory, radii, voxel_spacing, zvoxel_spacing)

do_branches(grid,testtree)
do_branches(grid,tree_two)



PIC = np.concatenate((np.zeros((head.shape[0],head.shape[1],head.shape[2]-grid.shape[2])),grid),axis=2)



final_head = np.where(PIC==1,PIC,head)

for i in range(head.shape[0]):
    #plt.figure(figsize=(15,15))
    plt.imshow(final_head[:,i,:]), plt.show()
#%%
nrrd.write("output_file/larger_mini_with_vessel.nrrd",final_head)