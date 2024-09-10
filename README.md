In this repository are Matlab and python codes. Matlab codes are for the transcranial photoacoustic simulations and the python codes for generating a digital phantom of a human head. 

How Matlab codes work

1) The simulation of the acoustic forward problem can be done with both codes for the forward simulation: forward_simulation_just_vessel_p0.m and forward_simulation_whole_p0.m. 
The difference is the used initial pressure distribution. In both simulations the inputs are: i) initial pressure distribution, ii) medium speed of sound map, iii) medium density map, iv) medium attenuation coefficient map.


2) For the reconstruction of the time series data from 1) are different assumptions for the speed of sound in the time reversal method available: i) constant speed of sound (reconstruction_constant_sos.m), ii) a layer for the skull and the rest constant (reconstruction_layer_sos.m), iii) ground truth speed of sound map (reconstruction_perfect_sos.m). 

For all the reconstructions the input are the time series data from 1).
 
ii) needs also the speed of sound map with a layer for the skull as additional input. 

iii) needs the ground truth speed of sound map, already used in 1) as additional input. 

How python codes work

- The python codes in this repository include: i) a vesselness filter developed by Frangi et al. (https://doi.org/10.1007/BFb0056195), ii)perlin noise generation and iii) a vessel tree generator developed by Iyer et al. (https://doi.org/10.1038/s41598-023-44633-2) with an voxalizer.
- i) needs an image as input, in this work it is an MRI
- ii) needs as an input a digital phantom in which the grey and white matter should grow and the repository (https://github.com/Maharshi-Pandya/Perlin-Noise-Implementation.git) needs to be copied
- iii) for using this code the repository (https://github.com/kritiyer/vessel_tree_generator.git) needs to be copied and their dependencies fulfilled
- Also are needed for the codes: matplotlib.pyplot, numpy, random, time, nrrd, scipy.ndimage, os, skimage.filters.

