# SegmentFibers
Stress Fiber segmentation algorithm

## Generate_FiberImage.m
This script generates a synthetic fiber image composed of random fibers with customizable 
* image size `x` and `y`
* number of fibers `nfib`
* fiber widths (mean and std) 
* fiber lengths (mean and std)
* `localalignment_flag` can be set to introduce pockets of local alignment

The script outputs two image files and a *.mat file all with a random 5 letter name tag
* **AAAAA_binary.tif** is a binary image of the fibers
* **AAAAA_gauss.tif** is a gray scale, blurred image of the fibers
* **AAAAA_fiberinfo.mat** is a struct variable with the fiber information

Here is an example of a randomly aligned fiber image:
![Randomly Aligned](/randomalign_example/PYWDF_gauss.png)

Here is an example of a locally aligned fiber image
![Randomly Aligned](/localalign_example/IGOTA_gauss.png)
