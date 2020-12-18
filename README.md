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
![Locally Aligned](/localalign_example/IGOTA_gauss.png)

## segmentfibers.m
This script is the function for segmenting fibers. The method for segmentation is based on two published papers:

**We implement a fiber tracing method developed by**

> Rogge, H., Artelt, N., Endlich, N. & Endlich, K. Automated segmentation
> and quantification of actin stress fibres undergoing experimentally
> induced changes.* J. Microsc. (2017).

**We reconstruct the fiber pieces using a pairing algorithm developed by**

> Zhang, Z., Xia, S. & Kanchanawong, P. An integrated enhancement and
> reconstruction strategy for the quantitative extraction of actin stress
> fibers from fluorescence micrographs. BMC Bioinformatics. (2017).

## SegmentFibers_validation.m
This script takes in an image and runs the segmentation algorithm.
Here is an example of two segmentation results:
![Segmentation example](SegmentationExample.png)