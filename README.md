# SegmentFibers
Generate an image with fibers of known width and length distributions for testing fiber segmentation algorithms. Segment fibers based on methods developed by [Rogge+ 2017](https://doi.org/10.1111/jmi.12593) and [Zhang+ 2017](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1684-y).

## Testing the algorithm

### Generate_FiberImage.m
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


### segmentfibers.m
This script is the function for segmenting fibers. The method for segmentation is based on two published papers:

**We implement a fiber tracing method developed by**

> Rogge, H., Artelt, N., Endlich, N. & Endlich, K. Automated segmentation
> and quantification of actin stress fibres undergoing experimentally
> induced changes.* J. Microsc. (2017).

**We reconstruct the fiber pieces using a pairing algorithm developed by**

> Zhang, Z., Xia, S. & Kanchanawong, P. An integrated enhancement and
> reconstruction strategy for the quantitative extraction of actin stress
> fibers from fluorescence micrographs. BMC Bioinformatics. (2017).


### SegmentFibers_validation.m
This script takes in an image and runs the segmentation algorithm.
Here is an example of two segmentation results:
![Segmentation example](SegmentationExample.png)


## Applying segmentfibers to sample data
The segmentfibers algorithm was used to segment basal cell stress fibers from human bronchial epithelial cells.

### test_images

Three sets of data are included in the folder `test_images`. Each set includes a *.tif image file of the basal stress fibers and a *.mat file of the [previously segmented basal cell boundaries](https://github.com/tkphung/CellSegmentation).

The *.tif files are
* Maximum intensity projection from F-actin z-stack images of pseudostratified epithelium
* Naming convention: Donor_Timepoint(hr)_Control/Pressure Treatment & Well #_Field of View_Region.tif
* Example: U13_72_P1_5_SF.tif  
    * Donor: U13
    * Timepoint: 72 hour after treatment
    * Treatment: Pressure, Well #1
    * Field of View ID: #5
    * Region: Stress Fiber (SF)

The corresponding *.mat files are logical matrices featuring the basal cell body segmentation from the same sample.

### Shell_segment_SF.m

The script `Shell_segment_SF.m` includes blocks of code using the segmentation functions and visualizing results.

## Authors

Thien-Khoi N. Phung
[@tkphung](https://twitter.com/tkphung)


## Acknowledgments

* [Rogge+ 2017](https://doi.org/10.1111/jmi.12593)
* [Zhang+ 2017](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1684-y).
