Scripts to quantify the fluorescence intensity of the Crbn regulatory sequences.

Z-stack images from Zeiss 880 confocal were save as maximum intensity projection using saveMaxProjectionfromCzi.ijm script.
They were also saved as sum projection using saveSumProjectionfromCzi.ijm script

Fluorescence intensity in the muscle nuclei was quantified using intensityMeasurement.ijm. A mask of the muscle nuclei was created using the max intensity projection. The mask was transfered to the sum projection image and used to quantify the total fluorescence intensity (integrated intensity).

The integrated fluorescence intensity was set to zero when no muscle nuclei was detected by the script intensityMeasurement.ijm. A visual inspection was also performed.
