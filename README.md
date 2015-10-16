# Photometry
This script converts the PTFIDE files into magnitudes and upperlimits,
making a number of corrections. These are based on based on the document
"http://web.ipac.caltech.edu/staff/fmasci/home/miscscience/forcedphot.pdf",
which is to be used in the PTFIDE format files produced from August 2015
and onwards, i.e. from version "forcepsffitdiff.pl v3.0".

The corrections included are:
- A correction for residual offset in the historical BASELINE. This is done
  so that any potential flux from the transient in the reference image is
  corrected for and be able to get the most accurate photometry. This is
  only important for transients. By historical we define 30 epochs earlier
  than the peak (here defined as the maximum flux measurement).
- UNCERTAINTY VALIDATION: two methods
- COURSE CHECK
- QUALITY CHECK


Corrections not included (Systematics):
- Incorrect PSF-template estimation
- Photometric zero-point calibrations
- Astrometric calibrations(determining PSF-placement)
- Error in supplied target position.

Disclaimer: This script is optimised for Ia supernova science and should
be changed for other science goals, especially for non-transients.

TO run the code: 
Put your PTFIDE code in a folder "forced_new" in this version of the code. Please also modify the code to accomondate your needs. 
