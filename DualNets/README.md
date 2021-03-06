# DualNets
Updated DualNets package (previously written by Yikai and Ruhollah): A gene expression-mRNA stability dualnetwork master regulator identifier. 
Based on the ARACNe (Algorithm for the Reconstruction of Accurate Cellular Networks) package form the Califano Lab, rewritten in python code.

Original paper: https://pubmed.ncbi.nlm.nih.gov/16723010/

Updates:
* Memory-safe methods. Previous version holds all the data in RAM, which will inevitably crash computer as network grows. 
* Added cpu parallel methods to speed up MI calculations. 
* Added binary search methods to speed up thresholding.
* Split MI calculation component from DPI cutoff to allow users to select appropriate thresholds to limit false positive rates.
* Removed R portion of code and replaced with python.




