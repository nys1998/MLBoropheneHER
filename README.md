# Machine Learning Predictor for Hydrogen Evolution Reaction on Borophene
Supplementary Data and Code file for the paper "Machine Learning with Interpretable Local Descriptors for Hydrogen Evolution Reaction on Borophene"

The input files required for the generation of the dataset is included in directories named after the borophene phases. The files are named after the strain type, for examples strain of 2% on the x-axis and -1% on the y-axis is named "a2b-1.in". To parse and generate the dataset, run GenerateDescriptor.py, which has dependencies on the functions in DescriptorFunction.py and Input_Parser.py.

If the user is just interested to run the K-Fold CV, the premade dataset "training9Sep.dat" is available. KFold.py is not dependent on any other code in this repository and can be run independently. 
