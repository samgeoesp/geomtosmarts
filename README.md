# GeomToSmarts

Recently, self-supervised approaches have utilised reaction smarts to generate a feature representation to be used for downstream ML tasks.

The purpose of this code us to convert 3D TS geometries from Gaussian optimisation files into reaction SMARTS. This ultimately will give access to 3D transition structure geometries across multiple papers to be quickly added to a machine readable format.

This code works by taking a TS structures that has an associated frequency (Gaussian .out file) and separates the reactants. Alongside these, the frequency is used to generate the product. These geometries are then converted into SMILES which are in turn converted to SMARTS. This code utilises functionality from DIASSEP to generate the structures.

*** Package still needs to be uploaded to the distribution archives. ***

