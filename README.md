# [geomtosmarts](src/geomtosmarts.py)

Recently, self-supervised approaches have utilised reaction smarts to generate a feature representation to be used for downstream ML tasks.

The purpose of this code us to convert 3D TS geometries from Gaussian output files into reaction SMARTS. This ultimately will give access to 3D transition structure geometries across multiple papers to be quickly added to a machine readable format.

This code works by taking a TS structures that has an associated frequency (Gaussian .out file) and separates the reactants. Alongside these, the frequency is used to generate the product. These geometries are then converted into SMILES which are in turn converted to SMARTS. 

## Install 

***Package still needs to be uploaded to the distribution archives.***

## Usage:

### Command line

```python geomtosmarts.py -f ../examples/diels_alder.out```

### Imported Package

```smiles, smarts = geo('file.out', r=False, k=False, n=False)```

*More information in [example.py](src/example.py)*

### Options

```
-k           Argument for keeping the created .xyz files. Defaults as False to not keep them.
-n           Argument for saving with filenames. Defaults as False to not save them.
-r           Argument for saving Reaction SMARTS .png file. Defaults as False to not save the image.
```

