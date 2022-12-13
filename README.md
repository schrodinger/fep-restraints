# AB-FEP with Restraints

Tools to add restraints to AB-FEP runs. Includes a workflow to determine reference structures and restraint widths from plain MD simulations.

## Installation

Currently there is no installation package. The scripts assume that this repo has been cloned into `~/dev/`

```
mkdir -p ~/dev
cd ~/dev
git clone git@gitlab.schrodinger.com:martin.voegele/abfep-restraints.git
```

## Usage

This repo contains code to analyze plain MD simulation and set up restraints based on their analysis. 
If you already know the reference structure and parameters of the restraints, skip the next two sections and go straight to "Adding Restraints to AB-FEP"

### Analysis of Plain MD Simulations

We want to compare active and inactive states, find clusters in the joint ensemble of multiple simulations, and calculate the RMSF within each cluster.
Active and inactive ensembles are compared using their distributions along all C-alpha distances of binding-pocket residues. K-means clustering is then performed on the n most important principal components with variable number of clusters k.

#### Preparation

To start prepare the following input files:
- A CSV file with the name, topology file, trajectory file, and label for each simulation to be included in the analysis. For a template, see [input.csv](https://gitlab.schrodinger.com/martin.voegele/abfep-restraints/-/blob/main/example/input.csv). Column names and their order must match the template.
- A CSV file with the corresponding selections of the chain and residue numbers for the binding pocket (used in the clustering analysis) and a selection of atoms to which restraints should be applied. For a template, see [selections.csv](https://gitlab.schrodinger.com/martin.voegele/abfep-restraints/-/blob/main/example/selections.csv). Again, column names and their order must match the template. 

Make sure that the selection of the binding pocket leads to sets of corresponding atoms across all simulations! If the PDB structures are numbered correctly, the residue numbers should be the same (but this is not always the case). A good way to determine which residues to consider part of the binding pocket, select the ligand in a PDB structure (or in multiple structures of the same protein) and expand the selection by 4 or 5 Angstroms, list all residue numbers that are selected in any simulation, and remove numbers of residues that are not present in one or more simulations.  

#### Analysis Run

The main analysis is performed all in one command, for example run
```
$SCHRODINGER/run ~/dev/abfep-restraints/analyze_cluster_rmsf.py -i input.csv -s selections.csv -o results -n 3 4 -k $(seq 1 8) --showstart
```
to perform a full analysis with k-means clustering on the 3 most important principal components and then again on the 4 most important principal components, each for k from 1 to 8 clusters, and save everything in the folder __results__.

#### Output and Interpretation

The folder __results__ will contain the following subfolders:
 - 1-distances: numerical data of all distances between CA atoms in the binding pocket.
 - 2-comparison: plots that show the distributions along CA distances that deviate the most.
 - 3-pca: plots of the principla component analysis. Square symbols represent the starting structure (if the flag `--showstart` was set) and circles represent the average of each state (active or inactive).
 - 4-clustering: scatter plots for each clustering run and elbow plots for all clustering runs for each value of n. Circles represent cluster centers.
 - 5-rmsf: RMSF plots for all clusters in each clustering run. Residues in the binding pocket are highlighted.
 - 6-sorted (only if the flag `-w` was set)

### Starting Structures and Ligand Modeling

### Adding Restraints to AB-FEP

This section assumes that the reference structure for the restraints is stored in a cms file and the width of the restraints in the field "sigma" of the corresponding atom. 


## Support and Contributing

For questions or ideas, open an issue, ping Martin Voegele (Desmond Team, NYC), or write an e-mail (martin.voegele@schrodinger.com).

## Acknowledgements

This repo includes a modified version of membrane_cms2fep.py. 
The script add_sigma_restraints.py uses code to add restraints to msj files by Mikolai Fajer.
The analysis workflow uses functions from PENSA (M. Voegele, MIT license).

## Project Status

This repo is still under development and being tested. It might change without warning.


