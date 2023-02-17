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

Warning: With the flag `-w`, the code tries to sort trajectory frames by cluster. Currently this only works if all trajectories have the same topology and otherwise causes the ccode to fail.


#### Output and Interpretation

The folder __results__ will contain the following subfolders:
 - 1-distances: numerical data of all distances between CA atoms in the binding pocket.
 - 2-comparison: plots that show the distributions along CA distances that deviate the most.
 - 3-pca: plots of the principal component analysis. Square symbols represent the starting structure (if the flag `--showstart` was set) and circles represent the average of each state (active or inactive).
 - 4-clustering: scatter plots for each clustering run and elbow plots for all clustering runs for each value of n. Circles represent cluster centers.
 - 5-rmsf: RMSF plots for all clusters in each clustering run. Residues in the binding pocket are highlighted. The RMSF is also stored in corresponding structure files for each cluster. We use the following tags for different representative structures: avg = average structure, cluster center (recommended for restraints), ref = cluster centroid, frame closest to the center, top = topology file, starting file of the simulation.
 - 6-sorted (only if the flag `-w` was set): Trajectory frames sorted by cluster (each cluster is one trajectory). Warning: The option `-w` currently only works if all trajectories have the same topology.

We can use the results to determine representative reference structures and widths of the restraints that avoid overlap between active and inactive ensembles.

The comparison of distances and the PCA show us how much the simulations started from active and inactive structures overlap and how much they deviate from the starting structure. If the starting structure is not within one of the dominant clusters, we have to think about whether this was caused by inaccurate modeling (protonation states, ions,...) or by artifacts in the starting structure (deformations caused by crystal contacts, stabilizing nanobodies,...).

The cluster plots show us which clustering parameters provide reasonable results. Unfortunately there is no strict criterion but a good clustering is usually at the "elbow" of the plot of the sum of the squared distances over the number of clusters and it should visibly "make sense" in the PCA plot. Choose a set of parameters n (number of PCA components used for the clustering) and k (number of clusters) and s (random seed) and remember them later.

Pick the RMSF plot for the set of clustering parameters that you picked from the clustering and check the RMSF, especially in the binding pocket. 
Compare it to the overlap in the plots from the distance comparison. 
Choose a factor by which the RMSF has to be rescaled such that the resulting distribution along a significant number of distances will not overlap (again, there is no strict rule for this). 
Ideally, the restraints will allow enough unrestrained motion to not perturb free energy calculation but not so much that the ensembles will overlap. 
For each cluster, use the corresponding "_avg.cms" file to define the restraints (see below for how to do that).


### Starting Structures and Ligand Modeling

Choose a simulation frame to start the ABFEP workflow in which the protein has sufficiently relaxed (especially for membrane proteins). As we use restraints, it does not matter to much which frame you use as long as it is not a crazy outlier. We provide a script to extract a specific frame from a set of simulations. For example, if you want to extract the tenth frame from all your template simulations, run
```
$SCHRODINGER/run ~/dev/abfep-restraints/extract_starting_frames.py -i example/input.csv -n 10
```
If you want to use a different frame for each trajectory, split up the input file and run the script separately for each trajectory.

The extracted frames can then be converted to "_pv.mae" files as needed for FEP using
```
$SCHRODINGER/run ~/dev/abfep-restraints/cms2fep.py frame.cms -o fep_pv.mae -ligand "{Ligand ASL}"
```
Then you can model the ligands into this structure or align the experimental structure in which you have already modeled them and copy them over.


### Adding Restraints to AB-FEP

This section assumes that the reference structure for the restraints is stored in a cms file and the width of the restraints in the field "sigma" of the corresponding atom. 

From the file with your modeled input poses, you can generate a job directory for a restrained AB-FEP job, for example:
```
$SCHRODINGER/run ~/dev/abfep-restraints/setup_restrained_abfep.py \
	input-poses_pv.maegz \
	-r pca-kmeans_n04_s42_k05_cluster00_rmsf_avg.cms \
	-j input-poses_n04_s42_k05_cluster00 \
	--scaling-factor 0.5 \
	--md-force-const 0.1 \
	--fep-force-const 1.0 \
	--align-sel 'protein and chain A and at.ptype CA and ( res. 54 - 75 or res. 86 - 148 or res. 165 - 240 or res. 254 - 282 or res. 298 - 320 )' \
	--ffbuilder \
	--host bolt_cpu \
	--subhost bolt_gpu \
	--project dev_GPU \
	--maxjob 0 \
	--retries 5 \
	--ffhost bolt_cpu \
	--md-sim-time 2000 \
	--fep-sim-time 10000
```
where the option `-j` determines the name of the directory and the job name. 
The option `-r` provides the restraint file and `--scaling-factor` the factor by which the width of the restraints from the cms file (e.g., from the RMSF as described above) is to be rescaled. 
The restraints file has to be aligned to the start structure and the flag `--align-sel` provides the ASL for the atoms at which to align them.
This selection must only include residues that are present in both the start file and the restraints file. It is the user's duty to make sure residue numbers, chain IDs etc. correspond to each other.

Then `cd` into the job directory and submit the job by running the bash script. 
You can also make manual adjustments there.


## Support and Contributing

For questions or ideas, open an issue, ping Martin Voegele (Desmond Team, NYC), or write an e-mail (martin.voegele@schrodinger.com).


## Acknowledgements

This repo includes a modified version of membrane_cms2fep.py. 
The script add_sigma_restraints.py uses code to add restraints to msj files by Mikolai Fajer.
The analysis workflow uses functions from PENSA (M. Voegele, MIT license).


## Project Status

This repo is still under development and being tested. It might change without warning.


