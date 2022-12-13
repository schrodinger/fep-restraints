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

### Starting Structures and Ligand Modeling

### Adding Restraints to AB-FEP

This section assumes that the reference structure for the restraints is stored in a cms file and the width of the restraints in the field "sigma" of the corresponding atom. 


## Support and Contributing

For questions or ideas, open an issue, ping Martin Voegele (Desmond Team, NYC), or write an e-mail (martin.voegele@schrodinger.com).

## Acknowledgements

This repo includes a modified version of membrane_cms2fep.py. 
The script add_sigma_restraints.py repurposes code to add restraints to msj files by Mikolai Fajer.
The analysis workflow uses functions from PENSA (M. Voegele, MIT license).

## Project Status

This repo is still under development and being tested. It might change without warning.


