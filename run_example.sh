BOLTDIR='/nfs/working/scidev/voegele/GPCR-Functional-Response/'
mkdir -p example

# Download the example trajectory if it is not there
if [ -d example/Sim01 ]; then
	echo "Example trajectory detected."
else
	echo "No example trajectory detected. Attempting download."
	scp -r boltio:$BOLTDIR/B2AR/MDSim/*_MDSim/Sim01 example/
fi

$S/run analyze_cluster_rmsf.py -i example/input.csv -s example/selections.csv -k 3 4 -w --showstart --step 5 --rmsf-across-all-trajectories

