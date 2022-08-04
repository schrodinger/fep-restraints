BOLTDIR='/nfs/working/scidev/voegele/GPCR-Functional-Response/'
mkdir -p example

# Download the example trajectory if it is not there
if [ -d example/Sim01 ]; then
	echo "Example trajectory detected."
else
	echo "No example trajectory detected. Attempting download."
	scp -r boltio:$BOLTDIR/B2AR/MDSim/*_MDSim/Sim01 example/
fi

$S/run extract_starting_frames.py -i example/input.csv -n 10
 

