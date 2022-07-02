# Download the example trajectory if it is not there
if [ -d example/Sim01 ]; then
	echo "Example trajectory detected."
else
	echo "No example trajectory detected. Attempting download."
	mkdir -p example
	scp -r boltio:/nfs/working/scidev/voegele/GPCR-Functional-Response/B2AR/MDSim/7DHI_Salbutamol_MDSim/Sim01 example/
fi

# Write selected frames from a trajectory
echo "Extracting frames."
EXAMPLE='example/Sim01/7DHI_Salbutamol_MDSim'
mkdir -p example/frames
$S/run read_frames_for_fep.py -c ${EXAMPLE}.cms -t ${EXAMPLE}_trj -n 5 55 555 -o example/frames

# Sanitize the frames
#echo "Sanitizing frames."
#$S/run sanitize_cms_for_reuse.py example/frames/*.cms

# Convert the frames to input files for FEP+
echo "Converting frames."
for FRAME in example/frames/*.cms; do
	OUT_FRAME=$(echo $FRAME | sed 's/\.cms//g' )
	$S/run membrane_cms2fep.py -o ${OUT_FRAME}_pv.mae $FRAME
done
