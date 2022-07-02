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
$S/run extract_frames_as_cms.py -c ${EXAMPLE}.cms -t ${EXAMPLE}_trj -n 5 55 555 -o example/frames

# Convert the frames to MAE files suitable for FEP+
echo "Converting frames."
for FRAME in example/frames/*.cms; do
	OUT_FRAME=$(echo $FRAME | sed 's/\.cms//g' )
	$S/run membrane_cms2fep.py -o ${OUT_FRAME}_pv.mae $FRAME
done

