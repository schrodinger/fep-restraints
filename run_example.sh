# Download the example trajectory if it is not there
if [ -d example/Sim01 ]; then
	echo "Example trajectory detected."
else
	echo "No example trajectory detected. Attempting download."
	mkdir -p example
	scp -r boltio:/nfs/working/scidev/voegele/GPCR-Functional-Response/B2AR/MDSim/7DHI_Salbutamol_MDSim/Sim01 example/
fi

EXAMPLE='example/Sim01/7DHI_Salbutamol_MDSim'
TOP=${EXAMPLE}.cms
TRJ=${EXAMPLE}_trj

# Read desired C-alpha distances from a trajectory
echo "Extracting C-alpha distances."
mkdir -p example/results
CHN="R"
RES="109 113 114 117 118 193 199 203 204 207 286 289 290 293 308 312 316"
DIS="example/results/ca-dist_7DHI_Sim01.csv"
$S/run read_ca-distances_from_trajectory.py \
	-c ${TOP} -t ${TRJ} -l ${CHN} -r ${RES} -o ${DIS}

# Perform clustering in principal component space of the C-alpha distances
echo "Calculating clusters via k-means"
CLUSTERS="example/results/clusters_on_pca"
$S/run cluster_features_pca_k-means.py -i $DIS -o $CLUSTERS

# Write selected frames from a trajectory
echo "Extracting frames."
mkdir -p example/frames
$S/run extract_frames_as_cms.py -c ${EXAMPLE}.cms -t ${EXAMPLE}_trj -n 5 55 555 -o example/frames

# Convert the frames to MAE files suitable for FEP+
echo "Converting frames."
for FRAME in example/frames/*.cms; do
	OUT_FRAME=$(echo $FRAME | sed 's/\.cms//g' )
	$S/run membrane_cms2fep.py -o ${OUT_FRAME}_pv.mae $FRAME
done

