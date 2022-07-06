# Download the example trajectory if it is not there
if [ -d example/Sim01 ]; then
	echo "Example trajectory detected."
else
	echo "No example trajectory detected. Attempting download."
	mkdir -p example
	scp -r boltio:/nfs/working/scidev/voegele/GPCR-Functional-Response/B2AR/MDSim/7DHI_Salbutamol_MDSim/Sim01 example/
fi

EXAMPLE='example/Sim01/7DHI_Salbutamol_MDSim'
TOP=${EXAMPLE}-out.cms
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

# Write clusters as trajectories
echo "Writing clusters as trajectories"
EX_NAME='clusters_on_pca_n03_s42_k04'
mkdir -p example/cl_traj
$S/run extract_clusters_as_trj.py \
	-c example/Sim01/7DHI_Salbutamol_MDSim-out.cms \
	-t example/Sim01/7DHI_Salbutamol_MDSim_trj \
	-n example/results/${EX_NAME}_ca-dist_7DHI_Sim01.csv \
	-o example/cl_traj/${EX_NAME}_ca-dist_7DHI \
	-d example/results/${EX_NAME}_summary.csv

# Extract centroids
$S/run extract_centroids.py \
        -c example/Sim01/7DHI_Salbutamol_MDSim.cms \
        -t example/Sim01/7DHI_Salbutamol_MDSim_trj \
        -o example/cl_traj/${EX_NAME}_ca-dist_7DHI \
        -d example/results/${EX_NAME}_summary.csv

echo "Converting centroids"
for FRAME in example/cl_traj/*centroid*.cms; do
       OUT_FRAME=$(echo $FRAME | sed 's/\.cms//g' )
       $S/run membrane_cms2fep.py -o ${OUT_FRAME}_pv.mae $FRAME
done

# Calculate RMSF of each cluster
SEL="(res.num ${RES}) AND (atom.ptype \" CA \") AND (chain.name ${CHN})"
CMS=example/cl_traj/${EX_NAME}_ca-dist_7DHI.cms
for CLUSTER in 00 01 02 03; do
	REF="example/cl_traj/${EX_NAME}_ca-dist_7DHI_centroid${CLUSTER}.cms"
	TRJ="example/cl_traj/${EX_NAME}_ca-dist_7DHI_cluster${CLUSTER}.xtc"
	OUT="example/cl_traj/${EX_NAME}_ca-dist_7DHI_cluster${CLUSTER}"
	$S/run calculate_rmsf_from_trajectory.py \
		-c "$CMS" -t "$TRJ" -s "$SEL" -o "$OUT"\
		--ref_file "$REF" --ref_sel "$SEL" 
done

# Write selected frames from a trajectory
#echo "Extracting frames."
#mkdir -p example/frames
#$S/run extract_frames_as_cms.py -c ${EXAMPLE}.cms -t ${EXAMPLE}_trj -n 5 55 555 -o example/frames

# Convert the frames to MAE files suitable for FEP+
#echo "Converting frames."
#for FRAME in example/frames/*.cms; do
#	OUT_FRAME=$(echo $FRAME | sed 's/\.cms//g' )
#	$S/run membrane_cms2fep.py -o ${OUT_FRAME}_pv.mae $FRAME
#done

