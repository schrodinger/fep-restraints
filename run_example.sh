# Download the example trajectory if it is not there
if [ -d example/Sim01 ]; then
	echo "Example trajectory detected."
else
	echo "No example trajectory detected. Attempting download."
	mkdir -p example
	scp -r boltio:/nfs/working/scidev/voegele/GPCR-Functional-Response/B2AR/MDSim/*_MDSim/Sim01 example/
fi
scp -r boltio:/nfs/working/scidev/voegele/GPCR-Functional-Response/B2AR/analysis/definitions/restraints_*.txt example/

# The systems to analyze
SYSTEMS="2RH1_Carazolol 3P0G_BI-167107 5JQH_Carazolol 7DHI_Salbutamol"

# Read desired C-alpha distances from a trajectory
echo "Extracting C-alpha distances."
mkdir -p example/results
for SYS in $SYSTEMS; do
	TOP=example/Sim01/${SYS}_MDSim-out.cms
	TRJ=example/Sim01/${SYS}_MDSim_trj
	CHN=$(grep $SYS example/restraints_bindingpocket_chains.txt | sed "s/${SYS}: //g")
	RES=$(grep $SYS example/restraints_bindingpocket_residues.txt | sed "s/${SYS}: //g")
	DIS="example/results/ca-dist_${SYS}.csv"
	$S/run read_ca-distances_from_trajectory.py \
		-c ${TOP} -t ${TRJ} -l ${CHN} -r ${RES} -o ${DIS}
done

# Perform clustering in principal component space of the C-alpha distances
echo "Calculating clusters via k-means"
ALL_DIS=""
for SYS in $SYSTEMS; do
	ALL_DIS="$ALL_DIS example/results/ca-dist_${SYS}.csv"
done
CLUSTERS="example/results/clusters_on_pca"
$S/run cluster_features_pca_k-means.py -i $ALL_DIS -o $CLUSTERS -w

# Write clusters as trajectories
echo "Writing clusters as trajectories"
mkdir -p example/cl_traj
EX_NAME='clusters_on_pca_n03_s42_k04'
for SYS in $SYSTEMS; do
	$S/run extract_clusters_as_trj.py \
	-c example/Sim01/${SYS}_MDSim-out.cms \
	-t example/Sim01/${SYS}_MDSim_trj \
	-n example/results/${EX_NAME}_ca-dist_${SYS}.csv \
	-o example/cl_traj/${EX_NAME}_ca-dist_${SYS} \
	-d example/results/${EX_NAME}_summary.csv
done

# Extract centroids
echo "Extracting centroids."
EX_NAME='clusters_on_pca_n03_s42_k04'
CMS=""
TRJ=""
DEF=example/results/${EX_NAME}_summary.csv
OUT=example/cl_traj/${EX_NAME}_ca-dist
for SYS in $SYSTEMS; do
	CMS="$CMS example/Sim01/${SYS}_MDSim-out.cms"
	TRJ="$TRJ example/Sim01/${SYS}_MDSim_trj"
done
$S/run extract_centroids.py -c $CMS -t $TRJ -d $DEF -o $OUT

# Create FEP+ input files from centroids
echo "Converting centroids."
for FRAME in example/cl_traj/*centroid*.cms; do
       OUT_FRAME=$(echo $FRAME | sed 's/\.cms//g' )
       $S/run membrane_cms2fep.py -o ${OUT_FRAME}_pv.mae $FRAME
done



# --- OLD ---
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

