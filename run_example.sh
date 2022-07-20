BOLTDIR='/nfs/working/scidev/voegele/GPCR-Functional-Response/'
mkdir -p example

# Download the example trajectory if it is not there
if [ -d example/Sim01 ]; then
	echo "Example trajectory detected."
else
	echo "No example trajectory detected. Attempting download."
	scp -r boltio:$BOLTDIR/B2AR/MDSim/*_MDSim/Sim01 example/
fi

# Download the example simulation directory if it is not there
if [ -d example/ABFEP_k3c0 ]; then
        echo "Example simulation directory detected."
else
        echo "No example simulation directory detected. Attempting download."
	scp -r boltio:$BOLTDIR/B2AR/ABFEP_centroids_fcvar/ABFEP_k3c0 example/
fi

# Download all definitions
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
	$S/run read_ca_distances_from_trajectory.py \
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


EX_NAME='clusters_on_pca_n03_s42_k04'

mkdir -p example/rmsf-in-orig

# Calculate RMSF of each cluster
for CLUSTER in 00 01 02 03; do

	echo " -- Cluster $CLUSTER -- "

	# Get the reference system 
	REFTOP="example/cl_traj/${EX_NAME}_ca-dist_centroid${CLUSTER}.cms"
	ORIGIN=$( $S/run get_centroid_origin.py -d example/results/${EX_NAME}_summary.csv -c $CLUSTER )
	REFSYS=$( echo $ORIGIN | sed 's/ca-dist_//g' | sed 's/\.csv//g' )
	# Align to the binding pocket
	REFCHN=$(grep $REFSYS example/restraints_bindingpocket_chains.txt | sed "s/${REFSYS}: //g")
	REFRES=$(grep $REFSYS example/restraints_bindingpocket_residues.txt | sed "s/${REFSYS}: //g")
	REFSEL="(res.num $REFRES) AND (atom.ptype \" CA \") AND (chain.name $REFCHN)"

	# Get the trajectory file(s) for this cluster in the origin of its centroid and the corresponding topologies and selections
	ALLTRJ=""
	ALLTOP=""	
	for TRJ in example/cl_traj/${EX_NAME}_ca-dist_${REFSYS}_cluster${CLUSTER}.xtc; do
		# Get the corresponding trajectory and system name
		TOP=$( echo $TRJ | sed "s/_cluster${CLUSTER}.xtc/.cms/g" )
		SYS=$( echo $TRJ | sed "s/_cluster${CLUSTER}.xtc//g" | sed "s/example\/cl_traj\/${EX_NAME}_ca-dist_//g" )
		# Align to the binding pocket
		CHN=$(grep $SYS example/restraints_bindingpocket_chains.txt | sed "s/${SYS}: //g")
		RES=$(grep $SYS example/restraints_bindingpocket_residues.txt | sed "s/${SYS}: //g")
		echo "(res.num $RES) AND (atom.ptype \" CA \") AND (chain.name $CHN)" >> selections_align.tmp
		echo $(grep $SYS example/restraints_subset_selections.txt | sed "s/${SYS}: //g") >> selections_write.tmp
		ALLTRJ="$ALLTRJ $TRJ"
		ALLTOP="$ALLTOP $TOP"
	done

	# Run the RMSF calculation
	OUTSEL=$(grep $REFSYS example/restraints_subset_selections.txt | sed "s/${REFSYS}: //g")
	$S/run calculate_rmsf_from_trajectories.py \
		-c $ALLTOP -t $ALLTRJ -s selections_align.tmp -w selections_write.tmp \
		--ref_file $REFTOP --ref_sel_align "$REFSEL" --ref_sel_write "$OUTSEL" \
		-o example/rmsf-in-orig/${EX_NAME}_cluster${CLUSTER} --align_avg 

	rm selections_align.tmp
	rm selections_write.tmp

done

# Write RMSF-based restraints from a reference file to MSJ. 

DIR="example/b01f50_ABFEP_k3c0_restraints"
SEL="protein AND backbone AND a.ptype CA"
REF="example/rmsf-in-orig/clusters_on_pca_n03_s42_k04_cluster02_rmsf_avg.cms"

rm -rf $DIR
cp -r example/ABFEP_k3c0 $DIR

for MSJ in $DIR/*.msj; do
        OLD=$( echo $MSJ | sed 's/\.msj/\.msj.orig/g' )
        mv $MSJ $OLD
        $S/run ~/dev/abfep-restraints/write_restraints_from_mae_to_msj.py \
                $REF $OLD $MSJ -a "$SEL" -f 50 
done

