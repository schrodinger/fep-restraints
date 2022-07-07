
EX_NAME='clusters_on_pca_n03_s42_k04'

mkdir -p example/rmsf

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

	# Get all trajectory files for this cluster and the corresponding topologies and selections
	ALLTRJ=""
	ALLTOP=""	
	
	for TRJ in example/cl_traj/${EX_NAME}_ca-dist_*_cluster${CLUSTER}.xtc; do
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
	$S/run calculate_rmsf_from_trajectories.py \
		-c $ALLTOP -t $ALLTRJ \
		-s selections_align.tmp \
		-w selections_write.tmp \
		--ref_file $REFTOP --ref_sel "$REFSEL" \
		-o example/rmsf/${EX_NAME}_cluster${CLUSTER} 

	rm selections_align.tmp
	rm selections_write.tmp

done
