$SCHRODINGER/run ~/dev/fep-restraints/submit_to_graphdb.py \
	example/OPRD1-I_sim01-frame051_morphine-scaffold.maegz \
	graphdb-job-absolute-binding.yaml \
	example/pca-kmeans_n04_s42_k05_cluster00_rmsf_avg.cms \
	--asl 'protein and chain A and at.ptype CA and ( res. 54 - 75 or res. 86 - 148 or res. 165 - 240 or res. 254 - 282 or res. 298 - 320 )' \
	--sf 1.0 \
	--fc 10.0

