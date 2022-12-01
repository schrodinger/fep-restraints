rm -rf input_by_ligand
rm -rf OPRD1-I_sim01-frame051_morphine-scaffold_restrained

$S/run individual_input.py example/OPRD1-A_sim01-frame051_morphine-scaffold.maegz input_by_ligand 

$S/run setup_restrained_abfep.py \
	example/OPRD1-I_sim01-frame051_morphine-scaffold.maegz \
	-r example/pca-kmeans_n04_s42_k05_cluster00_rmsf_avg.cms \
	-j OPRD1-I_sim01-frame051_morphine-scaffold_restrained \
	--scaling-factor 0.5 \
	--md-force-const 0.1 \
	--fep-force-const 1.0 \
	--align-sel 'protein and chain A and at.ptype CA and ( res. 54 - 75 or res. 86 - 148 or res. 165 - 240 or res. 254 - 282 or res. 298 - 320 )' \
	--ffbuilder \
	--host bolt_cpu \
	--subhost bolt_gpu \
	--project dev_GPU \
	--maxjob 0 \
	--retries 5 \
	--ffhost bolt_cpu \
	--md-sim-time 2000 \
	--fep-sim-time 10000

#--oplsdir example/Opioid-Receptor-Ligands_FFB_oplsdir/custom_2022_1.opls \

