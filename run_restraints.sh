# Write constant restraints to MSJ.

DIR="example/b01f50_ABFEP_k3c0_bin"
SEL="protein AND backbone AND not a.ele H"

rm -rf $DIR
cp -r example/ABFEP_k3c0 $DIR

for MSJ in $DIR/*.msj; do 
	OLD=$( echo $MSJ | sed 's/\.msj/\.msj.orig/g' )
	mv $MSJ $OLD
	$S/run ~/dev/abfep-restraints/write_constant_restraints_to_msj.py \
		-m $OLD -a "$SEL" -t 'posre_fhbw' -k 50 -sigma=0.1 -o $MSJ
done


# Write RMSF-based restraints from a reference file to MSJ. 

DIR="example/b01f50_ABFEP_k3c0_mikolai"
SEL="protein AND backbone AND a.ptype CA"
REF="example/rmsf/clusters_on_pca_n03_s42_k04_cluster02_rmsf_avg.cms"

rm -rf $DIR
cp -r example/ABFEP_k3c0 $DIR

for MSJ in $DIR/*.msj; do
        OLD=$( echo $MSJ | sed 's/\.msj/\.msj.orig/g' )
        mv $MSJ $OLD
        $S/run ~/dev/abfep-restraints/write_restraints_from_mae_to_msj.py \
                $REF $OLD $MSJ -a "$SEL" -f 50 
done


# Write constant restraints from a reference file to MSJ. 

DIR="example/b01f50_ABFEP_k3c0_constant-ref"
SEL="protein AND backbone AND a.ptype CA"
REF="example/rmsf/clusters_on_pca_n03_s42_k04_cluster02_rmsf_avg.cms"

rm -rf $DIR
cp -r example/ABFEP_k3c0 $DIR

for MSJ in $DIR/*.msj; do
        OLD=$( echo $MSJ | sed 's/\.msj/\.msj.orig/g' )
        mv $MSJ $OLD
        $S/run ~/dev/abfep-restraints/write_restraints_from_mae_to_msj.py \
                $REF $OLD $MSJ -a "$SEL" -f 50 -b 0.1
done



