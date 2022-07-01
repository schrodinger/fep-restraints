if [ -d example/Sim01 ]; then
	echo "Example trajectory detected."
else
	echo "No example trajectory detected. Attempting download."
	mkdir -p example
	scp -r boltio:/nfs/working/scidev/voegele/GPCR-Functional-Response/B2AR/MDSim/7DHI_Salbutamol_MDSim/Sim01 example/
fi

