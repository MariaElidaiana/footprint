export CONDA_DIR=/cvmfs/des.opensciencegrid.org/fnal/anaconda2/
export PATH=$CONDA_DIR/bin:$PATH
source activate $CONDA_DIR/envs/default

echo "Environment is ready."
echo "We will launch a jupyter notebook no browser."
echo "To access the notebook page, create an ssh tunel in your local machine using"
echo "ssh -N -f -L 8888:localhost:8888 des41.fnal.gov &"
echo "then point your browser to"
echo "http://localhost:8888"

jupyter notebook --no-browser

 