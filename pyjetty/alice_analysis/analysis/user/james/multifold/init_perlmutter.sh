#!/usr/bin/bash

# Get command line option to determine whether we need to install the virtual environment, or just enter it
for i in "$@"; do
  case $i in
    --install)
      INSTALL=TRUE
      shift
      ;;
    -*|--*)
      echo "Unknown option $i"
      exit 1
      ;;
    *)
      ;;
  esac
done
if [ ! -z ${INSTALL} ]; then
    echo
    # Install miniconda
    #wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    #bash Miniconda3-latest-Linux-x86_64.sh -b
    source /global/homes/j/jdmull/miniconda3/bin/activate

    # Create virtual env and install conda root: https://github.com/conda-forge/root-feedstock/
    conda remove -n myenv --all
    conda create --name myenv
    conda activate myenv
    conda install root=6.26.10=py39h4661b88_1 --channel conda-forge

    # Install additional packages
    pip install --user \
    keras-tuner==1.1.2 \
    matplotlib==3.5.1 \
    networkx==2.7.1 \
    numpy==1.23.5 \
    pandas==1.4.1 \
    pyyaml==6.0 \
    scikit-learn==1.0.2 \
    seaborn==0.11.2 \
    silx==1.1.0 \
    tensorflow==2.11.0 \
    torch==1.11 \
    torch-geometric==2.0.4 \
    torch-scatter==2.0.9 \
    torch-sparse==0.6.13 \
    uproot==4.3.7
else
    # Activate virtual env
    source /global/homes/j/jdmull/miniconda3/bin/activate
    conda activate myenv
fi

export OUTPUT_DIR='/global/cfs/cdirs/alice/jdmull/multifold'