#### https://coptr.readthedocs.io/en/latest/installation.html

# Clone code from Git
git clone https://github.com/tyjo/coptr
cd coptr/

# Load conda
module load miniconda/3

# Set conda channel priority to flexible
conda config --show
conda config --env --set channel_priority flexible

# Creates an environment called coptr
mamba env create -f coptr.yml

# Reset conda channel priority to strict (default)
conda config --env --set channel_priority strict

# Activate the environment and install via pip
conda activate coptr
pip install .

# Test
coptr

