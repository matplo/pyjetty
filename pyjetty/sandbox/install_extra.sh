# https://github.com/plotly/plotly.py
#pip install jupyterlab==0.35 "ipywidgets>=7.2"
#pip install plotly==3.10.0
pip install "jupyterlab>=0.35" "ipywidgets>=7.2"
pip install "plotly>=3.10.0"

# Avoid "JavaScript heap out of memory" errors during extension installation
# (OS X/Linux)
export NODE_OPTIONS=--max-old-space-size=4096
# (Windows)
# set NODE_OPTIONS=--max-old-space-size=4096

# Jupyter widgets extension
jupyter labextension install @jupyter-widgets/jupyterlab-manager@0.38 --no-build

# FigureWidget support
jupyter labextension install plotlywidget@0.11.0 --no-build

# offline iplot support
jupyter labextension install @jupyterlab/plotly-extension@0.18.2 --no-build

# JupyterLab chart editor support (optional)
jupyter labextension install jupyterlab-chart-editor@1.1 --no-build

# Build extensions (must be done to activate extensions since --no-build is used above)
jupyter lab build

# To install ipympl with pip:
pip install ipympl
# If using JupyterLab
# Install nodejs: https://nodejs.org/en/download/
jupyter labextension install @jupyter-widgets/jupyterlab-manager
jupyter labextension install jupyter-matplotlib

# Unset NODE_OPTIONS environment variable
# (OS X/Linux)
unset NODE_OPTIONS
# (Windows)
# set NODE_OPTIONS=

# static image export
# using the to_image and write_image functions in the plotly.io package.
# This functionality requires the installation of the plotly orca command line utility and the psutil Python package.
pip install psutil

# https://github.com/plotly/orca
# npm install -g electron@1.8.4 orca
npm install -g electron@1.8.4 orca
