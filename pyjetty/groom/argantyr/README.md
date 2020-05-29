# On analysis of Argantyr

## Note: notable changes to some code

- read all branches of the event tree
https://github.com/matplo/pyjetty/blob/8e983165123f7feade3ccc887907dc083a294465/pyjetty/mputils/data_io.py#L82
- changed also (expose all branches to the user): https://github.com/matplo/pyjetty/blob/8e983165123f7feade3ccc887907dc083a294465/pyjetty/mputils/data_io.py#L18
- using this (modified) - all attributes except anything with 'Particle': https://github.com/matplo/pyjetty/blob/8e983165123f7feade3ccc887907dc083a294465/pyjetty/mputils/data_io.py#L57
- one can do this:
- https://github.com/matplo/pyjetty/blob/0347992ba31000595db0275e685871e36f727abe/pyjetty/groom/embed_groomers.py#L248
- https://github.com/matplo/pyjetty/blob/0347992ba31000595db0275e685871e36f727abe/pyjetty/groom/embed_groomers.py#L293

## centrality bins:

number of participant cut:
90pc+ < 6
80pc+ < 15
70pc+ < 31
60pc+ < 56
50pc+ < 90
40pc+ < 133
30pc+ < 186
20pc+ < 248
10pc+ < 325

number of charged particles:
90pc+ < 277
80pc+ < 555
70pc+ < 1158
60pc+ < 2272
50pc+ < 3945
40pc+ < 6239
30pc+ < 9293
20pc+ < 13245
10pc+ < 18467

- some code used here
- https://github.com/matplo/pyjetty/blob/master/pyjetty/groom/argantyr/centrality.py
- https://github.com/matplo/pyjetty/blob/master/pyjetty/groom/argantyr/centrality_percentiles.ipynb
- to make histograms/files input for the notebook: tdraw_cfg.py tdraw_centrality.cfg -r input_file:centrality_output.root
- the notebook should give the figure:
