#!/usr/bin/env python

import pandas as pd

df = pd.DataFrame({	'num_legs': [2, 4, 8, 0], 
					'num_wings': [2, 0, 0, 0], 
					'num_specimen_seen': [10, 2, 1, 8]}, 
					index=['falcon', 'dog', 'spider', 'fish'])

df['num_legs'].sample(n=3, random_state=1)
df.sample(frac=0.5, replace=True, random_state=1)
df.sample(n=2, weights='num_specimen_seen', random_state=1)
