import itertools
import pandas as pd
import numpy as np
x = pd.read_csv("tablenorm_recenterends_HepG2_Rep1_20.txt", sep='\t',index_col=0)
y = pd.read_csv("tablenorm_recenterends_HepG2_Rep2_20.txt", sep='\t',index_col=0)
xs = np.array(x)
xs = xs.flatten()
ys = np.array(y)
ys = ys.flatten()

import matplotlib.pyplot as plt
plt1 = plt.scatter(
    x = xs[0:2241],
    y = ys[0:2241]
)

plt.title("replicate comparison (same tile)")
plt.xlabel("median normalized MPRA activity for tile T, rep1")
plt.ylabel("median normalized MPRA activity for tile T, rep2")

plt.savefig("graph1.png")