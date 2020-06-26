import itertools
import pandas as pd
import numpy as np
x = pd.read_csv("tablenorm_recenterends_HepG2_Rep1_20.txt", sep='\t',index_col=0)
xs = np.array(x)
xs = xs.flatten()

import matplotlib.pyplot as plt
plt1 = plt.scatter(
    x = xs[0:2240],
    y = xs[1:2241]
)

plt.title("consecutive tile comparison (same replicate)")
plt.xlabel("median normalized MPRA activity for tile T, rep1")
plt.ylabel("median normalized MPRA activity for tile T+1, rep2")

plt.savefig("graph2.png")