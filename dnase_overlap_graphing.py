import pandas as pd
import matplotlib.pyplot as plt

# data = pd.read_csv("hepg2_295.txt", sep='\t', names = ["base", "perc_overlap"])

# data.plot(kind='scatter',x='base',y='perc_overlap',color='red')
# plt.title("Hepg2 Overlap")
# plt.savefig("hepg2_overlap.png")

# data = pd.read_csv("h1hesc_295.txt", sep='\t', names = ["base", "perc_overlap"])

# data.plot(kind='scatter',x='base',y='perc_overlap',color='red')
# plt.title("H1hesc Overlap")
# plt.savefig("h1hesc_overlap.png")

# data = pd.read_csv("k562_295.txt", sep='\t', names = ["base", "perc_overlap"])

# data.plot(kind='scatter',x='base',y='perc_overlap',color='red')
# plt.title("K562 Overlap")
# plt.savefig("k562_overlap.png")

data = pd.read_csv("huvec_295.txt", sep='\t', names = ["base", "perc_overlap"])

data.plot(kind='scatter',x='base',y='perc_overlap',color='red')
plt.title("Huvec Overlap")
plt.savefig("huvec_overlap.png")