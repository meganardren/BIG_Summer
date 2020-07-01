import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool

fn = "coords_ScaleUpDesign1_hg19.txt"
# open middles_of_sharpr_regions_fn for writing, creates cell_type specific file
middles_of_sharpr_regions_fn = open("Hepg2.bed","w")
with open(fn) as f: 
        for line in f:
        # for each line in the file, split it up to its columns
            region_id, chrom, start, end = line.split('\t')
        # split up region_id
            cell_type, region_id_in_state, tile_pos, chrom, center = region_id.split('_')
        # ... now convert to bed and write to a different file
        #only add Hepg2 celltype to Hepg2.bed file
            if (cell_type == "Hepg2"): #select specific cell type
            #write to file in bed format
                L = [chrom, '\t', start, '\t', end]
                middles_of_sharpr_regions_fn.writelines(L)
# close middles_of_sharpr_regions_fn for writing
middles_of_sharpr_regions_fn.close()

#open middles_of_sharpr_regions_fn for reading
middles_of_sharpr_regions_fn = open("Hepg2.bed","r")
# open denase_bed_fn for reading
dnase_bed_fn = open("wgEncodeUwDnaseHepg2PkRep1.narrowPeak","r")

#create dnase_beds bedtool
dnase_beds = BedTool(dnase_bed_fn)
#create region_beds bedtool
region_beds = BedTool(middles_of_sharpr_regions_fn)
#create intersection bedtool
hits = region_beds.intersect(dnase_beds, c=True) # see pybedtools documentation
#make hits a dataframe
hits_df = hits.to_dataframe()

#close files
f.close()
middles_of_sharpr_regions_fn.close()
dnase_bed_fn.close()

print(hits_df)

true_hits = (hits_df['name'] == 1)
hits_df[true_hits]

perc_hits = len(hits_df[true_hits])/len(hits_df)
perc_hits