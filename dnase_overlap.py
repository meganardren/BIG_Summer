import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool

def compute_mipoint_overlap_percentage(offset,target_cell_type):
    fn = "coords_ScaleUpDesign1_hg19.txt"
    # open middles_of_sharpr_regions_fn for writing, creates cell_type specific file
    middles_of_sharpr_regions_fn = open("target_cell_type.bed","w")
    with open(fn) as f: 
            for line in f:
            # for each line in the file, split it up to its columns
                region_id, chrom, start, end = line.split('\t')
            # split up region_id
                cell_type, region_id_in_state, tile_pos, chrom, center = region_id.split('_')
            # ... now convert to bed and write to a different file
            #only add Hepg2 celltype to Hepg2.bed file
                if (cell_type == target_cell_type): #select specific cell type
                #write to file in bed format
                    #calculate middle position and assign start position)
                    int_start = int(start)
                    int_end = int(end)
                    start_pos = int_start+offset
                    #assign end position
                    end_pos = start_pos+1 
                    #write to file in bed format
                    start_pos = str(start_pos)
                    end_pos = str(end_pos)
                    L = [chrom, '\t', start_pos, '\t', end_pos,'\n']
                    middles_of_sharpr_regions_fn.writelines(L)
    # close middles_of_sharpr_regions_fn for writing
    middles_of_sharpr_regions_fn.close()

    #open middles_of_sharpr_regions_fn for reading
    middles_of_sharpr_regions_fn = open("target_cell_type.bed","r")
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
    #only ones that overlapped
    true_hits = (hits_df['name'] == 1)
    #calculate percentage of hits
    perc_hits = len(hits_df[true_hits])/len(hits_df)
    #return percent overlap
    return perc_hits

    #close files
    f.close()
    middles_of_sharpr_regions_fn.close()
    dnase_bed_fn.close()

   
    
if __name__ == "__main__":
    mid_hits_hepg2 = compute_mipoint_overlap_percentage(147,'Hepg2') #1.0
    mid_hits_h1hesc = compute_mipoint_overlap_percentage(147,'H1hesc') #0.13893129770992366
    mid_hits_k562 = compute_mipoint_overlap_percentage(147,'K562') #0.15572519083969466
    mid_hits_huvec = compute_mipoint_overlap_percentage(147,'Huvec') #0.16030534351145037
    
#     hepg2_295 = open("hepg2_295.txt","w")
#     for i in range(295):
#         hits_perc = compute_mipoint_overlap_percentage(i,'Hepg2')
#         hits_perc = str(hits_perc)
#         i = str(i)
#         hepg2_295.write(i)
#         hepg2_295.write('\t')
#         hepg2_295.write(hits_perc)
#         hepg2_295.write('\n')
#     hepg2_295.close()

# h1hesc_295 = open("h1hesc_295.txt","w")
# for i in range(295):
#     hits_perc = compute_mipoint_overlap_percentage(i,'H1hesc')
#     hits_perc = str(hits_perc)
#     i = str(i)
#     h1hesc_295.write(i)
#     h1hesc_295.write('\t')
#     h1hesc_295.write(hits_perc)
#     h1hesc_295.write('\n')
# h1hesc_295.close()

# k562_295 = open("k562_295.txt","w")
# for i in range(295):
#     hits_perc = compute_mipoint_overlap_percentage(i,'K562')
#     hits_perc = str(hits_perc)
#     i = str(i)
#     k562_295.write(i)
#     k562_295.write('\t')
#     k562_295.write(hits_perc)
#     k562_295.write('\n')
# k562_295.close()

# huvec_295 = open("huvec_295.txt","w")
# for i in range(295):
#     hits_perc = compute_mipoint_overlap_percentage(i,'Huvec')
#     hits_perc = str(hits_perc)
#     i = str(i)
#     huvec_295.write(i)
#     huvec_295.write('\t')
#     huvec_295.write(hits_perc)
#     huvec_295.write('\n')
# huvec_295.close()