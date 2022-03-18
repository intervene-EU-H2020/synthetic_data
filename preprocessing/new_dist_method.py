import pandas as pd

for chr in range(1, 22+1):
    rsid_map = 'data/inputs/processed/rsid_map_list_chr{}.txt'.format(chr)
    genetic_map = 'data/inputs/raw/genetic_maps/chr{}.interpolated_genetic_map'.format(chr)
    rsid_df = pd.read_csv(rsid_map, sep='\t')
    genetic_df = pd.read_csv(genetic_map, header=None, delim_whitespace=True).rename(columns={0:'rsid', 1:'bp', 2:'cm'})
    dist_df = pd.merge(rsid_df, genetic_df, on='rsid', how='left')
    dist_df = dist_df[['id_hg38', 'cm']].rename(columns={'id_hg38':'Variant', 'cm':'Distance'})
    dist_df['Index'] = range(1, len(dist_df)+1)
    dist_df = dist_df[['Index', 'Variant', 'Distance']]
    dist_df['Distance'] = dist_df['Distance'].interpolate()
    dist_df['Distance'] = dist_df['Distance'].backfill()
    dist_df.to_csv('data/inputs/processed/1KGPhase3.chr{}.hapmap.distfile'.format(chr), index=None)