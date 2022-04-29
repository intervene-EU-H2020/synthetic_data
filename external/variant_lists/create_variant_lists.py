import pandas as pd

path = "external/variant_lists/1KGPhase3_hm3_hg19_hg38_mapping_cached.tsv"

df = pd.read_csv(path, sep='\t')

# preprocessing (for extracting variants from raw data)

for i in range(22):
    chr = i+1
    tmp1 = df[df['chr']==chr][['chr', 'pos_hg38', 'a1', 'a2']]
    tmp1['chr'] = 'chr' + tmp1['chr'].astype(str)
    tmp2 = tmp1.copy()
    tmp2['a1'] = tmp1['a2']
    tmp2['a2'] = tmp1['a1']
    tmp = pd.concat([tmp1, tmp2]).sort_values(by='pos_hg38').drop_duplicates()
    tmp.to_csv('external/variant_lists/hapmap_variant_list_tabsep_chr{}.txt'.format(chr), sep='\t', header=None, index=None)

# postprocessing (for creating lists from extracted data)

datapath = "data/inputs/raw/1KG+HGDP/1KG+HGDP.chr{}.hapmap.final.vcf"

def get_snpids(datafile):
    snps = []
    with open(datafile) as fp:
        lines = fp.readlines()
        for line in lines:
            if not line.startswith('#'):
                snps.append(line.split()[2])
    return snps

num_variants = 0
for i in range(22):
    chr = i+1
    snpids = get_snpids(datapath.format(chr))
    snplist = pd.read_csv('external/variant_lists/hapmap_variant_list_tabsep_chr{}.txt'.format(chr), sep='\t', header=None)
    snplist['id'] = snplist[0] + ':' + snplist[1].astype(str) + ':' + snplist[2] + ':' + snplist[3]
    snplist = snplist[snplist['id'].isin(snpids)].sort_values(by=1)[['id']]
    num_variants+=len(snplist)
    # write variant list to file
    with open('external/variant_lists/hapmap_variant_list_chr{}.txt'.format(chr), 'w') as f:
        for item in snplist['id'].tolist():
            f.write("%s\n" % item)
    # map to rsid
    tmp1 = df[df['chr']==chr][['chr', 'pos_hg38', 'a1', 'a2', 'rsid']]
    tmp1['chr'] = 'chr' + tmp1['chr'].astype(str)
    tmp2 = tmp1.copy()
    tmp2['a1'] = tmp1['a2']
    tmp2['a2'] = tmp1['a1']
    tmp = pd.concat([tmp1, tmp2]).sort_values(by='pos_hg38').drop_duplicates()
    tmp['id'] = tmp['chr'] + ':' + tmp['pos_hg38'].astype(str) + ':' + tmp['a1'] + ':' + tmp['a2']
    tmp = tmp[['rsid', 'id', 'pos_hg38']]
    snplist = pd.merge(snplist, tmp, how='left').sort_values(by='pos_hg38').drop(columns=['pos_hg38'])
    assert len(snplist[snplist['rsid'].isna()])==0
    # check how many rsids matched
    assert len(tmp['rsid'].unique()) >= len(snplist['rsid'].unique())
    print("chr{}:\t{}% of hapmap3 variants retained ({} excluded)".format(chr, round(len(snplist['rsid'].unique())/len(tmp['rsid'].unique())*100,2), len(tmp['rsid'].unique())-len(snplist['rsid'].unique())))
    snplist = snplist.rename(columns={'id':'id_hg38'})
    snplist.to_csv('external/variant_lists/rsid_map_list_chr{}.txt'.format(chr), sep='\t', index=None)

print(num_variants)