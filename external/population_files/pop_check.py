import pandas as pd

data_1kg = pd.read_csv("1KGPhase3.pop.panel", delim_whitespace=True)
data_hgdp = pd.read_csv("igsr-human genome diversity project.tsv", sep='\t')
data_other = pd.read_csv("NIHMS798259-supplement-supp_datatable1.csv")

def get_sex(sex):
    if sex=='female' or sex=='XX':
        return 2
    elif sex=='male' or sex=='XY':
        return 1
    else:
        return None

data_1kg = data_1kg[['FamilyID', 'SampleID', 'FatherID', 'MotherID', 'Sex', 'Population', 'Superpopulation']]
data_hgdp = data_hgdp[['Sample name', 'Sex', 'Superpopulation name', 'Population name']].rename(columns={'Sample name':'SampleID', 'Superpopulation name':'Superpopulation', 'Population name':'Population'})
data_hgdp['Sex'] = data_hgdp['Sex'].apply(lambda x : get_sex(x))
data_other = data_other[['Sample ID (Illumina)', 'Region', 'Population ID', 'Genetic sex assignment']].rename(columns={'Sample ID (Illumina)':'SampleID', 'Region':'Superpopulation', 'Population ID':'Population', 'Genetic sex assignment':'Sex'})
data_other['Sex'] = data_other['Sex'].apply(lambda x : get_sex(x))
df = pd.concat([data_1kg, data_hgdp, data_other])

fp = open("../../data/inputs/processed/1KG+HGDP/1KG+HGDP.chr22.hapmap.final.recode.vcf")
samples = None
for i, line in enumerate(fp):
    if i == 13:
        samples = line.split()[9:]

fp.close()

len(samples)

mer_df = pd.merge(df, pd.DataFrame({'SampleID':samples}), how='right')

# pop_map = {'AFR':'Africa', 'Africa (HGDP)':'Africa', 'AMR': 'America', 'America (HGDP)':'America', 'Central South Asia (HGDP)':'Central South Asia', 'CentralAsiaSiberia':'Central Asia Siberia', 'EAS':'East Asia', 'EUR':'Europe', 'East Asia (HGDP)':'East Asia', 'EastAsia':'East Asia', 'Europe (HGDP)':'Europe', 'Middle East (HGDP)':'Middle East', 'Oceania (HGDP)':'Oceania', 'Oceania (SGDP),Oceania (HGDP)':'Oceania', 'SAS':'South Asia', 'SouthAsia':'South Asia', 'WestEurasia':'West Eurasia'}
# mer_df['Superpopulation'] = mer_df['Superpopulation'].apply(lambda x : pop_map[x] if x in pop_map else x)

#mer_df[mer_df['Superpopulation'].isna()]

mer_df.to_csv('merged_pop.tsv', index=None, sep='\t')
