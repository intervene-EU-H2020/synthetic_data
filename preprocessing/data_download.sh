#!/bin/bash

data_dir=data/inputs/raw/1KG+HGDP
bcftools=bcftools
variant_list_path=../../../../external/variant_lists

echo "INFO: downloading data to" $data_dir

# Create directory for the data
mkdir -p ${data_dir}

# Select download tool
dltool=""
if command -v axel &> /dev/null; then
    dltool="axel"
else
    echo "INFO: 'axel' tool not found. Use 'axel' for accelerated downloads."
    echo "INFO: More information: https://github.com/axel-download-accelerator/axel"
    echo "INFO: Falling back to 'wget'."
    dltool="wget"
fi

# Download 1KG+HGDP data
(cd ${data_dir};
for chr in $(seq 1 22); do
    variant_list="${variant_list_path}/hapmap_variant_list_tabsep_chr${chr}.txt" # contains list of hapmap variants we want to keep
    bcf_in="hgdp.tgp.gwaspy.merged.chr${chr}.merged.bcf"
    vcf_hapmap="1KG+HGDP.chr${chr}.hapmap.vcf.gz"
    vcf_tmp="1KG+HGDP.chr${chr}.hapmap.tmp.vcf.gz"
    vcf_final="1KG+HGDP.chr${chr}.hapmap.final.vcf"
    # download data and index files
    ${dltool} https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes/hgdp.tgp.gwaspy.merged.chr${chr}.merged.bcf
    ${dltool} https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes/hgdp.tgp.gwaspy.merged.chr${chr}.merged.bcf.csi
    # extract hapmap variants
    ${bcftools} view -O z --types snps -R ${variant_list} -o ${vcf_hapmap} ${bcf_in}
    # change format from GT:PS to just GT
    ${bcftools} annotate -x ^FORMAT/GT ${vcf_hapmap} -o ${vcf_tmp}
    # add ID field
    ${bcftools} annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' ${vcf_tmp} -o ${vcf_final}
    gzip ${vcf_final}
    # remove intermediate files
    rm ${bcf_in}
    rm ${bcf_in}.csi
    rm ${vcf_hapmap}
    rm ${vcf_tmp}
done
)

# Download genetic mapping for bp/cM conversion
(cd ${data_dir};
mkdir -p genetic_maps
cd genetic_maps
for chr in $(seq 1 22); do
    ${dltool} https://github.com/joepickrell/1000-genomes-genetic-maps/raw/master/interpolated_from_hapmap/chr${chr}.interpolated_genetic_map.gz
    gunzip chr${chr}.interpolated_genetic_map.gz
done
)

# Download age of mutation maps
(cd ${data_dir};
mkdir -p mutation_maps
cd mutation_maps
for chr in $(seq 1 22); do
    ${dltool} https://human.genome.dating/bulk/atlas.chr${chr}.csv.gz
    gunzip atlas.chr${chr}.csv.gz
done
)
