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
    vcf_final="1KG+HGDP.chr${chr}.hapmap.final.vcf"
    # download data and index files
    ${dltool} https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes/hgdp.tgp.gwaspy.merged.chr${chr}.merged.bcf
    ${dltool} https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes/hgdp.tgp.gwaspy.merged.chr${chr}.merged.bcf.csi
    # extract hapmap variants
    ${bcftools} view -O z -R ${variant_list} -o ${vcf_hapmap} ${bcf_in}
    # change format from GT:PS to just GT
    ${bcftools} annotate -x ^FORMAT/GT ${vcf_hapmap} -o ${vcf_final}
    # remove intermediate files
    rm ${bcf_in}
    rm ${bcf_in}.csi
    rm ${vcf_hapmap}
done
)

# TODO download population lists
(cd ${data_dir} && ${dltool} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt)
mv ${data_dir}/20130606_g1k_3202_samples_ped_population.txt ${data_dir}/1KGPhase3.pop.panel

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
