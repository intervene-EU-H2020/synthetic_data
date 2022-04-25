#!/bin/bash

data_dir=$1
echo "INFO: downloading data to" $data_dir

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

# Create directory for the data
mkdir -p ${data_dir}

for chr in $(seq 1 22); do
    target_name="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf"
    (cd ${data_dir};
    if [ ! -f ${target_name}.gz ]; then 
        ${dltool} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/${target_name}.gz; 
    fi
    )
    (cd ${data_dir};
    mv ${target_name}.gz 1KGPhase3.chr${chr}.vcf.gz
    )
done

# Download population information for 1000 Genomes dataset
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
