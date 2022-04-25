#!/bin/bash

data_dir=$1

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
