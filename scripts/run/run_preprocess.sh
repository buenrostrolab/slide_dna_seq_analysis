# set parameters

#data_dir=data/human_colon_cancer_3_dna
#data_dir=data/human_colon_cancer_methyl
data_dir=data/dev
#sample=human_colon_cancer_methyl_4x
#sample=210212_36
sample=S1_36_S4
bins=reference/genomic_bins/hg19_1Mb_bins.txt
bin_size=1000000
bin_size_name='1Mb'

# generate count and frag files

python3 scripts/preprocess/get_sparse_counts.py $data_dir/$sample.bead_locations.csv $bins $bin_size $data_dir/$sample.reads.txt > $data_dir/$sample.sparse_counts_$bin_size_name.txt &
#python scripts/preprocess/get_sparse_frags.py $data_dir/$sample.bead_locations.csv $bins $bin_size $data_dir/$sample.reads.txt > $data_dir/$sample.sparse_frags_$bin_size_name.txt &
