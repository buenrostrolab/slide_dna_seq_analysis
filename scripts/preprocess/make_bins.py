import sys
import natsort
import math

chrom_sizes_file = open(sys.argv[1],'r')
bin_size = int(sys.argv[2])

print("\t".join(["chr", "bin_start", "bin_end", "bin_len", "chr_ind", "bin_ind"]))

# store chromosome sizes

chrom_sizes = {}

for line in chrom_sizes_file:

    column = line.rstrip().split()
    chrom = column[0]

    if "_" not in chrom: 
        chrom_sizes[chrom] = int(column[1])

# loop through chromosomes

chrom_ind = 1
bin_ind = 1
for chrom in natsort.natsorted(chrom_sizes):

    if chrom == "chrM":
        continue

    chrom_size = chrom_sizes[chrom]
    num_bins = int(math.ceil((float(chrom_size)/bin_size)))

    for i in range(num_bins):
        bin_start = max(1,i*bin_size)
        bin_end = min(chrom_size,(i+1)*bin_size)
        print("\t".join([chrom, str(bin_start), str(bin_end), str(bin_end-bin_start), str(chrom_ind), str(bin_ind)]))
        #print("\t".join([chrom, str(max(1,i*bin_size)), str(min(chrom_size,(i+1)*bin_size)), str(chrom_ind), str(bin_ind)]))
        bin_ind += 1

    chrom_ind += 1

# add chrM to the end

chrom = "chrM"

chrom_size = chrom_sizes[chrom]
num_bins = int(math.ceil((float(chrom_size)/bin_size)))

for i in range(num_bins):
    bin_start = max(1,i*bin_size)
    bin_end = min(chrom_size,(i+1)*bin_size)
    print("\t".join([chrom, str(bin_start), str(bin_end), str(bin_end-bin_start), str(chrom_ind), str(bin_ind)]))
    #print("\t".join([chrom, str(max(1,i*bin_size)), str(min(chrom_size,(i+1)*bin_size)), str(chrom_ind), str(bin_ind)]))
    bin_ind += 1
