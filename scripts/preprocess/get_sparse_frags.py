import sys
import natsort
import gzip
import math

barcode_file = open(sys.argv[1],'r')
bin_file = open(sys.argv[2],'r')
bin_size = int(sys.argv[3])
frag_file = open(sys.argv[4],'r')
#frag_file = gzip.open(sys.argv[4],'rb')

complement = {"A":"T", "C":"G", "G":"C", "T":"A"}

barcodes = {}
bins = {}
counts = {}

barcode_ind = 1

def reverse_complement(seq):

    out = ""
    rev = seq[::-1]
    for i in range(len(rev)):
        out += complement[rev[i]]
    return out

# load barcodes

for line in barcode_file:

    if "barcodes" in line:
        continue

    column = line.rstrip().split(",")
    barcode = reverse_complement(column[0])
    barcodes[barcode] = barcode_ind
    barcode_ind += 1

# load bins

for line in bin_file:

    if "bin_start" in line:
        continue

    column = line.rstrip().split() 
    chrom = column[0]
    bin_start = int(column[1])
    bin_end = int(column[2])
    bin_ind = int(column[5])

    bins[chrom+":"+str(int(math.ceil(float(bin_end)/bin_size)))] = bin_ind 

for line in frag_file:

    column = line.rstrip().split()
    chrom = column[0]
    frag_start = int(column[1])
    barcode = column[2]
    
    if barcode not in barcodes:
        continue

    barcode_ind = barcodes[barcode]
 
    bin_start = int(math.ceil(float(frag_start)/bin_size))
    
    if chrom+":"+str(bin_start) not in bins:
        continue
    
    start_ind = bins[chrom+":"+str(bin_start)]

    if (barcode_ind,start_ind) not in counts:
        counts[(barcode_ind,start_ind)] = 1
    else:
        counts[(barcode_ind,start_ind)] += 1

for (barcode_ind, bin_ind) in sorted(counts):

    count = int(counts[barcode_ind,bin_ind])
    for i in range(count):
        print("\t".join([str(barcode_ind), str(bin_ind), "1"]))





