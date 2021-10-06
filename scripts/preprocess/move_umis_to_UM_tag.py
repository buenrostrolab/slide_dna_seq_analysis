import sys
import pysam

bam_file = pysam.AlignmentFile(sys.argv[1], "rb")
out_file = pysam.AlignmentFile(sys.argv[2], "wb", template=bam_file)

for read in bam_file.fetch():

        header = read.query_name
        umi = header.split("+")[1]
        read.set_tag('UM',umi)
        out_file.write(read)

sys.exit()

bam = open(sys.argv[1], 'r')

for line in bam:

	if line[0] == "@":
		print line.rstrip()
		continue

	column = line.rstrip().split()
	umi = column[0].split("+")[1]
	column[0] = column[0].split("_")[0]

	column.append("UM:Z:"+umi)
	print "\t".join(column)
