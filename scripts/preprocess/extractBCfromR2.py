import sys
from itertools import izip # only works in python 2, better for large files because it doesn't create a list
import gzip

r1 = gzip.open(sys.argv[1], 'rb')
r2 = gzip.open(sys.argv[2], 'rb')
r3 = gzip.open(sys.argv[3], 'rb')
o1 = gzip.open(sys.argv[4], 'wb')
o2 = gzip.open(sys.argv[5], 'wb')
o3 = gzip.open(sys.argv[6], 'wb')

line_count = 0

r1_all = ""
r3_all = ""

for line1, line2, line3 in izip(r1, r2, r3):

    line_count+=1

    r1_all += line1
    r3_all += line3

    # dumb fastq reader

    if line_count % 4 == 1:       # sequence header
        r2_1 = line2.rstrip()
    elif line_count % 4 == 2:     # sequence
        r2_2 = line2.rstrip()
    elif line_count % 4 == 3:     # quality header
        r2_3 = line2.rstrip()
    elif line_count % 4 == 0:     # quality
        r2_4 = line2.rstrip()

        # extract barcode

        idx = r2_2.find('ACGCTGAAGA')

        if idx > 15 and idx < 24:
            
            r2_2_fix = r2_2[idx-15:idx-15+7]
            r2_4_fix = r2_2[idx-15:idx-15+7]
            
            r2_2_fix += r2_2[idx+10:idx+10+7]
            r2_4_fix += r2_4[idx+10:idx+10+7]

            if len(r2_2_fix) == 14:

                print >> o1, r1_all.rstrip()
                print >> o3, r3_all.rstrip()

	        print >> o2, r2_1
	        print >> o2, r2_2_fix
	        print >> o2, r2_3
	        print >> o2, r2_4_fix

        r1_all = ""
        r3_all = ""

sys.exit()

