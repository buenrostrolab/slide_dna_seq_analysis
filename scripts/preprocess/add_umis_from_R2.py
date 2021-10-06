import sys
from itertools import izip # only works in python 2, better for large files because it doesn't create a list
import gzip

r1 = gzip.open(sys.argv[1], 'rb')
r2 = gzip.open(sys.argv[2], 'rb')
r3 = gzip.open(sys.argv[3], 'rb')

o1 = gzip.open(sys.argv[4], 'wb')
o2 = gzip.open(sys.argv[5], 'wb')

flag = 0

for line1, line2, line3 in izip(r1, r2, r3):
   
	#print line1, line2, line3

	# tricky to get all headers but filter quality score lines that begin with "@", may need to adjust
	if line1[0:1] == "@" and " " in line1:

		l1 = line1.rstrip().split()
		l2 = line2.rstrip().split()
		l3 = line3.rstrip().split()

		if l1[0] == l2[0] and l2[0] == l3[0]:
			flag = 1

	elif flag == 1:

		umi = line2.rstrip()

		print >> o1, "_".join(l1) + "+" + umi
		print >> o2, "_".join(l3) + "+" + umi

		print >> o1, line1.rstrip()
		print >> o2, line3.rstrip()


		flag = 0

	else:
		print >> o1, line1.rstrip()
		print >> o2, line3.rstrip()
