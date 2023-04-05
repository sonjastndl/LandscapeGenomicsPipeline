import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import csv
from csv import reader, writer
import numpy as np

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################
parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--samples", dest="SAMP", help="samplename file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

samples=[]
with open(options.SAMP, 'r') as f:
    reader = csv.reader(f)
    for line in reader:
        print(line)
        samples.extend(line)


with open(options.IN, 'r', encoding='latin-1') as f:
    reader = csv.reader(f, delimiter=',')
    data = [row[29:] for i, row in enumerate(reader) if i > 0 and row[0] in samples]

data_transposed = np.transpose(data)

with open(options.OUT, 'w', newline='') as f:
    writer = csv.writer(f, delimiter=' ')
    writer.writerows(data_transposed)
