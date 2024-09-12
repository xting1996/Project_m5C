#####################################################
##xiaotingzhang@pku.edu.cn
##2019-12-29
#####################################################
##import
import collections
import argparse
import sys
parser = argparse.ArgumentParser(description="Add CCA tail to tRNA")
parser.add_argument("-Input_raw_tRNA","--raw_tRNA",
		nargs="?",type=argparse.FileType('r'),default=sys.stdin,
		help="Input rae tRNA file which needs to add CCA tail")
parser.add_argument("-Output_CCA_tRNA","--CCA_tRNA",
		nargs="?",type=argparse.FileType("w"),default=sys.stdout,
		help="Output added CCA tRNA file")
args = parser.parse_args()

input_file = args.raw_tRNA
output_file = args.CCA_tRNA

seq=collections.OrderedDict()

for line in input_file:
	if line.startswith('>'):
		#name = "_".join(line.split())
		name = line.split()[0]
		seq[name]=""
	else:
		seq[name] += line.replace("\n","").replace("U","T")

#output = open(output_file,'wb')

for l in seq.keys():
	seq[l] = seq[l]+"CCA"
	output_file.write(l+"\n"+seq[l]+"\n")
