import argparse
import subprocess
import os

# set input parameters
parser = argparse.ArgumentParser(description="A wrapper for running Unitiger for more values of k and abundance. \n Ver. 1.0")

group = parser.add_argument_group("Required arguments")
group.add_argument("-r", dest="readFileName", help="FASTQ file(s) (separated by comma)", required=True)
group.add_argument("-o", dest="outputFileName", help="output file", required=True)

group = parser.add_argument_group("Optional arguments")
group.add_argument("-a", dest="mina", help="min a", default=3)
group.add_argument("-A", dest="maxa", help="max A", default=3)
group.add_argument("-k", dest="mink", help="min k", default=15)
group.add_argument("-K", dest="maxk", help="max k", default=0)
group.add_argument("-t", dest="threads", help="number of threads (0 = all cores)", default = 0)

args = parser.parse_args()

mink = int(args.mink)
maxk = int(args.maxk)
mina = int(args.mina)
maxa = int(args.maxa)

unitigerPath = str(os.path.realpath(__file__)).rpartition('/')[0] + "/Unitiger"

if maxk == 0:
	readFile = file(args.readFileName.split(",")[0], "r")
	line = readFile.readline()
	line = readFile.readline()
	maxk = len(line) - 10

metricsFile = dict()
writtenHeader = dict()

for a in range(mina, maxa + 1):
	metricsFile[a] = file(args.outputFileName + ".mink" + str(mink) + ".maxk" + str(maxk) + ".a" + str(a) + ".metrics.csv" , "w")
	writtenHeader[a] = False

for k in range(mink, maxk + 1):
	for a in range(mina, maxa + 1):
		print "RUNNING: " + unitigerPath + " -r " + args.readFileName + " -o " + args.outputFileName + " -k " + str(k) + " -a " + str(a) + " -t " + str(args.threads) + " -s "
		subprocess.call(unitigerPath + " -r " + args.readFileName + " -o " + args.outputFileName + " -k " + str(k) + " -a " + str(a) + " -t " + str(args.threads) + " -s ", shell=True)
		
		currentMetricsFile = file(args.outputFileName + ".k" + str(k) + ".a" + str(a) + ".metrics.csv" , "r")
		lines = currentMetricsFile.readlines()
		if not writtenHeader[a]:
			metricsFile[a].write(lines[0])
			writtenHeader[a] = True
		metricsFile[a].write(lines[1])
		metricsFile[a].flush()
		currentMetricsFile.close()
		
		os.remove(args.outputFileName + ".k" + str(k) + ".a" + str(a) + ".metrics.csv")

for a in range(mina, maxa + 1):
	metricsFile[a].close()


