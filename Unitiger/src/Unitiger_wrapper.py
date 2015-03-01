import argparse
import subprocess
import os

# set input parameters
parser = argparse.ArgumentParser(description="A wrapper for running Unitiger for more values of k and abundance. \n Ver. 1.0")

group = parser.add_argument_group("Required arguments")
group.add_argument("-u", dest="unitigerPath", help="path to Unitiger executable", required=True)
group.add_argument("-r", dest="readFile", help="FASTQ file(s) (separated by comma)", required=True)
group.add_argument("-o", dest="outputFile", help="output file", required=True)
group.add_argument("-k", dest="mink", help="min k", default=5)
group.add_argument("-K", dest="maxk", help="max k", default=5)

group = parser.add_argument_group("Optional arguments")
group.add_argument("-a", dest="mina", help="min a", default=3)
group.add_argument("-A", dest="maxa", help="max A", default=3)
group.add_argument("-t", dest="threads", help="number of threads (0 = all cores)", default = 1)

args = parser.parse_args()

mink = int(args.mink)
maxk = int(args.maxk)
mina = int(args.mina)
maxa = int(args.maxa)

metricsFile = file(args.outputFile + ".mink" + str(mink) + ".maxk" + str(maxk) + ".metrics.csv" , "w")
writtenHeader = False

for k in range(mink, maxk + 1):
	for a in range(mina, maxa + 1):
		print "RUNNING: " + args.unitigerPath + " -r " + args.readFile + " -o " + args.outputFile + " -k " + str(k) + " -a " + str(a) + " -t " + args.threads + " -s "
		subprocess.call(args.unitigerPath + " -r " + args.readFile + " -o " + args.outputFile + " -k " + str(k) + " -a " + str(a) + " -t " + args.threads + " -s ", shell=True)
		
		currentMetricsFile = file(args.outputFile + ".k" + str(k) + ".a" + str(a) + ".metrics.csv" , "r")
		lines = currentMetricsFile.readlines()
		if not writtenHeader:
			metricsFile.write(lines[0])
			writtenHeader = True
		metricsFile.write(lines[1])
		metricsFile.flush()
		currentMetricsFile.close()
		
		os.remove(args.outputFile + ".k" + str(k) + ".a" + str(a) + ".metrics.csv")

metricsFile.close()


