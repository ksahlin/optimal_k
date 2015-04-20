
import argparse
import sys
import warnings

try:
	import matplotlib
	matplotlib.use('agg')
	import matplotlib.pyplot as plt
	import numpy as np
except ImportError:
	print("Couldn''t import either numpy or matplotlib. You need tohave them installed. \
		Or manually look in the csv files generated to look for best combination of k and a.")
	sys.exit()




def main(csvfile, plot=False):
	lines = open(csvfile,'r').readlines()[1:]
	vals = [line.strip().split(',') for line in lines]
	points = []
	for val in vals:
		try:
			point = (int(val[0]), float(val[8]))
			points.append(point)
		except ValueError:
			pass

	#points = map(lambda x: ( int(x[0]) , float(x[8])), vals)
	try:

		points = np.array(points)
		# get x and y vectors
		x = points[:,0]
		y = points[:,1]
	except IndexError:
		print -1, -1
		return -1, -1

	# go through all polynomial fits and take the one with highest fit that 
	# does not throw a warning

	with warnings.catch_warnings():
		warnings.filterwarnings('error')
		for degree in range(2,40):
			try:
				# calculate polynomial
				z = np.polyfit(x, y, degree)
				f = np.poly1d(z)
			except Warning: 
				# print 'Warning raised for degree {0}'.format(degree)
				final_degree = degree -1
				break

	# print "final degree= ", final_degree

	# calculate polynomial
	z = np.polyfit(x, y, final_degree)
	f = np.poly1d(z)			

	# calculate new x's and y's
	x_new = np.linspace(x[0], x[-1], len(points))
	y_new = f(x_new)
	# print len(y_new)
	# print max(y_new)
	max_x = x_new[y_new.argmax()]

	if plot:
		plt.plot(x,y,'o', x_new, y_new)
		plt.xlim([x[0]-1, x[-1] + 1 ])
		plt.savefig(plot+'.png')

	print max_x, max(y_new)
	return max_x, max(y_new)
if __name__ == '__main__':

	parser = argparse.ArgumentParser()#prog="Infer variants with simple p-value test using theory of GetDistr - proof of concept.")
	parser.add_argument('csv', type=str, help='csv file with predictions. ')
	parser.add_argument('--plot', type=str, default=False, help=' prefix to plot with fitted curve')

	args = parser.parse_args()
	main( args.csv , args.plot)

