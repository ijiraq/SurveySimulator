#!/usr/bin/env python
"""Compute the probabilty of the Anderson-Darling statistic
given a model distribution and sample.

This program uses a bootstrap approach to compute the probability
disrtibution of the A-D statistic for a given model distribution.  The
result returned is 1-CDF for the A-D statistic of the sample.  Values
of 0.01 reject the hypothesis that the sample is from the model at the
99% confidence level.

ad.py --help to see command line arguements and usage."""

import sys,os,logging
import optparse, re
import numpy

def ad_bootstrap(s:numpy.array, m:numpy.array, n_bootstraps:int = 1000) -> float:
	"""See if s is from m via AD using a bootstrap.

	s: sorted ndarray of sample data
	m: sorted ndarray of model data
	n_bootstraps: number of bootstrap iterations to use (default 1000)
	"""

	rng = numpy.random.default_rng()
	if m[-1]<s[-1] or s[0]<m[0]:
		raise ValueError("Model and sample do not overlap")
	n_sample_pts = len(s)
	sdf = numpy.array([ad(rng.choice(m, n_sample_pts), m) for i in range(n_bootstraps)])
	sdf.sort()
	x = ad(s, m)
	p = float((sdf>x).sum())/len(sdf)
	return p

def ad(s:numpy.array, m:numpy.array) -> float:
	"""Compute A-D for s from m where s and m are sorted arrays"""
	S = 0
	n_sample_pts = len(s)
	for j in range(n_sample_pts):
		s1 = numpy.log(cdf(s[j],m))
		s1 += numpy.log(1-cdf(s[n_sample_pts-1-j],m))
		s1 = (2*(j+1)-1)*s1
		S+=s1
	S /= n_sample_pts
	return -n_sample_pts-S


def cdf(x:float, y:numpy.array) -> float:
	"""Given a sorted array of values (y) compute the CDF at value (x)"""
	n = len(y)
	return (n-(y>x).sum())/n

	
if __name__ == "__main__":

	parser = optparse.OptionParser()
	parser.add_option("--col1",action="store", type="int", default=1,
			  help="Column in model data to act as source, default is 1")
	parser.add_option("--col2",action="store", type="int", default=1,
			  help="Column in observed data to act as source default is 1")
	parser.add_option("--nboots",action="store",type="int",default=500,
			  help="Number of bootstrap iterations")
	
	parser.add_option("--verbose","-v",action="store_true",help="Provide verbose feedback")
	parser.add_option("--debug",action="store_true",help="Provide debugging feedback")

	parser.usage="%prog [options] model.data observations.data "

	(opt, files)=parser.parse_args()

	if len(files) != 2 :
		parser.print_help()
		sys.exit(0)

### Configure the error message logging levels
	logging.basicConfig()
	logger=logging.getLogger()
	if opt.verbose:
		logger.setLevel(logging.INFO)
	if opt.debug:
		logger.setLevel(logging.DEBUG)

### Check that the input files are accessible.
	for file in files:
		if not os.access(file,os.R_OK):
			logger.error("Can't access the file: %s" % (file))
			sys.exit(-1)

### Open the file, split the lines on white space, load the 
### model in to data['model'] and the observations into data['obs']
### Assumes floating point numbers. 


### setup a data structure to read the data into, using a loop..

	order=['model','obs']
	filename={'model': files[0], 'obs': files[1]}
	col={'model': opt.col1-1, 'obs': opt.col2-1}
	data={'model': [], 'obs': []}
	for dataset in order:
		lines=open(filename[dataset])
		row=0
		for line in lines:
			row+=1
			cols=line.split()
			if (len(cols)-1 < col[dataset]):
				logger.error("Fewer than %d colums in file %s at row %d: skipping row" % (col[dataset],filename[dataset],row))
				continue
			if not re.match("^\d*(.(\d*))?$",cols[col[dataset]]):
				logger.error("skipping row %d of file %s: column %d  contains bad values" % (row,filename[dataset],col[dataset]))
				continue
			try:
				data[dataset].append(float(cols[col[dataset]].strip()))
			except:
				logger.error("Failed trying convert %s to float (skipping: file %s, row %d)" %(cols[col[dataset]],filename[dataset],row),exc_info=1)
				continue
		logger.info("Read %d rows from column %d of file %s" % (len(data[dataset]),col[dataset],filename[dataset]))

	print ad_bootstrap(data['obs'],data['model'],opt.nboots)



