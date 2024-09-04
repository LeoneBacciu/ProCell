from __future__ import print_function
import sys, os
from fstpso import FuzzyPSO
from .procell_core import Simulator
from .fitness import rebin

from pylab import  *
from scipy.linalg import norm
import getopt
import datetime, time
from subprocess import check_output
from collections import defaultdict


_SQRT2 = np.sqrt(2)

def dummy(): return 0

def _prepare_files_for_GPU(names, prop, mean, st):
	# each cell population = one row
	# proportion [space] mean [space] std
	with open("__temporary__", "w") as fo:
		for name in names:
			row = "%s %s %s\n" % (prop[name], mean[name], st[name])
			fo.write(row)
	

def _launch_GPU_simulation(executable, initial_histo, model_file, time_max, PHI, names):

	"""
	print "I RECEIVED"
	print executable
	print initial_histo
	print model_file
	print time_max
	print PHI
	print names
	"""

	command = [executable, "-h", initial_histo, "-c", model_file, "-t", str(time_max), "-p", str(PHI), "-r"]
	#print (" ".join(command))
	ret = check_output(command)
	ret = ret.decode("utf-8")
	split_rows = ret.split('\n')
	result = defaultdict(dummy)
	types  = defaultdict(list)


	for row in split_rows:
		row = row.strip("\r")
		#print (row)
		try:
			tokenized_row = list(map(float, row.split("\t")))
		except ValueError:
			continue

		fluorescence = tokenized_row[0]
		total = tokenized_row[1]
		result[fluorescence]=total

		for n, amount in enumerate(tokenized_row[2:]):
			types[fluorescence].append([names[n]]*int(amount))

	return result, types


class LogFSTPSO(FuzzyPSO):

	"""
	def UpdateLocalBest(self, verbose=False, semiverbose=True):
		with open(self.dump_fitness, "a") as fo: fo.write(str(self.G.CalculatedFitness)+"\n")
		with open(self.dump_struct, "a") as fo: fo.write("\t".join(map(str, self.G.X))+"\n")
		super(LogFSTPSO, self).UpdateLocalBest(verbose=verbose, semiverbose=semiverbose)
	"""

def hellinger1(p, q):
    return norm(np.sqrt(p) - np.sqrt(q)) / _SQRT2


def fitness_evaluate(CP, TGT, N=1000, draw=False):
	ratio = sum(CP.frequencies)/sum(TGT.T[1])
	TGT.T[1] *= ratio
	
	# lower and upper fluorescence values
	upper1 = max(CP.fluorescences)
	lower1 = min(CP.fluorescences)
	upper2 = max(TGT.T[0])
	lower2 = min(TGT.T[0])
	lower = min(lower1,lower2)
	upper = max(upper1,upper2)

	bins =  logspace(ceil(log10(lower)), ceil(log10(upper)), N)
	res1 = zeros(N, dtype=int)
	res2 = zeros(N, dtype=int)

	pos = 0
	for k,v in zip(CP.fluorescences, CP.frequencies):
		while(k>bins[pos]):
			pos+=1
		res1[pos] += v
	res2, bins = rebin(TGT, lower, upper, N)

	if draw:
		from matplotlib.pyplot import bar, legend, xscale, xlim, xlabel, ylabel, show		
		bar(bins, res2, bins/N*5, linewidth=0, label="Scaled simulation", color="green", alpha=.7)
		bar(bins, res1, bins/N*5,  linewidth=0, label="Experimental data (10 days)", color="red", alpha=.7)
		legend(loc="best")
		xscale("symlog")
		xlim(bins[0],bins[-1])
		xlabel("Fluorescence")
		ylabel("Frequency")
		show()

	#fitness = sum(abs(res1-res2))
	fitness = hellinger1(res1, res2)
	return fitness, stats.ks_2samp(res1, res2)


def fitness_gui(p, arguments, return_dictionaries=False):

	gui = arguments['form']
	names = gui._population_names
	N = len(names)

	use_gpu = (gui._path_to_GPU_procell is not None)

	def rel2abs_prob(V):
		if len(V)==0: V = [1] # for single-population cases

		putative_proportions = [V[0]]
		consumed = V[0]
		for p in V[1:N-1]:
			space_left = 1.-consumed
			#print "Consumed: %f\tLeft: %f" % (consumed, space_left)
			new_slice = space_left*p
			#print "New slice: %f" % new_slice
			putative_proportions.append( new_slice )
			consumed += new_slice
		putative_proportions.append(1.-consumed)
		return putative_proportions

	putative_proportions = rel2abs_prob(p[0:N-1])
	proportions = dict(zip(names, putative_proportions))

	putative_means = []
	pos = N-1
	for v in gui._population_means:
		if v!="-":
			putative_means.append(p[pos])
			pos+=1
		else:
			if use_gpu:
				putative_means.append(-1) 					# used by GPU version
			else:
				putative_means.append(sys.float_info.max)	# used by CPU version

	mean_div      = dict(zip(names, putative_means))

	putative_std = []
	for v in gui._population_std:
		if v!="-":
			putative_std.append(p[pos])
			pos+=1
		else:
			if use_gpu:
				putative_std.append(-1)						# used by GPU version
			else:
				putative_std.append(0)						# used by CPU version

	std_div = dict(zip(names, putative_std))	

	print ("Testing the following parameterization:")
	print (" - Proportions:", proportions)
	print (" - Mean:", mean_div)
	print (" - STD:", std_div)

	phi = float(gui.fluorescencethreshold.value())
	time_max = float(gui.simulationtime.value() ) 
	lower = float(gui.lowerbin.value())
	higher = float(gui.higherbin.value())
	calcbins = int(gui.bins.value())
	PHI      = float(gui.fluorescencethreshold.value())

	if gui._path_to_GPU_procell is not None:

		_prepare_files_for_GPU(names, proportions, mean_div, std_div)

		try:
			result_simulation, result_simulation_types = _launch_GPU_simulation(
				gui._path_to_GPU_procell, 
				gui._initial_histo_path, 
				"__temporary__", 
				time_max, 
				PHI,
				names
				)
		except:
			print ("ERROR: simulation on GPU crashed for unknown reasons.")

		#print "done"

	else: 

		Sim = Simulator()

		result_simulation, result_simulation_types = Sim.simulate(
				path=gui._initial_histo_path, 
				types=names, 
				proportions=proportions,
				div_mean=mean_div,
				div_std=std_div, 
				time_max=time_max, 
				verbose=False, 
				phi=phi, 
				distribution="gauss"
			)

	# rebinning
	sorted_res = array(sorted(result_simulation.items()))
	result_simulation_rebinned, bins_result_simulation_rebinned = rebin(sorted_res, lower, higher, N=calcbins)
	target_rebinned, bins_target_rebinned   	= rebin(gui._target_histo, lower, higher, N=calcbins)
	ratio = 1.0*sum(result_simulation_rebinned)/sum(target_rebinned)
	fitness = hellinger1(result_simulation_rebinned/ratio, target_rebinned)

	print (" - Fitness:", fitness, "\n")

	if return_dictionaries==False: 	return fitness
	return proportions, mean_div, std_div



"""
def fitness(p):

	global C

	if   C._modelname == "model1":
		# two parameters (mu_s, sigma_s)
		p_s = 1.
		mu_s = p[0]
		sigma_s = p[1]
		proportions = {"s": p_s}
		mean_div = {"s": mu_s}
		std_div = {"s": sigma_s}

	elif C._modelname == "model2":
		# three parameters (p_q, mu_s, sigma_s)
		p_q = p[0]
		p_s = 1.-p_q
		mu_s = p[1]
		sigma_s = p[2]
		proportions = {"q": p_q, "s": p_s}
		mean_div = {"q": sys.float_info.max, "s": mu_s}
		std_div = {"q": 0, "s": sigma_s}

	elif C._modelname == "model3":
		# six parameters (p_q, p_s, mu_f, mu_s, sigma_f, sigma_s)
		p_q = p[0]
		p_s = (1.-p_q)*p[1]
		p_f = 1.-p_q-p_s
		mu_f = p[2]
		mu_s = p[3]
		sigma_f = p[4]
		sigma_s = p[5]
		proportions = {"q": p_q, "s": p_s, "f": p_f}
		mean_div = {"q": sys.float_info.max, "s": mu_s, "f": mu_f }
		std_div = {"q": 0, "s": sigma_s, "f": sigma_f }

	elif C._modelname == "model4":
		# five parameters (p_s, mu_f, mu_s, sigma_f, sigma_s)
		p_s = p[0]
		p_f = 1.-p_s
		mu_f = p[1]
		mu_s = p[2]
		sigma_f = p[3]
		sigma_s = p[4]
		proportions = {"s": p_s, "f": p_f}
		mean_div = {"s": mu_s, "f": mu_f }
		std_div = {"s": sigma_s, "f": sigma_f }

	elif C._modelname == "model5":
		# nine parameters (p_q, p_s, p_n, mu_f, mu_n, mu_s, sigma_f, sigma_n, sigma_s)
		p_q = p[0]
		p_s = (1.-p_q)*p[1]
		p_n = (1.-p_s-p_q)*p[2]
		p_f = (1.-p_s-p_q-p_n)
		mu_f = p[3]
		mu_n = p[4]
		mu_s = p[5]
		sigma_f = p[6]
		sigma_n = p[7]
		sigma_s = p[8]
		proportions = {"q": p_q, "s": p_s, "n": p_n, "f": p_f}
		mean_div = {"q": sys.float_info.max, "s": mu_s, "n": mu_n, "f": mu_f }
		std_div = {"q": 0, "s": sigma_s, "n": sigma_n, "f": sigma_f }
	
	Sim = Simulator()

	print ("Proportions:", proportions)
	print ("Mean:", mean_div)
	print ("STD:", std_div)

	time_max = float(gui.simulationtime.value() ) 
	lower = float(gui.lowerbin.value())
	higher = float(gui.higherbin.value())
	calcbins = int(gui.bins.value())
	PHI      = float(gui.fluorescencethreshold.value())


	result_simulation, result_simulation_types = Sim.simulate(
			path=C._H0, types=C._types, proportions=proportions,
			div_mean=mean_div, div_std=std_div, time_max=time_max, verbose=False, phi=PHI, distribution=C._distribution)

	sorted_res = array(sorted([ [a,b] for (a,b) in zip(result_simulation.keys(), result_simulation.values())]))
	result_simulation_rebinned, bins_result_simulation_rebinned = rebin(sorted_res, lower, higher, N=calcbins)
	target_rebinned, bins_target_rebinned   	= rebin(gui._target_histo, lower, higher, N=calcbins)

	ratio = 1.0*sum(result_simulation_rebinned)/sum(target_rebinned)

	#plot(result_simulation_rebinned/ratio)
	#plot(target_rebinned)
	#print sum(result_simulation_rebinned/ratio)
	#print sum(target_rebinned)
	#show()
	fitness = hellinger1(result_simulation_rebinned/ratio, target_rebinned)
	#fitness = fitness_evaluate(result_simulation_rebinned, target_rebinned)

	print ("Fitness:", fitness, "\n")

	return fitness

"""

class Calibrator(object):

	def __init__(self):
		self._Htarget = None
		self._H0 = None
		self._types = []
		self._modelname = None
		self._time_max = 0
		self.TGT = None
		self._distribution = "gauss"
		self._abort_variable = False

	def set_types(self, types):
		self._types = types

	def set_target_from_file(self, path):
		self._Htarget = path
		self.TGT = loadtxt(self._Htarget)

	def set_initial_histogram_from_file(self, path):
		self._H0 = path

	def set_model_name(self, modelname):
		self._modelname = modelname

	def set_output_dir(self, path):
		try:
			os.mkdir(path)
		except:
			pass
		finally:
			self._outputdir = path

	def set_time_max(self, time):
		self._time_max = time

	def init_dump(self, p1, p2):
		with open(p1, "w") as fo: fo.write("# dump fitness\n")
		with open(p2, "w") as fo: fo.write("# dump structure\n")

	def time_stamp(self):
		return datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d-%H%M')

	def calibrate(self, repetition=0, max_iter=10):

		# controls
		if self._types == []: 		raise Exception("Types of cells not specified")
		if self._H0 == None: 		raise Exception("Please specify the initial experimental histogram")
		if self._Htarget == None: 	raise Exception("Please specify the target experimental histogram")
		if self._time_max<=0: 		raise Exception("Please specify a maximum time for estimation")
		if self._distribution=="": 	raise Exception("Please specify the type of distribution")

		print (" * Initial controls passed, loading FST-PSO")
		print ("   using distribution: %s" % self._distribution)

		FP = LogFSTPSO()
		FP.set_swarm_size(50)
		FP.dump_fitness = self._outputdir+os.sep+"%s_rep%d_%s_fitness_%s.txt" % (self._modelname, repetition, self._distribution, str(self.time_stamp())) 
		FP.dump_struct = self._outputdir+os.sep+"%s_rep%d_%s_structure_%s.txt" % (self._modelname, repetition, self._distribution, str(self.time_stamp())) 

		self.init_dump(FP.dump_fitness, FP.dump_struct)

		#FP.enabled_settings = ["inertia", "cognitive", "social", "maxvelocity"]
		FP.disable_fuzzyrule_minvelocity()
	
		if self._modelname=="model1": 
			if self._distribution == "gauss":
				FP.set_search_space([
					[30., 504.], [1e-10, 30.]
				])
			elif self._distribution == "gamma":
				FP.set_search_space([
					[6., 30.], [5., 15.]
				])

		elif self._modelname=="model2":
			if self._distribution == "gauss":
				FP.set_search_space([
					[1e-10, 1.], [30., 504.], [1e-10, 30.]
				])			
			elif self._distribution == "gamma": 
				FP.set_search_space([
					[1e-10, 1.], [6., 30.], [5., 15.]
				])

		elif self._modelname=="model3":
			if self._distribution == "gauss":
				FP.set_search_space([
					[1e-10, 1.], [1e-10, 1.], 
					[30., 63.], [63., 504.], 
					[1e-10, 20.], [1e-10, 30.]
				])
			elif self._distribution == "gamma":
				FP.set_search_space([
					[1e-10, 1.], [1e-10, 1.], 
					[6., 9.], [9., 30.], 
					[5., 15.], [5,20]					
				])

		elif self._modelname=="model4":
			if self._distribution == "gauss":
				FP.set_search_space([
					[1e-10, 1.], 
					[30., 63.], [63., 504.], 
					[1e-10, 20.], [1e-10, 30.]
				])			
			elif self._distribution == "gamma":
				FP.set_search_space([
					[1e-10, 1.], 
					[6., 9.], [9., 30.], 
					[5., 15.], [5,20]					
				])

		elif self._modelname=="model5":
			if self._distribution == "gauss":
				FP.set_search_space([
					[1e-10, 1.], [1e-10, 1.], [1e-10, 1.],
					[30., 63.], [63., 200.], [200., 504.], 
					[1e-10, 20.], [1e-10, 25.], [1e-10, 30.]
				])
			elif self._distribution == "gamma":
				FP.set_search_space([
					[1e-10, 1.], [1e-10, 1.], [1e-10, 1.], 
					[6., 9.], [9., 15.], [15., 30.], 
					[5., 15.], [5., 20.], [5., 25.]					
				])
		else:
			raise Exception("Problem: model not supported")


		FP.set_fitness(fitness, skip_test=True)
		
		return FP.solve_with_fstpso(creation_method={'name': 'uniform'}, 			
			max_iter=max_iter, 
			dump_best_fitness=FP.dump_fitness, 
			dump_best_solution=FP.dump_struct)

	def calibrate_gui(self, repetition=0, max_iter=10, search_space=None, 
		swarm_size=50, form=None, append_time_stamp=False, loginit=False):

		# controls
		if self._types == []: 		raise Exception("Types of cells not specified")
		if self._H0 == None: 		raise Exception("Please specify the initial experimental histogram")
		if self._Htarget == None: 	raise Exception("Please specify the target experimental histogram")
		if self._time_max<=0: 		raise Exception("Please specify a maximum time for estimation")
		if self._distribution=="": 	raise Exception("Please specify the type of distribution")

		print (" * Initial controls passed, loading FST-PSO")
		print ("   using distribution: %s" % self._distribution)

		if append_time_stamp:
			time_stamp = str(self.time_stamp())
		else:
			time_stamp = ""

		FP = LogFSTPSO()
		FP.dump_fitness = self._outputdir+os.sep+"%s_rep%d_%s_fitness_%s.txt" % (self._modelname, repetition, self._distribution, time_stamp) 
		FP.dump_struct = self._outputdir+os.sep+"%s_rep%d_%s_structure_%s.txt" % (self._modelname, repetition, self._distribution, time_stamp) 
		FP.set_search_space(search_space)
		FP.set_fitness(fitness_gui, skip_test=True, arguments = {'form': form})
		FP.set_swarm_size(swarm_size)

		FP.disable_fuzzyrule_minvelocity()

		if loginit:
			logtype = 'logarithmic'
		else:
			logtype = 'uniform'
		
		return FP.solve_with_fstpso(creation_method={'name': logtype}, 			
			max_iter=max_iter, 
			dump_best_fitness=FP.dump_fitness, 
			dump_best_solution=FP.dump_struct,
			#initial_guess_list = [[0.5, 25,30, 1,1]]
		)


if __name__ == "__main__":

	argv = sys.argv

	distribution = "gauss"

	initial_histogram = "./target_GFPpos/AML9_Exp3_t10d_UT1_BM_GFPpos.txt"
	target_histogram  = "./target_GFPpos/AML9_Exp3_t10d_DOX3_BM_GFPpos.txt"

	try:
		opts, args = getopt.getopt(argv[1:],"i:t:m:o:r:d:",["initial=", "target=", "model=", "output=", "repetition=", "distribution="])
	except getopt.GetoptError:
		print ('test.py -i <initial histogram> -t <target histogram> -m <model> -o <outputdir> -r <repetition> -d <distribution>')
		sys.exit(2)
	
	for opt, arg in opts:
		if opt == '-h':
			print ('test.py -i <initial> -t <target> -o <outputdir> -r <repetition> -d <distribution>')
			sys.exit()
		elif opt in ("-i", "--initial"):
			initial_histogram = arg
			print (" * Using initial histogram", initial_histogram)
		elif opt in ("-t", "--target"):
			target_histogram = arg
			print (" * Using target histogram", target_histogram)
		elif opt in ("-m", "--model"):
			model = arg
			print (" * Using model", model)
		elif opt in ("-o", "--output"):
			outputdir = arg
		elif opt in ("-r", "--repetition"):
			repetition = int(arg)
		elif opt in ("-d", "--distribution"):
			distribution = arg

	if   model == "model1":
		types = ["s"]
	elif model == "model2":
		types = ["q", "s"]
	elif model == "model3":
		types = ["q", "f", "s"]
	elif model == "model4":
		types = ["f", "s"]
	elif model == "model5":
		types = ["q", "f", "n", "s"]


	C = Calibrator()
	C._distribution = distribution

	C.set_model_name(model)
	C.set_output_dir(outputdir)
	C.set_types(types)
	C.set_time_max(10*24) #hours
	C.set_initial_histogram_from_file(initial_histogram)
	C.set_target_from_file(target_histogram)

	res = C.calibrate(repetition=repetition, max_iter=100)
	print (res)


