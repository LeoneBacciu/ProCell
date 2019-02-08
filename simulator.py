from cell import *
from cells_population import *
from numpy import *
from scipy.stats import entropy
from fitness import rebin
import gc
import collections


class CellsPopulationSimulator(object):

	def __init__(self, population=None, fluorescence_threshold=0):
		
		if population==None:
			print "ERROR: please specify a population"
			exit(3)

		self.fluorescence_threshold=fluorescence_threshold
		self.population = population
		self.final = {}


	def plot_output(self, ref=None, filename=None, file_distances=None,
	 plot_title="", time_max=24*10, save_precalc=None, show_initial=False, args=None):		
		
		from pylab import *
		import matplotlib 
		matplotlib.rcParams.update({'font.size': 18})

		fig = figure(figsize=(16,9), dpi=100)
		fig2 = figure(figsize=(16,9), dpi=100)
		axr = fig.add_subplot(111)
		axl = fig2.add_subplot(111)
		
		N = 1000

		max_frequency = self.population.fluorescences[-1]
		upper = max_frequency
		lower = 1e-20

		ratio = 1

		if ref!=None:

			t=0
			A = loadtxt(ref)
			ratio = sum(self.final.values())/sum(A.T[1])
			A, bins = rebin(A, lower, upper, N=N)

			# plot target
			for k,v in zip(bins, A):
				if t==0:
					axr.plot([k,k], [0,v], label="Target", alpha=0.5, color="green", linewidth=4)
					t=t+1
				else:
					axr.plot([k,k], [0,v], alpha=0.5, color="green", linewidth=4)

		# plot output simulation
		sorted_output = array(self.get_ordered_output().items())
		B, bins= rebin(sorted_output, lower, upper, N=N)

		if save_precalc!=None:
			self.save_precalc(save_precalc,A,B/ratio)
		
		KLD = self.calculate_distance(A, B/ratio)
		HD  = self.hellinger_distance(A, B/ratio)
		with open(file_distances, "w") as fo:		
			fo.write("# KLD\tHellinger\n")
			fo.write(str(KLD)+"\t"+str(HD))

		
		t = 0
		for k,v in zip (bins, B):
			if t==0:
				axr.plot([k,k],[0,v/ratio], c="r", label="Predicted distribution", alpha=0.5, linewidth=4)
				t+=1
			else:
				axr.plot([k,k],[0,v/ratio], c="r", alpha=0.5, linewidth=4)


		# initial population		
		if show_initial == True:
			t = 0
			for k,v in zip(self.population.fluorescences, self.population.frequencies):
				if t==0:
					axr.plot([k,k],[0,v], c="black", label="Initial histogram", alpha=0.3)
					t+=1
				else:
					axr.plot([k,k],[0,v], c="black", alpha=0.3)
		else:
			t = 0
			for k,v in zip(self.population.fluorescences, self.population.frequencies):
				if t==0:
					axl.plot([k,k],[0,v], c="black", label="Initial histogram", alpha=0.3)
					t+=1
				else:
					axl.plot([k,k],[0,v], c="black", alpha=0.3)


		axl.set_xscale("symlog")
		axr.set_xscale("symlog")
		axl.set_xlim(0, 10**log10(max_frequency))
		axr.set_xlim(0, 10**log10(max_frequency))
		axl.set_xlabel("Fluorescence GFP-A")
		axr.set_xlabel("Fluorescence GFP-A")
		axl.set_ylabel("Cells frequency")
		axr.set_ylabel("Cells frequency")
		axl.set_title(organ+", $t=0$ hours")
		ttl = axr.set_title( plot_title+ ", $t_{max}="+str(time_max)+"$ hours")		
		ttl.set_position([.5, 1.05])
		axl.legend(loc="best")
		axr.legend(loc="best")
		tight_layout()
		
		if filename==None:
			show()
		else:
			print " * saving", filename
			fig.savefig(filename)
		fig.clf()
		gc.collect()


	def save_precalc(self, output, A,B):
		savetxt(output+"A", A)
		savetxt(output+"B", B)

	def load_precalc(self, input):
		A = loadtxt(input+"A")
		B = loadtxt(input+"B")
		return A,B



	def calculate_distance(self, A, B):
		A = array(map(int,A))
		A = [1e-10 if x==0 else x for x in A]
		B = array(map(int,B))
		B = [1e-10 if x==0 else x for x in B]
		return entropy( A,qk=B )


	def hellinger_distance(self, A, B):
		B /= sum(B)
		A = A.astype(float)/sum(A)
		sim = 0
		sqrt2 = 1./sqrt(2.)
		tot = 0
		for v1, v2 in zip(A,B):
			tot += (sqrt(v1)-sqrt(v2))**2
		sim = sqrt2*tot
		return sim
		

	def get_ordered_output(self):
		od = collections.OrderedDict(sorted(self.final.items()))
		return od

	def save_output(self, filename=None):
		import time, datetime
		ts = time.time()
		ts = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H%M%S')
		od = self.get_ordered_output()
		if filename==None:
			filename='result'+ts+'.csv'
		with open(filename, 'w') as f:
			for k, v in od.iteritems():
				f.write(str(k)+"\t"+str(v)+"\n")


	def simulate(self, time_max=0, static=True, limit=0):

		while(len(self.population.stack)>0):

			cell_to_simulate = self.population.stack.pop()			

			# is it quiescient? then go directly to the end of the simulation
			if cell_to_simulate.timer == -1:				
				try:
					self.final[cell_to_simulate.fluorescence]+=1
				except:
					self.final[cell_to_simulate.fluorescence]=1

			else:

				# do I still have margin for a cell division?
				if cell_to_simulate.T + cell_to_simulate.timer < time_max:

					# take one step forward and divide
					cell_to_simulate.T += cell_to_simulate.timer					
					newfluorescence = cell_to_simulate.fluorescence*0.5

					# create new cells only if the fluorescence is high enough, else continue
					if newfluorescence>self.fluorescence_threshold:
						if static:
							nc1 = self.population.create_cell(fluorescence=newfluorescence, cell_type=cell_to_simulate.type, ID=self.population.last_ID+1)
							nc2 = self.population.create_cell(fluorescence=newfluorescence, cell_type=cell_to_simulate.type, ID=self.population.last_ID+2)							
							nc1.T = cell_to_simulate.T
							nc2.T = cell_to_simulate.T
							self.population.last_ID += 2
						self.population.stack.append(nc1)
						self.population.stack.append(nc2)
					else:
						del cell_to_simulate

				# simulation over for this cell: count and terminate
				else:
					try:
						self.final[cell_to_simulate.fluorescence]+=1
					except:
						self.final[cell_to_simulate.fluorescence]=1


class SimulationWrapper(object):

	def __init__(self, args):
		
		self.CP = CellsPopulation(
			initial_histogram=args['initial_histogram'],
			timers=args['timers'],	
			types_p=args['types_p'],
			sigma=args['sigma']
		)

		self.SIM = CellsPopulationSimulator(
			population=self.CP, 
			fluorescence_threshold=args['fluorescence_threshold']
		)	 

	def simulate(self, TIME_MAX):
		self.SIM.simulate(time_max=TIME_MAX)
		return self.SIM.get_ordered_output()


if __name__ == '__main__':

	pass

