from .cells_population import *
from scipy import stats


def rebin(series, lower, upper, N=1000):
	bins = logspace(ceil(log10(lower)), ceil(log10(upper)), N)
	res  = zeros(N, dtype=int)
	pos  = 0
	for k,v in zip(series.T[0], series.T[1]):
		while(k>bins[pos]):
			pos+=1
			if pos>=N:
				print("ERROR: bins exceeded, this should never happen, aborting.")
				print(lower,upper,N)
				exit()
		res[pos] += v
	return res, bins
	

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

	fitness = sum(abs(res1-res2))
	return fitness, stats.ks_2samp(res1, res2)


if __name__ == '__main__':	

	pass