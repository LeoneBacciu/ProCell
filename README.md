# ProCell

ProCell is a modeling and simulation framework to investigate cell proliferation dynamics that, differently from other approaches, takes into account the inherent stochasticity of cell division events.

ProCell uses as input:
- a histogram of initial cell fluorescences (e.g., GFP signal in the population); 
- the number of different sub-populations, along with their proportions;
- the mean and standard deviation of division time for each population;
- a fluorescence minimum threshold; 
- a maximum simulation time T, expressed in hours.

The output produced by ProCell is a histogram of GFP fluorescence after time T.

## Installing and using ProCell

This part is still under development. I am setting up the GITHUB repository and currently creating a PyPI package to simplify the installation process. 

Meanwhile, you can just download the source code and launch ProCell's GUI (python gui/gui.py). A tutorial about modeling, calibration and simulation will be published soon.

## More info about ProCell

If you need additional information about ProCell please write to: nobile@disco.unimib.it.
