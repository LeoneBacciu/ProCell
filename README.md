# ProCell

ProCell is an [award winning](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8791485) modeling and simulation framework designed to investigate cell proliferation dynamics that, differently from other approaches, takes into account the inherent stochasticity of cell division events.

ProCell manipulates raw data coming from flow cytometry experiments. Specifically, uses as input:
- a histogram of initial cell fluorescences (e.g., GFP signal in the population); 
- the number of different sub-populations, along with their proportions;
- the mean and standard deviation of division time for each population;
- a fluorescence minimum threshold; 
- a maximum simulation time T, expressed in hours.

The output produced by ProCell is a histogram of GFP fluorescence after time T.

## Installing and using ProCell

This part is still under development. I am setting up the GITHUB repository and currently creating a PyPI package to simplify the installation process. 

Meanwhile, you can just download the source code and launch ProCell's GUI (python gui/gui.py), [designed by Luca Zanini](https://www.behance.net/gallery/84353179/-UX-UI-PROCELL).

![ProCell's GUI](https://lh5.googleusercontent.com/wj8YpG-YqbE2wkeN2kJXKiKVox_Kj3GL5kt3X9dKRRuxqqw6fTgbUNNakCGrNspP-zeGJ1fqH6sFrgQOPC94IB7jRwUZ8syAlFst85Mji6QV4NMKaZDkfbhEtpYXKNRctWt4IaTc)

A tutorial about modeling, calibration and simulation will be published soon. You can find a list of ProCell's pre-requisites here:

https://docs.google.com/document/d/1OC0DDQQHAKs6sMOi7hWleWcwgU04h9RAsZzOye6dnAU/edit?usp=sharing

## GPU version of ProCell

We are developing a GPU-accelerated version of ProCell, able to strongly reduce the computational effort of simulations. The alpha implementation can be downloaded from: https://github.com/ericniso/cuda-pro-cell

## More info about ProCell

If you need additional information about ProCell please write to: nobile@disco.unimib.it.

## Citing ProCell

Nobile M.S., Vlachou T., Spolaor S., Bossi D., Cazzaniga P., Lanfrancone L., Mauri G., Pelicci P.G., Besozzi D.: Modeling cell proliferation in human acute myeloid leukemia xenografts, Bioinformatics, 2019 (in press)

Nobile M.S., Vlachou T., Spolaor S., Cazzaniga P., Mauri G., Pelicci P.G., Besozzi D.: ProCell: Investigating cell proliferation with Swarm Intelligence, Proceedings of the 16th IEEE International Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB 2019), Certosa di Pontignano, Siena, Tuscany, Italy, 2019 
