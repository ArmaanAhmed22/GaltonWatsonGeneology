# Galton Watson Geneology Simulator

Can simulate continuous-time Galton Watson processes whilst tracking the geneology associated to the population. The algorithm has the same time complexity (at the cost of increased space complexity on the order of the size of the Galton Watson tree).

## Features
- Simulate multi-type continuous-time Galton Watson processes (for offspring distributions with finite support)
- Can use the geneology to compute the site-frequency-spectrum and single-cell mutational burden
- Have various stopping conditions at which to save data (at times, population sizes, sizes of certain types of the population, Galton Watson ancestral tree sizes, etc).
- Can obtain pruned Galton-Watson trees as well as the reconstructed Galton-Watson tree from taking a sample.