# GenericSimAnnealing
I wanted an abstract implementation of simulated annealing to apply to a couple
of different optimization problems that I've been considering, and I thought
it was worth putting the implementation out there for anyone else to use. I'm 
such a template exists elsewhere, but I wanted the experience of implementing
it myself.

## Files
`generic_simulated_anneal.hpp` contains the abstract base class 
'AbstractOptimizationSolution' which the annealing function is written around,
and the function `genericSimulatedAnneal()` which performs the annealing.
This header is the only file needed to use the routine. Other directories are
example use cases (mostly problems I found interesting like DNA shotgun 
sequencing). 

## Usage
This section to be filled out as I fine-tune the code and decide on options to
include or not.

