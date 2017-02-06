# GenericSimAnnealing
I've been investigating several different problems that can be recast as 
optimization challenges; I wanted to try applying simulated annealing to these.
Since the point of my investigations was purely educational, it seemed 
worthwhile to do the implementation myself. I've done so in a manner that 
allows (hopefully) full flexibility in choice of cooling schedule, acceptance
probability function, and random number generation.

## Contents
The project consists of one file, `generic_simulated_anneal.hpp`. This file 
defines the following:

+ A struct GenSimAnnealParams, which holds the parameters for an annealing  
  run. This struct contains fields for, in order:
  
  * The maximum number of temperatures to be used (or the maximum number of  
    outer iterations, however you'd rather phrase it). 
  * The number of inner iterations to be performed at each temperature.
  * A cost reduction tolerance which will exit the annealing routine early if a  
    good enough value is found - this tests `curr_cost / initial_cost < tol` and  
    breaks if the inequality is true.
  * A boolean value indicating whether or not you'd like verbose output to  
    stdout. May be helpful in tuning your parameters.
    
+ The main simulated annealing routine. There are two versions provided. Both
  share the following template parameters:
  
  * A class `solution`. This represents a particular configuration in whatever
    the search space for your problem is. This class must implement two methods:
	`solution::cost()`, which returns an integral type representing the cost of
    this configuration, and `solution::perturb()`, which modifies the solution
    slightly in place. The value returned by the `cost()` function is that which
    the algorithm attempts to minimize.
  * A class `cost_t`. This is the type returned by the `solution` class' `cost()`
    method. It should be an integral type convertable to `double`.
  * A class `acceptance_t`; an object of this type shall be a function object
    implementing `double acceptance_t::operator()(cost_t, cost_t, double)`. 
    It will be called with a new cost, an old (lower) cost, and the current
    temperature; it should return a value between 0 and 1 giving the probability
    that the new solution is accepted.
  * A class `schedule_t`. This is a function object implementing
    `double schedule_t::operator()(unsigned)`, and will be called to determine 
    the temperature at each outer iteration of the method.
  
  The second version of the routine has a fifth parameter, a class 
  `generator_t`, which shall be a function object implementing a random number
  generator. When called with `operator()`, it should return a random number
  greater than 0, and it must implement a max() method that returns the maximum
  possible value of random number returned. 
 
  The default random number engine (first version) is `std::mt19937_64`, 
  initialized with current system time.
  
## Usage
The routine is called with one of the following signatures:

	genericSimulatedAnneal< ... > (const solution& s0,
	                               const GenSimAnnealParams& params,
	                               const acceptance_t& acceptance_f,
                                   const schedule_t& schedule);

	genericSimulatedAnneal< ... > (const solution& s0,
	                               const GenSimAnnealParams& params,
	                               const acceptance_t& acceptance_f,
	                               const schedule_t& schedule,
	                               generator_t& generator);

For best performance, make sure that calls to `solution::cost()` are as cheap
as possible, and ensure that data in the solution class that remains constant
when the solution is perturbed is not deep copied by the assignment operator.
There is no way around making one deep copy of the solution at each iteration,
so having e.g. a `std::vector` as a data member if it isn't changed by
perturbations would cause a lot of unnecessary overhead in copies.

The subdirectory "TravelingSalesman" has a sample implementation of solving 
the traveling salesman problem. 

### License
The code in `generic_simulated_anneal.hpp` is Copyright (c) Sean McBane, 2017
under the terms of the MIT License (see LICENSE for detailed legal statement). 
Example code in subdirectories is not copyrighted.

Sean McBane, 2/6/2017 <seanmcb94@gmail.com>
