# GenericSimAnnealing
I wanted to try using simulated annealing for a couple of different optimization
problems that I'd looked at lately. Since I've been working on these problems
for my own edification anyway, it seemed most profitable to roll my own version,
and since in my limited experience I've already found several different problems
that the algorithm might be applied to, a generic implementation seemed best. 
I thought surely someone else could use the code (though I'm also sure that 
other similar codes exist), and that's how this repository came to be.

## Contents
The project consists of one file, `generic_simulated_anneal.hpp`. This file 
defines the following:

+ A struct GenSimAnnealResults, which holds the results from an annealing run.  
  It has four fields, for the solution, the final value of the cost function,  
  the number of objective function evaluations, and the number of (outer)  
  iterations performed, respectively. This struct is templated on the solution  
  class and the cost function return type.
+ A struct GenSimAnnealParams, which holds the parameters for an annealing  
  run. This struct contains fields for, in order:
  
  * The maximum number of temperatures to be used (or the maximum number of  
    outer iterations, however you'd rather phrase it). 
  * The number of inner iterations to be performed at each temperature.
  * The temperature reduction factor alpha by which the current temperature is  
    multiplied at each outer iteration.
  * A cost reduction tolerance which will exit the annealing routine early if a  
    good enough value is found - this tests `curr_cost / initial_cost < tol` and  
    breaks if the inequality is true.
  * A boolean value indicating whether or not you'd like verbose output to  
    stdout. May be helpful in tuning your parameters.
    
+ A set of default GenSimAnnealParams, and
+ The main simulated annealing routine; genericSimulatedAnneal. The routine is  
  templated on a solution class, a cost function return type, and a user defined
  function returning the acceptance probability. More on that under 'Usage'.
      
  
## Usage
To use the annealing routine the user must define a class that represents 
somehow the solution to the optimization problem that you wish to solve. How 
this class is implemented does not matter to the algorithm, except for two 
required methods:

+ The class must supply a method cost() that evaluates the cost of the current  
  solution. The algorithm attempts to modify the solution to minimize the value  
  returned by this function.
+ The class must supply a method perturb() that makes a small, random change to  
  the solution *in place*. It could make a large change, or no change... but  
  that would defeat the purpose.
  
While these are technically the only requirements for the solution class, the 
user should note that a copy is made at every inner iteration before the 
solution is modified, so that if the change is not accepted the old solution is
still available. If your class's assignment operator does a deep copy of 
anything (such as an STL container), you will incur a large overhead for all the
copies/allocations. If such a member does NOT need to be modified, use a pointer
to the data in your class instead. See the traveling salesman example code for
this in practice. Secondly, it is important that the cost function be relatively
inexpensive - but that should be common sense.

You must also provide the third template argument, which is a function to 
compute the acceptance probability from two costs and the current temperature.
The function should expect the second cost to be greater than the first and 
return the probability that a solution with the second cost is chosen over the 
first. An important note first, however: *the temperature is always between 0 and 1*.
This is different from many other simulated annealing implementations, where the
temperature is often a function of the iteration (maybe user-defined) and the
probability function is always the same. In my implementation, T starts at 1.0 
and is multiplied by alpha at each iteration. The acceptance probability function
that you supply should take this into account.

For detailed examples of using the code, please see the subdirectories in this
repository and the comments in the source code.

### License
This code is copyright (c) Sean McBane under the terms of the MIT License.

Sean McBane, 3/30/2016 <seanmcb94@gmail.com>
