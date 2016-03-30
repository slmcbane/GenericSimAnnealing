/* Example use case for GenericSimAnneal, solving the classic NP-complete
 * problem of the traveling salesman who is trying to visit a set of cities
 * in a closed loop with the least possible distance covered. Please note
 * that the acceptance probability function here should not be taken as a
 * reference as to what makes a good one - it's just a guess at something
 * about right. It works decently.
 *
 * Note the use of a shared_ptr<vector<long> > to hold the coordinates -
 * this avoids the deep copy of the vectors that would be performed at each
 * inner iteration otherwise.
 *
 * This code is released as an example under no copyright whatsoever.
 */

#include "generic_simulated_anneal.hpp"
#include <ctime>
#include <vector>
#include <cassert>
#include <cmath>
#include <memory>

// A class representing a particular tour.
// xcoords and ycoords hold the x- and y- coordinates of all the cities.
// visited is a vector of the order in which cities are visited (beginning
// and ending with city # 0).
class Tour {
public:
  size_t numcities;
  std::shared_ptr<std::vector<long> > xcoords, ycoords;
  std::vector<size_t> visited;

  Tour(std::vector<long> xs, std::vector<long> ys);
  unsigned long cost() const;
  void perturb();
};

// Constructor.
Tour::Tour(std::vector<long> xs, std::vector<long> ys) {
  assert(xs.size() == ys.size());
  numcities = xs.size();

  xcoords = std::make_shared < std::vector<long>
      > (std::vector<long>(numcities));
  ycoords = std::make_shared < std::vector<long>
      > (std::vector<long>(numcities));

  visited = std::vector<size_t>(numcities + 1);
  for (size_t i = 0; i < numcities; ++i) {
    (*xcoords)[i] = xs[i];
    (*ycoords)[i] = ys[i];
    visited[i] = i;
  }
  visited[numcities] = 0;
}

// Computes the cost of the tour as given by the floor of the
// Pythagorean distance.
unsigned long Tour::cost() const {
  unsigned long dist = 0;
  std::vector<long>& x = *xcoords;
  std::vector<long>& y = *ycoords;
  for (size_t i = 0; i < numcities; ++i) {
    long dx = x[visited[i + 1]] - x[visited[i]];
    long dy = y[visited[i + 1]] - y[visited[i]];
    double d = sqrt(double(dx * dx + dy * dy));
    dist += (unsigned long) floor(d);
  }
  return dist;
}

// Swaps two random cities on the tour; this is a good enough perturbation
// function for an example case.
void Tour::perturb() {
  size_t city1 = 1
      + (unsigned long) ceil((double(rand())) / RAND_MAX * (numcities - 3));
  size_t city2 = city1;
  while (city2 == city1)
    city2 = 1
        + (unsigned long) ceil((double(rand())) / RAND_MAX * (numcities - 2));

  size_t temp = visited[city1];
  visited[city1] = visited[city2];
  visited[city2] = temp;
}

// p = exp((c1-c2)/T) is a common formulation; I introduced a constant factor
// to scale T since it is between 0 and 1 in the algorithm.
double acceptanceProbability(unsigned long c1, unsigned long c2, double T) {
  double a = 600.0;
  double cost1 = double(c1);
  double cost2 = double(c2);
  double p = exp((cost1 - cost2) / a / T);
  return p;
}

// And here is the demonstration of how the algorithm is then used - note
// that you should call srand at the beginning of your code.
int main() {
  srand(time(NULL));

  std::vector<long> xs =
      { 0, 194, 908, 585, 666, 76, 633, 963, 789, 117, 409, 257, 229, 334, 837,
          382, 921, 54, 959, 532, 934, 720, 117, 519, 933, 408, 750, 465, 790,
          983, 605, 314, 272, 902, 340, 827, 915, 483, 466, 451, 698 };
  std::vector<long> ys = { 0, 956, 906, 148, 196, 59, 672, 801, 752, 620, 65,
      747, 377, 608, 374, 841, 910, 903, 743, 477, 794, 973, 555, 496, 152, 52,
      3, 174, 890, 861, 790, 430, 149, 674, 780, 507, 187, 931, 503, 435, 569 };
  Tour mytour(xs, ys);

  // Initialize parameters with default values;
  GenSimAnnealParams params = DEFAULT_SIM_ANNEALING_PARAMS;

  std::cout << "Enter max temps: ";
  std::cin >> params.max_temps;
  std::cout << "Enter iterations per temperature: ";
  std::cin >> params.iters_per_temp;
  std::cout << "Enter alpha: ";
  std::cin >> params.alpha;

  GenSimAnnealResults<Tour, unsigned long> results = genericSimulatedAnneal<
      Tour, unsigned long, acceptanceProbability>(mytour, params);

  std::cout << "Tour length: " << results.finalcost << std::endl;
  std::cout << "Function evaluations: " << results.function_evals << std::endl;
  std::cout << "Computed tour: ";
  for (size_t i = 0; i < results.best.visited.size(); ++i)
    std::cout << results.best.visited[i] << " ";
  std::cout << std::endl;

  return 0;
}
