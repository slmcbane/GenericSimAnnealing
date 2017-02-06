/* 
 * Example use case for GenericSimAnneal, solving the classic NP-complete
 * problem of the traveling salesman who is trying to visit a set of cities
 * in a closed loop with the least possible distance covered. Please note
 * that the acceptance probability function and the cooling schedule here
 * have the standard forms - see 
 * <https://en.wikipedia.org/wiki/Simulated_annealing>; however, other forms
 * for these functions are easily implemented due to the generic nature of 
 * the implementation.
 *
 * Note the use of shared_ptr's to hold the vectors of coordinates and the
 * internal random number generator in the Tour class; this way the objects
 * are not deep copied at every iteration as they would be otherwise.
 *
 * This code is released as an example under no copyright whatsoever.
 */

#include "../generic_simulated_anneal.hpp"
#include <vector>
#include <cassert>
#include <cmath>
#include <chrono>
#include <memory>
#include <random>

using std::chrono::system_clock;

// A class representing a particular tour.
// xcoords and ycoords hold the x- and y- coordinates of all the cities.
// visited is a vector of the order in which cities are visited (beginning
// and ending with city # 0).
class Tour
{
private:
	std::shared_ptr<std::mt19937> generator;
public:
	size_t numcities;
	std::shared_ptr<const std::vector<long>> xcoords, ycoords;
	std::vector<unsigned> visited;

	Tour(const std::vector<long>& xs, const std::vector<long>& ys);
	unsigned long cost() const;
	void perturb();
};

Tour::Tour(const std::vector<long>& xs, const std::vector<long>& ys)
{
	assert(xs.size() == ys.size());
	numcities = xs.size();
	xcoords = std::make_shared<const std::vector<long>>(xs);
	ycoords = std::make_shared<const std::vector<long>>(ys);
	visited = std::vector<unsigned>(numcities + 1);
	for (unsigned i = 0; i < numcities; ++i)
		visited[i] = i;
	visited[numcities] = 0;
	generator = std::make_shared<std::mt19937>(
		system_clock::to_time_t(system_clock::now()));
}

// Computes the cost of the tour as given by the floor of the
// Pythagorean distance.
unsigned long Tour::cost() const
{
	unsigned long dist = 0;
	const std::vector<long>& x = *xcoords;
	const std::vector<long>& y = *ycoords;
	for (size_t i = 0; i < numcities; ++i) {
		long dx = x[visited[i + 1]] - x[visited[i]];
		long dy = y[visited[i + 1]] - y[visited[i]];
		double d = sqrt(double(dx * dx + dy * dy));
		dist += (unsigned long)(floor(d));
	}
	return dist;
}

// Swaps two random cities on the tour; this is a good enough perturbation
// function for an example case.
void Tour::perturb()
{
	auto& gen = *generator;
	unsigned city1 = 1 + unsigned(floor((double(gen()) * (numcities - 1))
										/ gen.max()));
	unsigned city2 = city1;
	while (city2 == city1)
		city2 = 1 + unsigned(floor((double(gen()) * (numcities - 1))
								   / gen.max()));

	unsigned temp = visited[city1];
	visited[city1] = visited[city2];
	visited[city2] = temp;
}

// And here is the demonstration of how the algorithm is then used - remember
// to call srand at the beginning of your code; even if your perturbation 
// function does something else the annealing routine uses rand() to decide
// whether or not to accept an increase in cost.
int main()
{
	std::vector<long> xs =
	{ 0, 194, 908, 585, 666, 76, 633, 963, 789, 117, 409, 257, 229, 334, 837,
		382, 921, 54, 959, 532, 934, 720, 117, 519, 933, 408, 750, 465, 790,
		983, 605, 314, 272, 902, 340, 827, 915, 483, 466, 451, 698 };
	std::vector<long> ys = { 0, 956, 906, 148, 196, 59, 672, 801, 752, 620, 65,
		747, 377, 608, 374, 841, 910, 903, 743, 477, 794, 973, 555, 496, 152, 52,
		3, 174, 890, 861, 790, 430, 149, 674, 780, 507, 187, 931, 503, 435, 569 };
	Tour mytour(xs, ys);

	GenSimAnnealParams params;
	double alpha;

	std::cout << "Enter max temps: ";
	std::cin >> params.max_temps;
	std::cout << "Enter iterations per temperature: ";
	std::cin >> params.iters_per_temp;
	std::cout << "Enter alpha: ";
	std::cin >> alpha;
	params.verbose = false;
	params.cost_reduction_tol = 0.0;

	// This lambda is the cooling schedule; for this example, temperature is
	// reduced by a fixed multiplier alpha at each iteration.
	auto s = [alpha] (unsigned iter) { return 1.0 * pow(alpha, iter); };

	// This lambda is our acceptance probability function. Note that this is the 
	// same as the 'canonical' probability function often found in simulated
	// annealing implementations except for the addition of a scaling factor for
	// temperature (since T is always 1.0 or less).
	auto a = [] (unsigned long c1, unsigned long c2, double T) {
		return exp((double(c1) - double(c2)) / T / 600.0);
	};

	std::mt19937 generator(system_clock::to_time_t(system_clock::now()));

	// This is how the routine can be called if a lambda expression is used for
	// the acceptance function; it's very similar for the case of a functor or a 
	// normal function declaration.
	Tour result =
		genericSimulatedAnneal<Tour, unsigned long, decltype(a), decltype(s), std::mt19937>
		(mytour, params, a, s, generator);

	std::cout << "Tour length: " << result.cost() << std::endl;
	std::cout << "Computed tour: ";
	for (auto i : result.visited)
		std::cout << i << " ";
	std::cout << std::endl;

	return 0;
}
