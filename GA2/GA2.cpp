// GA2.cpp : Defines the entry point for the console application.
//

//#pragma warning(disable:4996)
//#define _SCL_SECURE_NO_WARNINGS


#include <iostream>
#include <algorithm>
#include <vector>

#include "DataReader.h"
#include "Chromosome.h"


using namespace std;

const int nofRestarts = 1000;
const bool runMultistart = false;
const bool runIteratedLocalSearch = false;

int main(int argc, char* argv[])
{

	vector<Node> nodes = DataReader::GetData("data.txt");
	Chromosome::_nodeList = nodes;

	auto optimalSolution = Chromosome(nodes.size(), true);
	cout << "Optimal soution of provided graph: " << optimalSolution._score << " isValid: " << optimalSolution.checkValidity() << endl;

	//multistart random initial solutions + local search
	if (runMultistart) {
		int best = numeric_limits<int>::max();
		for (int i = 0; i < nofRestarts; i++){
			auto solution = Chromosome(nodes.size());

			int oldScore = solution._score;
			int improvements = solution.swapNodesOpt();
			cout << oldScore << " " << improvements << " : " << solution._score << endl;

			if (solution._score < best) best = solution._score;
		}
		cout << "Multistart best solution score: " << best << endl;
	}

	//ILS
	//TODO perturbation sizes
	if (runIteratedLocalSearch) {
		for (int i = 0; i < 1; i++){
			Chromosome candidate = Chromosome(nodes.size());
			candidate.swapNodesOpt();

			cout << "starting score: " << candidate._score << endl;
			for (int attempts = 0; attempts < 10; attempts++){
				candidate.mutate(0.1);
				cout << "after mutation: " << candidate._score << endl;
				candidate.swapNodesOpt();
				cout << "result score: " << candidate._score << " valid? " << candidate.checkValidity() << endl;
			}

		}
	}

	

	if (true) {
		auto population = Chromosome::generateRandomPopulation(50, nodes.size());

		sort(population.begin(), population.end(), [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });

		int worstPopulationScore = population[population.size() - 1]._score;

		bool betterSolutionFound;
		do {
			betterSolutionFound = false;

			// tornament selection of parrents
			auto parentsA = Chromosome::GATournamentSelection(population, 2);
			auto parentsB = Chromosome::GATournamentSelection(population, 2);

			auto children = Chromosome::GAGenerateChildren(parentsA, parentsB);

			int bestChildScore = numeric_limits<int>::max();
			for (int i = 0; i < children.size(); i++){
				children[i] = children[i].swapNodesOpt();
				bestChildScore = min(bestChildScore, children[i]._score);
			}

			if (bestChildScore > worstPopulationScore){
				betterSolutionFound = true;

				//sort children
				sort(children.begin(), children.end(), [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });

				//add children to the parent population, both population and children are sorted so we can use merge with O(n+m) complexity
				vector<Chromosome> combinedPop;
				merge(population.begin(), population.end(), children.begin(), children.end(), combinedPop, [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });

				//get the best population.size() chromosomes
				population = vector<Chromosome>(combinedPop.begin(), combinedPop.begin() + population.size());

				//eletist selection
				


				// add children to population
				// sort population
				// truncate population
			}
			

		} while (betterSolutionFound);

	}


	cout << "done" << endl;
	#ifdef _WIN64
		cin.ignore();
    #endif

	return 0;
}

