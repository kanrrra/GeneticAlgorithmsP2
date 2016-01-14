// GA2.cpp : Defines the entry point for the console application.
//

//#pragma warning(disable:4996)
//#define _SCL_SECURE_NO_WARNINGS



#include "DataReader.h"
#include "Chromosome.h"

#include <iostream>
#include <algorithm>
#include <vector>

#include <ctime>

using namespace std;

const int nofExperiments = 30;
const int nofRestarts = 1000;
const bool runMultistart = true;
const bool runGA = false;
const bool runIteratedLocalSearch = false;
const int perturbationSize = 1;

int main(int argc, char* argv[])
{

	vector<Node> nodes = DataReader::GetData("data.txt");
	Chromosome::_nodeList = nodes;

	auto optimalSolution = Chromosome(nodes.size(), true);
	cout << "Optimal soution of provided graph: " << optimalSolution._score << " isValid: " << optimalSolution.checkValidity() << endl;
	
	//store optimal solution
	ofstream optimalSolutionFile;
	optimalSolutionFile.open("optimalSolution.txt");
	optimalSolutionFile << optimalSolution;
	optimalSolutionFile.close();


	//multistart random initial solutions + local search
	if (runMultistart) {
		vector<clock_t> msTimes;
		for (int experimentNo = 0; experimentNo < nofExperiments; experimentNo++){
			std::clock_t ms_start = std::clock();

			int best = numeric_limits<int>::max();
			for (int i = 0; i < nofRestarts; i++){
				auto solution = Chromosome(nodes.size());

				int oldScore = solution._score;
				int improvements = solution.swapNodesOpt();

				if (solution._score < best) best = solution._score;
			}
			std::clock_t ms_end = std::clock();
			msTimes.push_back(1000.0 * (ms_end - ms_start) / CLOCKS_PER_SEC);

			cout << "Multistart best solution score for experiment " << experimentNo << ": " << best << endl;
			cout << "CPU time used: " << 1000.0 * (ms_end - ms_start) / CLOCKS_PER_SEC << " ms\n";
		}
	}


	//ILS
	//TODO perturbation sizes
	if (runIteratedLocalSearch) {
		ofstream ILSsolution;
		ILSsolution.open("ILSsolution.txt");
		for (int i = 0; i < 10; i++){
			Chromosome candidate = Chromosome(nodes.size());
			candidate.swapNodesOpt();

			int oldScore = candidate._score;
			cout << "starting score: " << candidate._score << endl;
			do {
				oldScore = candidate._score;
				candidate.mutate(i);
				cout << "after mutation: " << candidate._score << endl;
				candidate.swapNodesOpt();
				cout << "result score: " << candidate._score << " valid? " << candidate.checkValidity() << endl;

				ILSsolution << candidate << endl;

			} while (candidate._score < oldScore);

		}
		ILSsolution.close();

	}

	if (runGA) {
		//store optimal solution
		ofstream GASolutionFile;
		GASolutionFile.open("GASolution.txt");

		auto population = Chromosome::generateRandomPopulation(50, nodes.size());

		for (int i = 0; i < population.size(); i++){
			population[i].swapNodesOpt();
		}

		sort(population.begin(), population.end(), [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });

		int worstPopulationScore = population[population.size() - 1]._score;
		int bestPopulationScore = population[0]._score;

		bool betterSolutionFound;
		do {
			
			betterSolutionFound = false;
			worstPopulationScore = population[population.size() - 1]._score;

			// tornament selection of parrents
			auto parentsA = Chromosome::GATournamentSelection(population, 2);
			auto parentsB = Chromosome::GATournamentSelection(population, 2);

			auto children = Chromosome::GAGenerateChildren(parentsA, parentsB);

			int bestChildScore = numeric_limits<int>::max();
			for (int i = 0; i < children.size(); i++){
				children[i].swapNodesOpt();
				bestChildScore = min(bestChildScore, children[i]._score);
			}

			if (bestChildScore < worstPopulationScore){
				betterSolutionFound = true;

				//sort children
				sort(children.begin(), children.end(), [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });

				//add children to the parent population, both population and children are sorted so we can use merge with O(n+m) complexity
				vector<Chromosome> combinedPop;
				std::merge(population.begin(), population.end(), children.begin(), children.end(), back_inserter(combinedPop), [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });

				//get the best population.size() chromosomes
				population = vector<Chromosome>(combinedPop.begin(), combinedPop.begin() + population.size());

				if (population[0]._score < bestPopulationScore){
					GASolutionFile << population[0] << endl;
					bestPopulationScore = population[0]._score;
				}
			}
			
			cout << "best: " << population[0]._score << " worst: " << population[population.size() - 1]._score << endl;
		} while (betterSolutionFound);

		GASolutionFile.close();
	}


	cout << "done" << endl;
	#ifdef _WIN64
		cin.ignore();
    #endif

	return 0;
}
//
//std::clock_t c_start = std::clock();
//std::clock_t c_end = std::clock();
//cout << "CPU time used: " << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n";