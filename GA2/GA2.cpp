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
#include <chrono>

using namespace std;

const int nofExperiments = 30;
const int nofRestarts = 1000;
const bool runMS = true;
const bool runGA = false;
const bool runILS = false;
const int perturbationSize = 1;

enum SearchType {
  MS, ILS, GA
};


struct ExperimentResult {
  int bestScore;
  clock_t cpuTime;
  chrono::milliseconds wallTime;
};

ExperimentResult summarizeExperiments(vector<ExperimentResult> exps){
	ExperimentResult result;

	result.cpuTime = 0;
	result.wallTime = chrono::milliseconds::zero();
	result.bestScore = 0;

	for (auto exp : exps){
		result.cpuTime += exp.cpuTime;
		result.bestScore += exp.bestScore;
		result.wallTime += exp.wallTime;
	}

	result.cpuTime /= exps.size();
	result.wallTime /= exps.size();
	result.bestScore /= exps.size();

	return result;
}

ostream& operator << (std::ostream &o, const ExperimentResult & rhs){
	o << rhs.bestScore << "," << rhs.cpuTime << "," << rhs.wallTime.count();
	return o;
}

ExperimentResult multiStart(vector<Node> nodes) {
	int best = numeric_limits<int>::max();
	for (int i = 0; i < nofRestarts; i++){
		auto solution = Chromosome(nodes.size());

		solution.swapNodesOpt();
//		cout << solution._score << ", ";

		if (solution._score < best) best = solution._score;
	}
//	cout << endl;
	ExperimentResult result;
	result.bestScore = best;
	return result;
}

ExperimentResult dynamicPathRelinking(int truncate, int ESSize, int globalIter, int localIter, int dth) {
	// - construct elite set (ES) with size b solutions
	// - sort ES by score

	for (int gi = 0; gi < globalIter; ++gi) { // or just set time limit
		for (int li = 0; li < localIter; ++li) {
			// - construct solution
			// - local search
			// - random select xj from ES
			// - get best solution from PR
			// - local search and save to y

			// - if the solution is better than the best in ES replace it
			// - elseif solution is better than the worst in ES and hamming distance is smaller than dth
				// find the closest solution in ES (by hamming distance) to y such that score y._score is better
				// replace them
				// sort ES
		}

		// recombine ES (not needed)
		// by PR till no new better solutions are found
	}

}

Chromosome pathRelink(Chromosome a, Chromosome b) {

}

ExperimentResult gaSearch(vector<Node> nodes, int popSize) {
	//store optimal solution
//	ofstream GASolutionFile;
//	GASolutionFile.open("GASolution.txt");

	auto population = Chromosome::generateRandomPopulation(popSize, nodes.size());

	for (int i = 0; i < population.size(); i++){
		population[i].swapNodesOpt();
	}

	sort(population.begin(), population.end(), [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });

	int worstPopulationScore;
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
				bestPopulationScore = population[0]._score;
			}
		}

		//cout << "best: " << population[0]._score << " worst: " << population[population.size() - 1]._score << endl;
	} while (betterSolutionFound);

//	GASolutionFile.close();
	ExperimentResult result;
	result.bestScore = bestPopulationScore;
	return result;
}

ExperimentResult iterativeLocalSearch(vector<Node> nodes, int perturbationSize) {
//	ofstream ILSsolution;
//	ILSsolution.open("ILSsolution.txt");
	int maxIterations = 50;
	double maxDeviation = 1.2;

	Chromosome candidate = Chromosome(nodes.size());
	candidate.swapNodesOpt();

	Chromosome bestSolution(candidate);

	do {
		candidate.mutate(perturbationSize);
		candidate.swapNodesOpt();

		//cout << "perturbation size: " << perturbationSize << " best score: " << bestSolution._score << " candidate score: " << candidate._score << " maxIter: " << maxIterations << " valid? " << candidate.checkValidity() << endl;

		if (candidate._score < bestSolution._score){
			bestSolution = candidate;
		}
		else if (candidate._score > bestSolution._score * maxDeviation){
			candidate = bestSolution;
		}
//		ILSsolution << candidate << endl;

		maxIterations--;
	} while (maxIterations > 0);

//	ILSsolution.close();
	ExperimentResult result;
	result.bestScore = bestSolution._score;
	return result;
}

vector<ExperimentResult> runExperiments(vector<Node> nodes, int count, SearchType type, int p1 = -1) {

	vector<ExperimentResult> results;

	for (int i = 0; i < count; ++i) {
		std::clock_t ms_start = std::clock();
		auto t1 = chrono::high_resolution_clock::now();
		ExperimentResult result;
		switch (type) {
			case SearchType::MS:
				cout << "running MS " << i << endl;
				result = multiStart(nodes);
			break;
			case SearchType::ILS:
				if (p1 == -1){
					throw runtime_error("ILS parameter cannot be negative!");
				}

				cout << "running ILS " << i << " and " << p1 << endl;
				result = iterativeLocalSearch(nodes, p1);
			break;
			case SearchType::GA:
				cout << "running GA " << i << " and " << p1 <<  endl;
				result = gaSearch(nodes, p1);
			break;
		}
		cout << "\tbest: " << result.bestScore << endl;
		std::clock_t ms_end = std::clock();
		auto t2 = chrono::high_resolution_clock::now();

		result.cpuTime = 1000.0 * (ms_end - ms_start) / CLOCKS_PER_SEC;
		result.wallTime = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
		results.push_back(result);
	}

	ofstream output;
	int timestamp = std::time(0);
	switch (type) {
		case SearchType::MS:
			output.open(string("results/") + to_string(timestamp) + "_MS1000.csv");
		break;
		case SearchType::ILS:
			output.open(string("results/") + to_string(timestamp) + "_ILS" + to_string(p1) + ".csv");
		break;
		case SearchType::GA:
			output.open(string("results/") + to_string(timestamp) + "_GA" + to_string(p1) + ".csv");
		break;
	}

	output << "score,cpu_time,wall_time" << endl;
	for (auto result : results) {
		output << result << endl;
	}
	output.close();
	return results;
}

int main(int argc, char* argv[])
{

	srand (time(NULL));

	vector<Node> nodes = DataReader::GetData("data.txt");
	Chromosome::_nodeList = nodes;

//	auto grcSolution = Chromosome::GRC();
//	cout << grcSolution._score << " " << grcSolution.checkValidity() << endl;




	auto optimalSolution = Chromosome(nodes.size(), true);
	cout << "Optimal soution of provided graph: " << optimalSolution._score << " isValid: " << optimalSolution.checkValidity() << endl;

	//store optimal solution
	ofstream optimalSolutionFile;
	optimalSolutionFile.open("optimalSolution.txt");
	optimalSolutionFile << optimalSolution;
	optimalSolutionFile.close();

	if (runMS) runExperiments(nodes, nofExperiments, SearchType::MS);
	if (runILS) {
		ofstream output;
		output.open(string("results/") + to_string(time(0)) + "_ILS_summary.csv");
		output << "perturbation,score,cpu_time,wall_time" << endl;
		for (int i = 2; i < 100; i++) {
			auto results = runExperiments(nodes, nofExperiments, SearchType::ILS, i);
			output << i << "," << summarizeExperiments(results) << endl;
		}
		output.close();
	}
	if (runGA) runExperiments(nodes, nofExperiments, SearchType::GA, 50);
	if (runGA) runExperiments(nodes, nofExperiments, SearchType::GA, 100);

	#ifdef _WIN64
		cin.ignore();
    #endif

	return 0;
}



//
//std::clock_t c_start = std::clock();
//std::clock_t c_end = std::clock();
//cout << "CPU time used: " << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n";