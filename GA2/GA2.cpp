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


const int nofExperiments = 1;
const int nofRestarts = 1000;
const bool runMS = false;
const bool runGA = false;
const bool runILS = false;
const bool runPR = true;
const int perturbationSize = 1;

enum SearchType {
  MS, ILS, GA, PR
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
	ofstream GarrySolutionFile;
	GarrySolutionFile.open("GarySolution_" + to_string(time(0)) + ".txt");

	int best = numeric_limits<int>::max();
	for (int i = 0; i < nofRestarts; i++){
		auto solution = Chromosome(nodes.size());

		GarrySolutionFile << solution << endl;

		solution.swapNodesOpt();
//		cout << solution._score << ", ";

		GarrySolutionFile << solution << endl;

		if (solution._score < best) best = solution._score;
	}

	GarrySolutionFile.close();

//	cout << endl;
	ExperimentResult result;
	result.bestScore = best;
	return result;
}

Chromosome pathRelink(Chromosome a, Chromosome b, double truncate = 1) {

//	cout << a << endl;
//	cout << b << endl;

	//vector<Chromosome> steps;
	//steps.push_back(a);

	int distance = Chromosome::distance(a, b);

	cout << "distance: " << distance << endl;
//	cout << "steps: " << endl;

	Chromosome current = a;
	Chromosome bestCandidate = a;
	int currentBestScore = a._score;

	// we want to swap all mismatching bits
	for (int k = 0; k < distance; ++k) {

		// generate all swaps of one bit
		vector<Chromosome> candidates;
		candidates.reserve(a._solution.size());

		for (int i = 0; i < a._solution.size(); ++i) {
			//if (steps.back()._solution[i] != b._solution[i]) {
			if (current._solution[i] != b._solution[i]) {
				//Chromosome n(steps.back());
				Chromosome n(current);
				n.flipNodeAtIdx(i);
				candidates.push_back(n);
			}
		}
		
		//find the smallest score and its idx
		int minScore = numeric_limits<int>::max();
		int minId = 0;
		for (int j = 0; j < candidates.size(); ++j) {
			if (candidates[j]._score < minScore) {
				minId = j;
				minScore = candidates[j]._score;
			}
		}

		// evaluate all
	//	for (auto c : candidates) {
	//		cout << c << endl;
	//	}

//		cout << candidates[minId] << endl;
		//steps.push_back(candidates[minId]);
		current = candidates[minId];
		if (current._score < currentBestScore){ 
			currentBestScore = current._score;
			bestCandidate = current;
		}
	}

	return bestCandidate;
	/*
	int minScore = numeric_limits<int>::max();
	int minId = 0;
	for (int j = 0; j < steps.size(); ++j) {
		if (steps[j]._score < minScore) {
			minId = j;
			minScore = steps[j]._score;
		}
	}
	

//	cout << endl << "best: " << endl;
//	cout << steps[minId] << endl;

	return steps[minId];*/
}

ExperimentResult dynamicPathRelinking(int ESSize, int globalIter, int localIter, int dth = 20, double truncate = 1) {
	// - construct elite set (ES) with size b solutions
	vector<Chromosome> es;
	for (int i = 0; i < ESSize; ++i) {
		auto chrom = Chromosome::GRC();
		chrom.swapNodesOpt();
		es.push_back(chrom);
	}
	sort(es.begin(), es.end(), [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });

	// - sort ES by score

	for (int gi = 0; gi < globalIter; ++gi) { // or just set time limit
		for (int li = 0; li < localIter; ++li) {
			cout << "iteration: " << gi << " / " << li << endl;
//			cout << "gcrCalls: " << Chromosome::gcrCalls << endl;
			// - construct solution
			Chromosome x = Chromosome::GRC();
			// - local search
			x.swapNodesOpt();
			// - random select xj from ES
			int j = rand() % ESSize;

			// - get best solution from PR
			// - local search and save to y
			Chromosome y = Chromosome::PathRelink(x, es[j]);
			y.swapNodesOpt();

			if (y._score < es[0]._score) { // - if the solution is better than the best in ES replace it
				es[0] = y;
			}
			else if (y._score < es[ESSize - 1]._score) { // - elseif solution is better than the worst in ES and hamming distance is smaller than dth
				bool addToEs = true;
				vector<int> distances(ESSize);
				for (int i = 0; i < ESSize; ++i) {
					distances[i] = Chromosome::distance(y, es[i]);
					if (distances[i] <= dth) {
						addToEs = false;
						break;
					}
				}
				if (addToEs) {
					// find the closest solution in ES (by hamming distance) to y such that score y._score is better
					int swapId;
					int swapDistance = numeric_limits<int>::max();;
					for (int i = 0; i < ESSize; ++i) {
						if (y._score < es[i]._score) {
							if (distances[i] < swapDistance) {
								swapId = i;
								swapDistance = distances[i];
							}
						}
					}
					// replace them
					es[swapId] = y;
					// sort ES
					sort(es.begin(), es.end(), [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });
				}
			}
		}
		/*
		int map[ESSize][ESSize];
		for (int k = 0; k < ESSize; ++k) {
			for (int i = 0; i < ESSize; ++i) {
				map[k][i] = 0;
			}
		}

		bool newSol = true;
		while (newSol) {
			newSol = false;

			for (int i = 0; i < ESSize; ++i) {
				for (int j = 0; j < ESSize / 2; ++j) {

					// not need to generate those that were generated before
					if (map[i][j] == 1) {
						continue;
					}
					else {
						map[i][j] = 1;
					}

					cout << "combining " << i << " " << j << endl;

					Chromosome y = pathRelink(es[i], es[j], truncate);
					y.swapNodesOpt();

					if (y._score < es[0]._score) { // - if the solution is better than the best in ES replace it
						es[0] = y;
						newSol = true;
					}

					else if (y._score < es[ESSize - 1]._score) { // - elseif solution is better than the worst in ES and hamming distance is smaller than dth
						bool addToEs = true;
						vector<int> distances(ESSize);
						for (int i = 0; i < ESSize; ++i) {
							distances[i] = Chromosome::distance(y, es[i]);
							if (distances[i] <= dth) {
								addToEs = false;
								break;
							}
						}
						if (addToEs) {
							// find the closest solution in ES (by hamming distance) to y such that score y._score is better
							int swapId;
							int swapDistance = numeric_limits<int>::max();;
							for (int i = 0; i < ESSize; ++i) {
								if (y._score < es[i]._score) {
									if (distances[i] < swapDistance) {
										swapId = i;
										swapDistance = distances[i];
									}
								}
							}
							// replace them
							es[swapId] = y;
							// sort ES
							sort(es.begin(), es.end(), [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });

							newSol = true;
						}
					}
				}
			}
		}*/
	}

//
//	for (int k = 0; k < ESSize; ++k) {
//		cout << es[k] << endl;
//	}

	ExperimentResult result;
	result.bestScore = es[0]._score;
	cout << "valid: " << es[0] << endl;
	return result;
}


ExperimentResult gaSearch(vector<Node> nodes, int popSize) {
	//store optimal solution
	ofstream GASolutionFile;
	GASolutionFile.open("GASolution_" + to_string(time(0)) + ".txt");

	auto population = Chromosome::generateRandomPopulation(popSize, nodes.size());

	for (int i = 0; i < population.size(); i++){
		population[i].swapNodesOpt();
	}

	sort(population.begin(), population.end(), [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });

	GASolutionFile << population[0] << endl;

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
				GASolutionFile << population[0] << endl;
			}
		}

		//cout << "best: " << population[0]._score << " worst: " << population[population.size() - 1]._score << endl;
	} while (betterSolutionFound);

	GASolutionFile.close();
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

vector<ExperimentResult> runExperiments(vector<Node> nodes, int count, SearchType type, int p1 = -1, int p2 = -1, int p3 = -1, int p4 = -1) {

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
			case SearchType::PR:
				cout << "running PR " << i << ":" << p1 << ":" << p2 << ":" << p3 << ":" << p4 << endl;
				result = dynamicPathRelinking(p1, p2, p3, p4);
			break;
		}
		std::clock_t ms_end = std::clock();
		auto t2 = chrono::high_resolution_clock::now();

		result.cpuTime = 1000.0 * (ms_end - ms_start) / CLOCKS_PER_SEC;
		result.wallTime = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
		results.push_back(result);

		cout << "\tbest: " << result.bestScore << endl;
		cout << "\tcpu ime: " << result.cpuTime << endl;
	}

	ofstream output;
	int timestamp = std::time(0);
	switch (type) {
		case SearchType::MS:
			output.open(string("results/") + to_string(timestamp) + "_MS10000.csv");
		break;
		case SearchType::ILS:
			output.open(string("results/") + to_string(timestamp) + "_ILS" + to_string(p1) + ".csv");
		break;
		case SearchType::GA:
			output.open(string("results/") + to_string(timestamp) + "_GA" + to_string(p1) + ".csv");
		break;
		case SearchType::PR:
			output.open(string("results/") + to_string(timestamp) + "_PR_" + to_string(p1) + "_" + to_string(p2) + "_" + to_string(p3) + "_" + to_string(p4) + ".csv");
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


	int better = 0;
	for (int j = 0; j < 100; ++j) {
		auto a = Chromosome::GRC();
		a.swapNodesOpt();
		auto b = Chromosome::GRC();
		b.swapNodesOpt();

//		cout << "a: " << a << endl;
//		cout << "b: " << b << endl << endl;

//		auto grcSolution = Chromosome::GRC();
//		cout << grcSolution._score << " " << grcSolution.checkValidity() << endl;
		auto prChild = Chromosome::PathRelink(a, b);

//		cout << endl << prChild << endl;
//		cout << "valid: " << prChild.checkValidity() << endl;
//		cout << "score: " << prChild._score << endl;
		prChild.swapNodesOpt();
		cout << j << ": a=" << a._score << " b=" << b._score << " pr=" << prChild._score << endl;

		if (prChild._score < a._score && prChild._score < b._score) {
			better++;
		}
	}

	cout << "improvement in " << better << " out of 100 cases" << endl;


	return 0;


	auto optimalSolution = Chromosome(nodes.size(), Chromosome::OPTIMAL);
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

	if (runPR) runExperiments(nodes, nofExperiments, SearchType::PR, 30, 10, 20, 40);

	#ifdef _WIN64
		cin.ignore();
    #endif

	return 0;
}



//
//std::clock_t c_start = std::clock();
//std::clock_t c_end = std::clock();
//cout << "CPU time used: " << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n";


/*
int nofTests = 10000;
cout << "random:" << endl;
clock_t c_start = std::clock();
double count = 0;
for (int i = 0; i < nofTests; i++){
count += Chromosome(nodes.size(), Chromosome::RANDOM)._score;
}
std::clock_t end = std::clock();
cout << "CPU time used: " << 1000.0 * (end - c_start) / CLOCKS_PER_SEC << " ms, avg score: " << count / nofTests << endl;

cout << "greedy:" << endl;
c_start = std::clock();
count = 0;
for (int i = 0; i < nofTests; i++){
count += Chromosome(nodes.size(), Chromosome::GREEDY)._score;
}
end = std::clock();
cout << "CPU time used: " << 1000.0 * (end - c_start) / CLOCKS_PER_SEC << " ms, avg score: " << count / nofTests << endl;

cin.ignore();
return 0;
*/
