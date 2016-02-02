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

#include <thread>
#include <future>

using namespace std;

int MAX_TIME = -1;
int MAX_ITER = -1;

const int ONE_THREAD = false;

const int nofExperiments = 30;
const bool runMS = false;
const bool runGA = false;
const bool runILS = false;
const bool runPR = true;

enum SearchType {
  MS, ILS, GA, PR
};
static const char * EnumStrings[] = { "MS", "ILS", "GA", "PR" };
static const char * GenTypeStrings[] = { "OPTIMAL", "RANDOM", "GREEDY" };


struct ExperimentResult {
  int bestScore;
  clock_t cpuTime;
  chrono::milliseconds wallTime;
  int iterations;
};

ExperimentResult summarizeExperiments(vector<ExperimentResult> exps){
	ExperimentResult result;

	result.cpuTime = 0;
	result.wallTime = chrono::milliseconds::zero();
	result.bestScore = 0;
	result.iterations = 0;

	for (auto exp : exps){
		result.cpuTime += exp.cpuTime;
		result.bestScore += exp.bestScore;
		result.wallTime += exp.wallTime;
		result.iterations += exp.iterations;
	}

	result.cpuTime /= exps.size();
	result.wallTime /= exps.size();
	result.bestScore /= exps.size();
	result.iterations /= exps.size();

	return result;
}

ostream& operator << (std::ostream &o, const ExperimentResult & rhs){
	o << rhs.bestScore << ", " << rhs.cpuTime << ", " << rhs.wallTime.count() << ", " << rhs.iterations;
	return o;
}

ExperimentResult multiStart_t(int size, Chromosome::GenerationType genType, int iterations, time_t duration){
	time_t endTime = time(0) + duration;
	if (duration < 0){
		endTime = numeric_limits<long>::max();
	}

	int best = numeric_limits<int>::max();
	int i;
	for (i = 0; i < iterations && time(0) < endTime; i++){
		auto solution = Chromosome(size, genType);
		solution.swapNodesOpt();

		if (solution._score < best) best = solution._score;
	}

	ExperimentResult result;
	result.bestScore = best;
	result.iterations = i;
	return result;
}

ExperimentResult multiStart(int size, Chromosome::GenerationType genType, int iterations = -1, time_t maxTime = -1) {
	if (iterations < 0) iterations = numeric_limits<int>::max();
	time_t endTime = time(0) + maxTime;
	if (maxTime < 0){
		endTime = numeric_limits<long>::max();
	}

	int best = numeric_limits<int>::max();
	int i;
	for (i = 0; i < iterations && time(0) < endTime; i++){
		auto solution = Chromosome(size, genType);
		solution.swapNodesOpt();

		if (solution._score < best) best = solution._score;
	}

	ExperimentResult result;
	result.bestScore = best;
	result.iterations = i;
	return result;
}

Chromosome pathRelink(Chromosome a, Chromosome b, double truncate = 1) {

//	cout << a << endl;
//	cout << b << endl;

	//vector<Chromosome> steps;
	//steps.push_back(a);

	int distance = Chromosome::distance(a, b);

//	cout << "distance: " << distance << endl;
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

ExperimentResult dynamicPathRelinking(int ESSize, int globalIter, int localIter, int dth, Chromosome::GenerationType genType, double truncate = 1, time_t maxTime = MAX_TIME) {
	// - construct elite set (ES) with size b solutions
	vector<Chromosome> es;
	for (int i = 0; i < ESSize; ++i) {
		Chromosome chrom(Chromosome::_nodeList.size(), genType);
		chrom.swapNodesOpt();
		es.push_back(chrom);
	}
	sort(es.begin(), es.end(), [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });

	// - sort ES by score

	time_t endTime = time(0) + maxTime;
	if (maxTime < 0){
		endTime = numeric_limits<long>::max();
	}

	int gi;
	for (gi = 0; gi < globalIter && time(0) < endTime; ++gi) { // or just set time limit
//		cout << "iteration: " << gi << endl;
		for (int li = 0; li < localIter; ++li) {
//			cout << "iteration: " << gi << " / " << li << endl;
//			cout << "gcrCalls: " << Chromosome::gcrCalls << endl;
			// - construct solution
			Chromosome x(Chromosome::_nodeList.size(), genType);

			// - local search
			x.swapNodesOpt();
			// - random select xj from ES
			int j = rand() % ESSize;

			// - get best solution from PR
			// - local search and save to y
			Chromosome y = Chromosome::PathRelink(x, es[j]);
//			Chromosome y = Chromosome::GACrossOver(x, es[j]);
			y.swapNodesOpt();

//			if (y._score > x._score) { //uncomment for GACrossover
//				y = x;
//			}

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
//					cout << "better ES member recombined" << endl;
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

		int map[30][30];
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

//					cout << "combining " << i << " " << j << endl;

					Chromosome y = Chromosome::PathRelink(es[i], es[j]);
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
		}
	}

//
//	for (int k = 0; k < ESSize; ++k) {
//		cout << es[k] << endl;
//	}

	ExperimentResult result;
	result.bestScore = es[0]._score;
	result.iterations = gi;
//	cout << "valid: " << es[0] << endl;
	return result;
}


ExperimentResult gaSearch(int size, int popSize, Chromosome::GenerationType genType, time_t maxTime = -1) {
	time_t endTime = time(0) + maxTime;
	if (maxTime < 0){
		endTime = numeric_limits<long>::max();
	}

	//store optimal solution
	//ofstream GASolutionFile;
	//GASolutionFile.open("GASolution_" + to_string(time(0)) + ".txt");

	auto population = Chromosome::generateRandomPopulation(popSize, size, genType);

	for (int i = 0; i < population.size(); i++){
		population[i].swapNodesOpt();
	}

	sort(population.begin(), population.end(), [](const Chromosome & a, const Chromosome & b) {return a._score < b._score; });

	//GASolutionFile << population[0] << endl;

	int worstPopulationScore;
	int bestPopulationScore = population[0]._score;

	int generations = 0;
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
				//GASolutionFile << population[0] << endl;
			}
		}

		generations++;

		//cout << "best: " << population[0]._score << " worst: " << population[population.size() - 1]._score << endl;
	} while (betterSolutionFound && time(0) < endTime);

	//GASolutionFile.close();
	ExperimentResult result;
	result.bestScore = bestPopulationScore;
	result.iterations = generations;
	return result;
}

//maxtime in seconds
ExperimentResult iterativeLocalSearch(int size, int perturbationSize, Chromosome::GenerationType genType, int maxIterations = 50, time_t maxTime = -1) {
//	ofstream ILSsolution;
//	ILSsolution.open("ILSsolution.txt");
	double maxDeviation = 1.2;

	Chromosome candidate = Chromosome(size, genType);
	candidate.swapNodesOpt();

	Chromosome bestSolution(candidate);

	time_t endTime = time(0) + maxTime;
	if (maxTime < 0){
		endTime = numeric_limits<long>::max();
	}

	if (maxIterations < 0){
		maxIterations = numeric_limits<int>::max();
	}

	int i = 0;
	while (maxIterations > i && time(0) < endTime) {
		candidate.mutate(perturbationSize);
		candidate.swapNodesOpt();

		//if the score is very bad -> use the old optimal solution to perturbate from
		if (candidate._score < bestSolution._score){
			bestSolution = candidate;
		}
		else if (candidate._score > bestSolution._score * maxDeviation){
			candidate = bestSolution;
		}
//		ILSsolution << candidate << endl;

		i++;
	}

//	ILSsolution.close();
	ExperimentResult result;
	result.bestScore = bestSolution._score;
	result.iterations = i;
	return result;
}

ExperimentResult runExperiment(int size, SearchType type, Chromosome::GenerationType genType, int p1, int p2, int p3, int p4){
	std::clock_t ms_start = std::clock();
	auto t1 = chrono::high_resolution_clock::now();
	ExperimentResult result;
	switch (type) {
	case SearchType::MS:
		//cout << "running MS " << i << endl;
		result = multiStart(size, genType, MAX_ITER, MAX_TIME);
		break;
	case SearchType::ILS:
		if (p1 == -1){
			throw runtime_error("ILS parameter cannot be negative!");
		}

		//cout << "running ILS " << i << " and " << p1 << endl;
		result = iterativeLocalSearch(size, p1, genType, MAX_ITER, MAX_TIME);
		break;
	case SearchType::GA:
		//cout << "running GA " << i << " and " << p1 << endl;
		result = gaSearch(size, p1, genType, MAX_TIME);
		break;
	case SearchType::PR:
		//cout << "running PR " << i << ":" << p1 << ":" << p2 << ":" << p3 << ":" << p4 << endl;
		result = dynamicPathRelinking(p1, p2, p3, p4, genType);
		break;
	}
	std::clock_t ms_end = std::clock();
	auto t2 = chrono::high_resolution_clock::now();

	result.cpuTime = 1000.0 * (ms_end - ms_start) / CLOCKS_PER_SEC;
	result.wallTime = chrono::duration_cast<chrono::milliseconds>(t2 - t1);

	return result;
}

vector<ExperimentResult> runExperiments(vector<Node> nodes, int count, SearchType type, Chromosome::GenerationType genType, int p1 = -1, int p2 = -1, int p3 = -1, int p4 = -1) {

	cout << "experiment parameters: " << p1 << ":" << p2 << ":" << p3 << ":" << p4 << endl;

	vector<ExperimentResult> results;
	int nodesSize = nodes.size();

	int nofThreads;
	if (ONE_THREAD) {
		nofThreads = 1;
	}
	else {
		nofThreads = max((unsigned int)1, std::thread::hardware_concurrency() - 1);
	}
	vector<future<ExperimentResult>> futures(nofThreads);
	
	int experimentsRun = 0;
	while (experimentsRun < count){
		cout << "running experiments " << experimentsRun << "..." << (experimentsRun + nofThreads) << " type: " << EnumStrings[type] << endl;

		for (int i = 0; i < nofThreads; i++){
			futures[i] = async(runExperiment, nodesSize, type, genType, p1, p2, p3, p4);
		}

		for (int i = 0; i < nofThreads; i++){
			auto result = futures[i].get();

			results.push_back(result);
		}

		experimentsRun += nofThreads;
	}

	ofstream output;
	int timestamp = std::time(0);
	string path = "results/" + to_string(timestamp) + "_" + EnumStrings[type] + "_" + GenTypeStrings[genType] + "_" + to_string(p1) + "_" + to_string(p2) + "_" + to_string(p3) + "_" + to_string(p4) + "_t" + to_string(MAX_TIME) + "_i" + to_string(MAX_ITER) + ".csv";
	output.open(path);

	cout << "creating result file: " << path << endl;

	output << "score,cpu_time,wall_time,iterations" << endl;
	for (auto result : results) {
		output << result << endl;
	}
	output.close();


	return results;
}

int main(int argc, char* argv[])
{
	Chromosome::GenerationType genType = Chromosome::GenerationType::GREEDY;
	if (argc > 1){
		if (string(argv[1]) == "greedy")
			genType = Chromosome::GenerationType::GREEDY;
		else {
			genType = Chromosome::GenerationType::RANDOM;
		}
		if (argc > 2)
			MAX_TIME = atoi(argv[2]);

	}

	srand (time(NULL));
	vector<Node> nodes = DataReader::GetData("data.txt");
	Chromosome::_nodeList = nodes;

	auto optimalSolution = Chromosome(nodes.size(), Chromosome::OPTIMAL);
	cout << "Optimal soution of provided graph: " << optimalSolution._score << " isValid: " << optimalSolution.checkValidity() << endl;

	//store optimal solution
	ofstream optimalSolutionFile;
	optimalSolutionFile.open("optimalSolution.txt");
	optimalSolutionFile << optimalSolution;
	optimalSolutionFile.close();

	vector<Chromosome::GenerationType> genTypes = {Chromosome::GenerationType::GREEDY, Chromosome::GenerationType::RANDOM};
	//for (auto genType : genTypes){
		cout << "using genType: " << GenTypeStrings[genType] << endl;

		if (runMS) runExperiments(nodes, nofExperiments, SearchType::MS, genType);

		if (runGA) runExperiments(nodes, nofExperiments, SearchType::GA, genType, 20);
		if (runGA) runExperiments(nodes, nofExperiments, SearchType::GA, genType, 50);
		if (runGA) runExperiments(nodes, nofExperiments, SearchType::GA, genType, 100);

		if (runILS) {
			/*ofstream output;
			output.open(string("results/") + to_string(time(0)) + "_ILS_summary.csv");
			output << "perturbation,score,cpu_time,wall_time" << endl;
			int i;
			for (i = 2; i <= 100; i++) {
				auto results = runExperiments(nodes, nofExperiments, SearchType::ILS, genType, i);
				output << i << "," << summarizeExperiments(results) << endl;
			}

			output.close();

			runExperiments(nodes, nofExperiments, SearchType::ILS, genType, 25);
			runExperiments(nodes, nofExperiments, SearchType::ILS, genType, 50);
			runExperiments(nodes, nofExperiments, SearchType::ILS, genType, 100);*/
			runExperiments(nodes, nofExperiments, SearchType::ILS, genType, 30);
		}

//		if (runPR) runExperiments(nodes, nofExperiments, SearchType::PR, genType, 3, 5, 200, 1);
//		if (runPR) runExperiments(nodes, nofExperiments, SearchType::PR, genType, 5, 5, 200, 1);
//		if (runPR) runExperiments(nodes, nofExperiments, SearchType::PR, genType, 10, 5, 200, 1);
//		if (runPR) runExperiments(nodes, nofExperiments, SearchType::PR, genType, 20, 5, 200, 1);
//		if (runPR) runExperiments(nodes, nofExperiments, SearchType::PR, genType, 30, 5, 200, 10);
//		if (runPR) runExperiments(nodes, nofExperiments, SearchType::PR, genType, 50, 5, 200, 1); // todo run solo
//		if (runPR) runExperiments(nodes, nofExperiments, SearchType::PR, genType, 100, 5, 200, 1); // todo run solo

//		if (runPR) runExperiments(nodes, nofExperiments, SearchType::PR, Chromosome::GenerationType::GREEDY, 30, 5, 200, 10);

		if (runPR) runExperiments(nodes, nofExperiments, SearchType::PR, Chromosome::GenerationType::GREEDY, 50, 5, 200, 1);
		if (runPR) runExperiments(nodes, nofExperiments, SearchType::PR, Chromosome::GenerationType::GREEDY, 100, 5, 200, 1);
		if (runPR) runExperiments(nodes, nofExperiments, SearchType::PR, Chromosome::GenerationType::GREEDY, 200, 5, 200, 1);

	//}

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



/*
ofstream randomTimes;
randomTimes.open("results/" + to_string(time(0)) + "randomTimes.csv");
int tests = 10000;

clock_t ms_start = clock();
for (int i = 0; i < tests; i++){
Chromosome test = Chromosome(nodes.size(), Chromosome::GenerationType::RANDOM);
test.swapNodesOpt();
randomTimes << test._score << endl;
}
randomTimes << endl << (double)(clock() - ms_start) / CLOCKS_PER_SEC << endl;

randomTimes.close();

ofstream greedyTimes;
greedyTimes.open("results/" + to_string(time(0)) + "greedyTimes.csv");

ms_start = clock();
for (int i = 0; i < tests; i++){
Chromosome test = Chromosome(nodes.size(), Chromosome::GenerationType::GREEDY);
test.swapNodesOpt();
greedyTimes << test._score << endl;
}

greedyTimes << endl << (double)(clock() - ms_start) / CLOCKS_PER_SEC << endl;

greedyTimes.close();

return 1;
*/