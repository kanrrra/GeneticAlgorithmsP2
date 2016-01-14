//#include "stdafx.h"
#include "Chromosome.h"
#include <stdexcept>
#include <numeric>

vector<char> Chromosome::defaultDistribution;
vector<Node> Chromosome::_nodeList;


vector<char> Chromosome::generateOptimalSolutionForAssignmentData(int size){
	vector<char> sol(size, 0);

	for (Node n : _nodeList){
		if (n._y > 0.48){
			sol[n._id] = 1;
		} else if (n._y > 0.39 && n._x > 0.7){
			sol[n._id] = 1;
		}
		if (n._x < 0.12 && n._y > 0.75 && n._y < 0.92){
			sol[n._id] = 0;
		}
		if (n._x > 0.6 && n._x < 0.7 && n._y > 0.5 && n._y < 0.6){
			sol[n._id] = 0;
		}
	}

	return sol;
}

vector<char> Chromosome::generateRandomSolution(int size){
	if (defaultDistribution.size() != size){
		vector<char> defaultDistribution0(size / 2, 0);
		vector<char> defaultDistribution1(size / 2, 1);

		defaultDistribution.reserve(size);

		defaultDistribution.insert(defaultDistribution.end(), defaultDistribution0.begin(), defaultDistribution0.end());
		defaultDistribution.insert(defaultDistribution.end(), defaultDistribution1.begin(), defaultDistribution1.end());
	}

	vector<char> solution = defaultDistribution;

	random_shuffle(solution.begin(), solution.end());

	return solution;
}

bool Chromosome::checkValidity(){
	int count = 0;
	for (int sol : _solution){
		if (sol == 1) count++;
		else if (sol != 0) throw domain_error("Invalid value exception");
	}

	return (count * 2) == _solution.size();
}

vector<Chromosome> Chromosome::generateRandomPopulation(int populationSize, int solutionSize){
	vector<Chromosome> population;
	population.reserve(populationSize);

	for (int i = 0; i < populationSize; i++){
		population.push_back(Chromosome(solutionSize));
	}

	return population;
}

int Chromosome::mutate(unsigned int perturbationSize){
	//static geometric_distribution<int> dist(p); // todo: try different distributions and parameters -> report in the paper
	static mt19937 generator;

	uniform_int_distribution<int> bitPosDist(0, _solution.size() - 1);

	//unsigned int nofbitflips = dist(generator);
	//cout << "bitflips: " << nofbitflips << endl;

	//nofbitflips = min((unsigned int)_solution.size(), nofbitflips);

	unsigned int nofbitflips = perturbationSize;

	vector<int> indices(nofbitflips*2);
	unsigned int counters[2] = { 0, 0 };

	for (unsigned int i = 0; i < nofbitflips*2; i++){
		while (true){
			int pos = bitPosDist(generator);

			int group = _solution[pos];
			if (counters[group] >= nofbitflips) continue;

			if (find(indices.begin(), indices.end(), pos) == indices.end()) {
				indices[i] = pos;
				counters[group]++;
				break;
			}
		}

		_solution[indices[i]] ^= 1;
	}

	int oldScore = _score;

	_score = calcScore();

	return oldScore - _score;
}

int Chromosome::swapNodesOpt(){
	int improvementsCount = 0;
	bool improvementFound;
	do {
		improvementFound = false;

		for (int i = 0; i < _solution.size() - 1; i++){

			int localScoreI = 0;
			//for (auto nb : _nodeList[i]._links){
			//	localScoreI += (_solution[i] != _solution[nb]);
			//}
			//cout << (int)_scoreContribution[i] << endl;


			localScoreI = _scoreContribution[i];
			localScoreI = 2 * localScoreI - _nodeList[i]._links.size();


			for (int j = i + 1; j < _solution.size(); j++){
				if (_solution[i] != _solution[j]){
					
					int localScoreJ = 0;
					//for (auto nb : _nodeList[j]._links){
					//	localScoreJ += (_solution[j] != _solution[nb]);
					//}

					localScoreJ = _scoreContribution[j];
					localScoreJ = 2 * localScoreJ - _nodeList[j]._links.size();

					int improvement = localScoreI + localScoreJ;

					//no improvement
					if (improvement <= 0){
						continue;
					}

					//check if the swapped nodes link to each other
					if (find(_nodeList[i]._links.begin(), _nodeList[i]._links.end(), j) != _nodeList[i]._links.end()){
						improvement -= 2;
					}

					//no improvement
					if (improvement <= 0){
						continue;
					}

					_score -= improvement;

					//swap nodes
					//_solution[i] ^= 1;
					//_solution[j] ^= 1;

					flipNodeAtIdx(i);
					flipNodeAtIdx(j);

					improvementsCount++;
					improvementFound = true;
					break;
				}
			}
		}


	} while (improvementFound);

	int testScore = calcScore();
	if (testScore != _score){
		throw runtime_error("bad score exception");
	}

	return improvementsCount;
}

void Chromosome::flipNodeAtIdx(int idx){
	char myColor = _solution[idx];
	for (int n : _nodeList[idx]._links){
		if (_solution[n] == myColor){
			_scoreContribution[n]++;
		}
		else {
			_scoreContribution[n]--;
		}
	}

	_scoreContribution[idx] = _nodeList[idx]._links.size() - _scoreContribution[idx];
	_solution[idx] ^= 1;
}

int Chromosome::calcScore(){
	_scoreContribution.resize(_nodeList.size());

	for (int i = 0; i < _nodeList.size(); i++){
		Node n = _nodeList[i];

		char myColor = _solution[n._id];
		int localScore = 0;
		for (int meighbourId : n._links){
			localScore += (myColor != _solution[meighbourId]);
		}
		_scoreContribution[i] = localScore;
	}

	int score = accumulate(_scoreContribution.begin(), _scoreContribution.end(), 0);

	if (score % 2) throw runtime_error("Score cannot be odd before division");

	score /= 2;

	return score;
}

Chromosome::Chromosome(int size, bool optimal)
{
	if (optimal) _solution = generateOptimalSolutionForAssignmentData(size);
 	else _solution = generateRandomSolution(size);
	
	_score = calcScore();
}

Chromosome::Chromosome(vector<char> & solution)
{
	_solution = solution;
	_score = calcScore();
}

Chromosome::~Chromosome()
{
}

vector<Chromosome> Chromosome::GATournamentSelection(vector<Chromosome> population, int tournamentSize, bool shufflePopulation) {
	if (shufflePopulation)
		random_shuffle(population.begin(), population.end());

	int nofTournaments = population.size() / tournamentSize;
	vector<Chromosome> selection;

	for (int i = 0; i < nofTournaments; i++){

		vector<Chromosome>::const_iterator start = population.begin() + i * tournamentSize;

		auto maxMember = max_element(start, start + tournamentSize, [](const Chromosome &lhs, const Chromosome &rhs)
		{
		  return lhs._score < rhs._score;
		});
		selection.push_back(*maxMember);
	}

	return selection;
}

vector<Chromosome> Chromosome::GAGenerateChildren(vector<Chromosome> parentsA, vector<Chromosome> parentsB){
	if (parentsA.size() != parentsB.size()) throw runtime_error("Parents don't have the same size");

	vector<Chromosome> children;
	children.reserve(parentsA.size());

	for (unsigned int i = 0; i < parentsA.size(); i++){
		children.push_back(GACrossOver(parentsA[i], parentsB[i]));
	}

	return children;
}

void Chromosome::invert() {
	for (int i = 0; i < _solution.size(); ++i) {
		_solution[i] ^= 1;

		_scoreContribution[i] = _nodeList[i]._links.size() - _scoreContribution[i];
	}
}

Chromosome Chromosome::GACrossOver(Chromosome parentA, Chromosome parentB) {
	int size = parentA._solution.size();
	int distance = 0;
	for (int j = 0; j < size; ++j) {
		distance += parentA._solution[j] != parentB._solution[j];
	}
	if (distance > size / 2) {
		parentA.invert();
	}

	vector<char> childSolution(size);
	int ones = 0, zeros = 0;
	for (int i = 0; i < size; ++i) {
		if (parentA._solution[i] == parentB._solution[i]) {
			childSolution[i] = parentA._solution[i];
			ones += parentA._solution[i];
			zeros += parentA._solution[i] ^ 1;
		}
		else {
			childSolution[i] = 2; // the placeholder
		}
	}

	random_device rd;
	mt19937 gen(rd());

	int addOnes = (size / 2) - ones;
	int addZeroes = (size / 2) - zeros;
	for (int k = 0; k < size; ++k) {
		if (childSolution[k] == 2) {
			bernoulli_distribution d(addOnes / (addOnes + addZeroes));
			char bit = d(gen);
			if (bit) {
				addOnes--;
			}
			else {
				addZeroes--;
			}
			childSolution[k] = bit;
		}
	}

	return Chromosome(childSolution);
}

ostream& operator<< (ostream & out, const Chromosome &sol){
	for (unsigned int i = 0; i < sol._solution.size(); i++){
		out << (int)sol._solution[i] << ",";
	}
	out << " " << sol._score;
	return out;
}
