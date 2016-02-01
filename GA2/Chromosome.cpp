//#include "stdafx.h"
#include "Chromosome.h"
#include <stdexcept>
#include <numeric>
#include <math.h>

vector<char> Chromosome::defaultDistribution;
vector<Node> Chromosome::_nodeList;
//int Chromosome::gcrCalls = 0;

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

bool Chromosome::checkValidity() const {
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
	int oldScore = _score;


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

		flipNodeAtIdx(indices[i]);
		//_solution[indices[i]] ^= 1;
	}


	//_score = calcScore();

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

					//_score -= improvement;

					//swap nodes
					//_solution[i] ^= 1;
					//_solution[j] ^= 1;

					//cout << "improvement: " << improvement << endl;

					flipNodeAtIdx(i);
					flipNodeAtIdx(j);


					improvementsCount++;
					improvementFound = true;
					break;
				}
			}
		}


	} while (improvementFound);

	//int testScore = calcScore();
	//if (testScore != _score){
	//	throw runtime_error("bad score exception");
	//}

	return improvementsCount;
}

void Chromosome::flipNodeAtIdx(int idx){
//	std::clock_t ms_start = std::clock();
	char myColor = _solution[idx];
	int scoreChange = 0;
	for (int n : _nodeList[idx]._links){
		if (_solution[n] == myColor){
			_scoreContribution[n]++;
			scoreChange++;
		}
		else {
			_scoreContribution[n]--;
			scoreChange--;
		}
	}

	//cout << "idx: " << idx << " " << scoreChange << endl;

	_score += scoreChange;

	_scoreContribution[idx] = _nodeList[idx]._links.size() - _scoreContribution[idx];
	_solution[idx] ^= 1;
//	std::clock_t ms_end = std::clock();
//	Chromosome::_flipTime +=  1000 * (ms_end - ms_start) / CLOCKS_PER_SEC;
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

Chromosome::Chromosome(int size, Chromosome::GenerationType gt)
{
	switch (gt)
	{
	case Chromosome::OPTIMAL:
		_solution = generateOptimalSolutionForAssignmentData(size);
		break;
	case Chromosome::RANDOM:
		_solution = generateRandomSolution(size);
		break;
	case Chromosome::GREEDY:
	default:
		_solution = Chromosome::GRCsolution();
		break;
	}
	
	_score = calcScore();
}

Chromosome::Chromosome(int size)
{
//	_solution = Chromosome::GRCsolution();
	_solution = generateRandomSolution(size);
	_score = calcScore();
}

Chromosome::Chromosome(vector<char> & solution)
{
	_solution = solution;
	_score = calcScore();
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
	int distance = Chromosome::distance(parentA, parentB);
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
	out << " | valid: " << sol.checkValidity();
	return out;
}

Chromosome Chromosome::GRC() {
	auto solution = Chromosome::GRCsolution();
	return Chromosome(solution);
}

vector<char> Chromosome::GRCsolution(double badConnectionWeight) {

//	Chromosome::gcrCalls++;
	int size = _nodeList.size();
	vector<char> solution(size);
	int firstPosition = rand() % size;
	solution[firstPosition] = 1;

	//int selectedCount = 1;

	for (int selectedCount = 1; selectedCount < size / 2; selectedCount++){
		int toSelect = size - selectedCount;

		vector<Candidate> candidates(toSelect);

		int canId = 0; // candidate id for the array, may be faster than vectors
		for (int i = 0; i < size; ++i) {
			if (solution[i] == 0) {
				int connections = 0;
				int connectionsBad = 0;
				for (auto link : _nodeList[i]._links) {
					if (solution[link] == 1) {
						connections++;
					}
					else {
						connectionsBad++;
					}
				}

				candidates[canId]._connections = connections;
				candidates[canId]._connectionsBad = connectionsBad;
				candidates[canId]._score = (1-badConnectionWeight) * connections - badConnectionWeight * connectionsBad;
				candidates[canId]._id = i;

				canId++;
			}
		}

		sort(candidates.begin(), candidates.end(), [](const Candidate & a, const Candidate & b) {return a._score > b._score; });

//		for (auto c : candidates) {
//			cout << c._score << " ";
//		}
//		cout << endl;

		double lowestConnectionCount = candidates[0]._score;
		
//		if (partSize == 0) {
//			partSize = 1;
//		}

		int partSize;
		for (partSize = 1; partSize < toSelect; partSize++){
			if (lowestConnectionCount > candidates[partSize]._score){
				break;
			}
		}
		/*
		int partSize = 1;
		while (lowestConnectionCount == candidates[partSize]._score && partSize < toSelect) {
			partSize++;
		}*/
		int addId = rand() % partSize;



		solution[candidates[addId]._id] = 1;
//		cout << candidates[addId]._id  << " " << candidates[addId]._connections << endl;
	}
//
//	for (char c : solution) {
//		cout << (int) c;
//	}
//	cout << endl;

	return solution;
}
Chromosome::Chromosome(const Chromosome &b) {
	_score = b._score;
	_solution = b._solution;
	_scoreContribution = b._scoreContribution;
}

int Chromosome::distance(const Chromosome &a, const Chromosome &b) {
	int distance = 0;
	for (int j = 0; j < a._solution.size(); ++j) {
		distance += a._solution[j] != b._solution[j];
	}
	return distance;
}

Chromosome Chromosome::PathRelink(const Chromosome & a, const Chromosome & guidingSolution) {
	int distance = Chromosome::distance(a, guidingSolution);
	int size = a._solution.size();

	Chromosome current(a);

	Chromosome best(a);
	int scoreBest = best._score;

	Chromosome bestOverall(best);
	int scoreBestOverall = best._score;


	for (int k = 0; k < distance / 2; ++k) {
		scoreBest = numeric_limits<int>::max();

//		cout << current << endl;

		int i,j;
		for (i = 0; i < size; ++i) {
			if (current._solution[i] == 0 && current._solution[i] != guidingSolution._solution[i]) {
				Chromosome tmp(current);
				tmp.flipNodeAtIdx(i);
//				cout << "found 0 at " << i << endl;
				vector<int> candidates;
				for (j = 0; j < size; ++j) {
					if (i != j && tmp._solution[j] == 1 && tmp._solution[j] != guidingSolution._solution[j]) {
						Chromosome tmp2(tmp);
						tmp2.flipNodeAtIdx(j);

						if (tmp2._score < scoreBest) {
							scoreBest = tmp2._score;
							best = tmp2;
						}

						candidates.push_back(j);

//						cout << tmp2 << endl;
//						cout << "found 1 at " << j << endl;
					}
				}

//				int ra = rand() % (int)candidates.size();
//
//				Chromosome tmp2(tmp);
//				tmp2.flipNodeAtIdx(ra);
//
//				if (tmp2._score < scoreBest) {
//					scoreBest = tmp2._score;
//					best = tmp2;
//				}

			}
		}

		if (best._score < scoreBestOverall) {
			scoreBestOverall = best._score;
			bestOverall = best;
		}

		current = best;

//		cout << current << endl;
	}

//	cout << "valid: " << bestOverall.checkValidity() << endl;

//	if (scoreBestOverall > a._score && scoreBestOverall > guidingSolution._score) {
//		cout << "Better!!" << endl;
//	}
	return bestOverall;
}


Chromosome Chromosome::PathRelink2(const Chromosome & a, const Chromosome & guidingSolution) {
	int distance = Chromosome::distance(a, guidingSolution);
	int size = a._solution.size();

	Chromosome current(a);

	Chromosome best(a);
	int scoreBest = best._score;

	Chromosome bestOverall(best);
	int scoreBestOverall = best._score;

//	cout << distance << endl;

	int bit = 0;

	for (int k = 0; k < distance / 2; ++k) {
		scoreBest = numeric_limits<int>::max();

//		cout << current << endl;

		int i,j;
		for (i = 0; i < size; ++i) {
			if (current._solution[i] == bit && current._solution[i] != guidingSolution._solution[i]) {
				std::clock_t ms_start = std::clock();
				Chromosome tmp(current);
				tmp.flipNodeAtIdx(i);
				std::clock_t ms_end = std::clock();
				Chromosome::_flipTime +=  ms_end - ms_start;


				if (tmp._score < scoreBest) {
					scoreBest = tmp._score;
					best = tmp;
				}

			}
		}

		bit ^= 1;

//		cout << "bit: " << bit << endl;

		if (bit == 0 && best._score < scoreBestOverall) {
			scoreBestOverall = best._score;
			bestOverall = best;
		}

		current = best;

//		cout << current << endl;
	}

//	cout << "valid: " << bestOverall.checkValidity() << endl;
//
//	if (scoreBestOverall > a._score && scoreBestOverall > guidingSolution._score) {
//		cout << "Better!!" << endl;
//	}
	return bestOverall;
}

Chromosome Chromosome::PathRelink3(const Chromosome & a, const Chromosome & guidingSolution) {
	int distance = Chromosome::distance(a, guidingSolution);
	int size = a._solution.size();

	Chromosome current(a);
	Chromosome best(a);
	int scoreBest = best._score;

	Chromosome bestOverall(best);
	int scoreBestOverall = best._score;

//	cout << distance << endl;

	int bit = 0;

	for (int k = 0; k < distance / 2; ++k) {
		scoreBest = numeric_limits<int>::max();

//		cout << current << endl;

		int i,j;
		for (i = 0; i < size; ++i) {
			if (current._solution[i] == bit && current._solution[i] != guidingSolution._solution[i]) {
//				std::clock_t ms_start = std::clock();
				int tmpScore = current.FakeFlipNodeAtIdx(i);
				if (tmpScore < scoreBest) {
					best = current;
					best.flipNodeAtIdx(i);
					scoreBest = best._score;
					if (tmpScore != best._score) {
						cout << "error" << endl;
					}
				}
//				std::clock_t ms_end = std::clock();
//				Chromosome::_flipTime +=  ms_end - ms_start;
			}
		}

		bit ^= 1;

		if (bit == 0 && best._score < scoreBestOverall) {
			scoreBestOverall = best._score;
			bestOverall = best;
		}

		current = best;

//		cout << current << endl;
	}

	return bestOverall;
}


//void Chromosome::flipNodeAtIdx(int idx){
//	char myColor = _solution[idx];
//	int scoreChange = 0;
//	for (int n : _nodeList[idx]._links){
//		if (_solution[n] == myColor){
//			_scoreContribution[n]++;
//			scoreChange++;
//		}
//		else {
//			_scoreContribution[n]--;
//			scoreChange--;
//		}
//	}
//
//	//cout << "idx: " << idx << " " << scoreChange << endl;
//
//	_score += scoreChange;
//
//	_scoreContribution[idx] = _nodeList[idx]._links.size() - _scoreContribution[idx];
//	_solution[idx] ^= 1;
//}
//int Chromosome::scoreChange(vector<char> &solution, int idx) {
//	return 0;
//}
Chromosome::~Chromosome() {

}

int Chromosome::FakeFlipNodeAtIdxScoreChange(int idx) {
	char myColor = _solution[idx];
	int scoreChange = 0;

	for (int n : _nodeList[idx]._links){
		if (_solution[n] == myColor){
			scoreChange++;
		}
		else {
			scoreChange--;
		}
	}
	return scoreChange;
}

int Chromosome::FakeFlipNodeAtIdx(int idx) {
	return _score + FakeFlipNodeAtIdxScoreChange(idx);
}

int Chromosome::FakeFlipNodeAtIdx(int idx, int flipped, int scoreChange) {
	char myColor = _solution[idx];
	char flippedColor = _solution[flipped];

//	cout << "flipped " << alreadyFlipped << endl;

	for (int n : _nodeList[idx]._links){
		if (_solution[n] == myColor){
			if (flipped == n) {
				scoreChange--;
			}
			else {
				scoreChange++;
			}
		}
		else {
			if (flipped == n) {
				scoreChange++;
			}
			else {
				scoreChange--;
			}
		}
	}

	return _score + scoreChange;
}


Chromosome Chromosome::PathRelink4(const Chromosome & a, const Chromosome & guidingSolution) {
	int distance = Chromosome::distance(a, guidingSolution);
	int size = a._solution.size();

	Chromosome current(a);

	Chromosome best(a);
	int scoreBest = best._score;

	Chromosome bestOverall(best);
	int scoreBestOverall = best._score;


	for (int k = 0; k < distance / 2; ++k) {
		scoreBest = numeric_limits<int>::max();

//		cout << current << endl;

		int i,j;
		for (i = 0; i < size; ++i) {
			if (current._solution[i] == 0 && current._solution[i] != guidingSolution._solution[i]) {
				int scoreChange = current.FakeFlipNodeAtIdxScoreChange(i);
//				Chromosome tmp(current);
//				tmp.flipNodeAtIdx(i);
//				cout << "found 0 at " << i << endl;
//				vector<int> candidates;

				std::clock_t ms_start = std::clock();

				for (j = 0; j < size; ++j) {
					if (i != j && current._solution[j] == 1 && current._solution[j] != guidingSolution._solution[j]) {

						int tmpScore = current.FakeFlipNodeAtIdx(j, i, scoreChange);

						if (tmpScore < scoreBest) {
							best = current;
							best.flipNodeAtIdx(i);
							best.flipNodeAtIdx(j);
							scoreBest = best._score;
//							if (tmpScore != best._score) {
//								cout << "error " << tmpScore << " " << best._score << endl;
//							}
						}

//						cout << tmp2 << endl;
//						cout << "found 1 at " << j << endl;
					}
				}

				std::clock_t ms_end = std::clock();
				Chromosome::_flipTime +=  ms_end - ms_start;

//				int ra = rand() % (int)candidates.size();
//
//				Chromosome tmp2(tmp);
//				tmp2.flipNodeAtIdx(ra);
//
//				if (tmp2._score < scoreBest) {
//					scoreBest = tmp2._score;
//					best = tmp2;
//				}

			}
		}

		if (best._score < scoreBestOverall) {
			scoreBestOverall = best._score;
			bestOverall = best;
		}

		current = best;

//		cout << current << endl;
	}

//	cout << "valid: " << bestOverall.checkValidity() << endl;

//	if (scoreBestOverall > a._score && scoreBestOverall > guidingSolution._score) {
//		cout << "Better!!" << endl;
//	}
	return bestOverall;
}