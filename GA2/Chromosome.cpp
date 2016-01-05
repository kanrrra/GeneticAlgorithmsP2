#include "stdafx.h"
#include "Chromosome.h"

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
		else if (sol != 0) throw new exception("Invalid value exception");
	}

	cout << count << endl;

	return (count * 2) == _solution.size();
}

vector<Chromosome> Chromosome::GenerateRandomPopulation(int populationSize, int solutionSize){
	vector<Chromosome> population;
	population.reserve(populationSize);

	for (int i = 0; i < populationSize; i++){
		population.push_back(Chromosome(solutionSize));
	}

	return population;
}

int Chromosome::swapNodesOpt(){
	int improvementsCount = 0;
	bool improvementFound = false;
	do {
		improvementFound = false;

		for (int i = 0; i < _solution.size() - 1; i++){

			int localScoreI = 0;
			for (auto nb : _nodeList[i]._links){
				localScoreI += (_solution[i] != _solution[nb]);
			}
			localScoreI = 2 * localScoreI - _nodeList[i]._links.size();

			for (int j = i + 1; j < _solution.size(); j++){
				if (_solution[i] != _solution[j]){
					
					int localScoreJ = 0;
					for (auto nb : _nodeList[j]._links){
						localScoreJ += (_solution[j] != _solution[nb]);
					}
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

					_solution[i] ^= 1;
					_solution[j] ^= 1;

					_score -= improvement;

					improvementsCount++;
					break;
				}
			}
		}


	} while (improvementFound);

	int testScore = calcScore();
	if (testScore != _score){
		throw new exception("bad score exception");
	}

	return improvementsCount;
}

int Chromosome::calcScore(){
	int score = 0;
	for (Node n : _nodeList){
		char myColor = _solution[n._id];
		for (int meighbourId : n._links){
			score += (myColor != _solution[meighbourId]);
		}
	}

	if (score % 2) throw new exception("Score cannot be odd before division");

	score /= 2;

	return score;
}

//Chromosome::Chromosome(Chromosome c){
//
//}

Chromosome::Chromosome(int size, bool optimal)
{
	if (optimal) _solution = generateOptimalSolutionForAssignmentData(size);
 	else _solution = generateRandomSolution(size);
	
	_score = calcScore();
}


Chromosome::~Chromosome()
{
}
