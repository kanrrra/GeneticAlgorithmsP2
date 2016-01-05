#include "stdafx.h"
#include "Chromosome.h"

vector<char> Chromosome::defaultDistribution;
vector<Node> Chromosome::_nodeList;

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
					if (improvement > 0){
						char tmp = _solution[i];
						_solution[i] = _solution[j];
						_solution[j] = tmp;

						_score -= improvement;
						improvementsCount++;
						break;
					}
				}
			}
		}


	} while (improvementFound);

	return improvementsCount;
}

void Chromosome::updateScore(){
	_score = 0;
	for (Node n : _nodeList){
		char myColor = _solution[n._id];
		for (int meighbourId : n._links){
			_score += (myColor != _solution[meighbourId]);
		}
	}
}

//Chromosome::Chromosome(Chromosome c){
//
//}

Chromosome::Chromosome(int size)
{
 	_solution = generateRandomSolution(size);
	updateScore();
}


Chromosome::~Chromosome()
{
}
