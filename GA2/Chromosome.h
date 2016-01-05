#pragma once

#include <vector>
#include <random>

#include <algorithm>
#include "Node.h"

using namespace std;

class Chromosome
{
public:
	int _score;
	vector<char> _solution;

	vector<char> generateRandomSolution(int size);

	static vector<Chromosome> GenerateRandomPopulation(int populationSize, int solutionSize);
	static vector<Node> _nodeList;

	int swapNodesOpt();

	Chromosome(int size);
	~Chromosome();

private:
	void updateScore();

	static vector<char> defaultDistribution;

};

