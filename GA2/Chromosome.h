#pragma once

#include <vector>
#include <random>

#include <algorithm>
#include "Node.h"

//todo remove useless include
#include <iostream>

using namespace std;

class Chromosome
{
public:
	int _score;
	vector<char> _solution;

	vector<char> generateRandomSolution(int size);
	vector<char> generateOptimalSolutionForAssignmentData(int size);
	bool checkValidity();

	static vector<Chromosome> GenerateRandomPopulation(int populationSize, int solutionSize);
	static vector<Node> _nodeList;

	int swapNodesOpt();
	int mutate(double p);

	Chromosome(int size, bool optimal = false);
	~Chromosome();

private:
	int calcScore();

	static vector<char> defaultDistribution;

};

