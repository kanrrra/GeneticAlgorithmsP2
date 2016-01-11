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

  	void invert();

	static vector<Chromosome> generateRandomPopulation(int populationSize, int solutionSize);
	static vector<Node> _nodeList;

  	static vector<Chromosome> GATournamentSelection(vector<Chromosome> population, int tournamentSize, bool shufflePopulation=true);
  	static vector<Chromosome> GAGenerateChildren(vector<Chromosome> parentsA, vector<Chromosome> parentsB);
  	static Chromosome GACrossOver(Chromosome parentA, Chromosome parentB);

	int swapNodesOpt();
	int mutate(double p);

	Chromosome(int size, bool optimal = false);
	Chromosome(vector<char> & solution);
	~Chromosome();


private:
	int calcScore();

	static vector<char> defaultDistribution;

};

