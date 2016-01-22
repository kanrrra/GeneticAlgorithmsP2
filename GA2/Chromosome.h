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
	enum GenerationType {
		OPTIMAL, RANDOM, GREEDY
	};

	int _score;
	vector<char> _solution;
	vector<char> _scoreContribution;

	vector<char> generateRandomSolution(int size);
	vector<char> generateOptimalSolutionForAssignmentData(int size);
	bool checkValidity();

  	void invert();

	static vector<Chromosome> generateRandomPopulation(int populationSize, int solutionSize);
	static vector<Node> _nodeList;

  	static vector<Chromosome> GATournamentSelection(vector<Chromosome> population, int tournamentSize, bool shufflePopulation=true);
  	static vector<Chromosome> GAGenerateChildren(vector<Chromosome> parentsA, vector<Chromosome> parentsB);
  	static Chromosome GACrossOver(Chromosome parentA, Chromosome parentB);
  	static Chromosome PathRelink(const Chromosome & a, const Chromosome & b);

	void flipNodeAtIdx(int idx);
	int swapNodesOpt();
	int mutate(unsigned int perturbationSize);

  	Chromosome(int size);
	Chromosome(int size, Chromosome::GenerationType gt);
	Chromosome(const Chromosome & b);
	Chromosome(vector<char> & solution);

	friend ostream& operator<< (ostream& out, const Chromosome& sol);

  	static Chromosome GRC();
  	static vector<char> GRCsolution(double badConnectionWeight = 0);
//  	static int gcrCalls;

	static int distance(const Chromosome & a, const Chromosome & b);
private:
	int calcScore();
  	static int scoreChange(vector<char> & solution, int idx);

	static vector<char> defaultDistribution;

  	struct Candidate {
			int _id;
			int _connections;
			int _connectionsBad;
			double _score;
	};

};

