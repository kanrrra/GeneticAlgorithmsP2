// GA2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>

#include "DataReader.h"
#include "Chromosome.h"

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	cout << "hello world" << endl;

	vector<Node> nodes = DataReader::GetData("data.txt");
	Chromosome::_nodeList = nodes;

	auto pop = Chromosome(nodes.size(), true);
	cout << pop._score << " isValid: " << pop.checkValidity() << endl;

	auto population = Chromosome::GenerateRandomPopulation(1000, nodes.size());
	
	int best = population[0]._score;
	for (auto pop : population){
		int oldScore = pop._score;
		int improvements = pop.swapNodesOpt();
		cout << oldScore << " " << improvements << " : " << pop._score << endl;

		if (pop._score < best) best = pop._score;
	}
	cout << "best: " << best << endl;
	 
	cin.ignore();
	return 0;
}

