// GA2.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include <iostream>

#include "DataReader.h"
#include "Chromosome.h"

using namespace std;

int main(int argc, char* argv[])
{

	vector<Node> nodes = DataReader::GetData("data.txt");
	Chromosome::_nodeList = nodes;

	auto pop = Chromosome(nodes.size(), true);
	cout << "optimal soution of provided graph: " << pop._score << " isValid: " << pop.checkValidity() << endl;


	for (int i = 0; i < 1; i++){
		Chromosome candidate = Chromosome(nodes.size());
		candidate.swapNodesOpt();

		cout << "starting score: " << candidate._score << endl;
		for (int attempts = 0; attempts < 10; attempts++){
			candidate.mutate(0.1);
			cout << "after mutation: " << candidate._score << endl;
			candidate.swapNodesOpt();
			cout << "result score: " << candidate._score << " valid? " << candidate.checkValidity() << endl;
		}

	}


	

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

