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


	if (false) {
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
	}

	if (false) {
		auto population = Chromosome::generateRandomPopulation(1000, nodes.size());

		int best = population[0]._score;
		for (auto pop : population){
			int oldScore = pop._score;
			int improvements = pop.swapNodesOpt();
			cout << oldScore << " " << improvements << " : " << pop._score << endl;

			if (pop._score < best) best = pop._score;
		}
		cout << "best: " << best << endl;
	}

	if (true) {
		auto population = Chromosome::generateRandomPopulation(50, nodes.size());

		// tornament selection of parrents

		do {

			auto parentsA = Chromosome::GATournamentSelection(population, 2);
			auto parentsB = Chromosome::GATournamentSelection(population, 2);

			auto children = Chromosome::GAGenerateChildren(parentsA, parentsB);
			// todo: local seachr on children
			// check if child is better than worst parent
			// add children to population
			// sort population
			// truncate population

		} while (false);

	}

	#ifdef _WIN64
		cin.ignore();
    #endif

	return 0;
}

