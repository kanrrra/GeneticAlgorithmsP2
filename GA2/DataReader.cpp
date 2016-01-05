#include "stdafx.h"
#include "DataReader.h"


vector<Node> DataReader::GetData(string path){
	ifstream datafile(path);

	vector<Node> nodes;

	string line;
	while (getline(datafile, line)){
		Node currentNode(line);
		nodes.push_back(currentNode);
	}

	return nodes;
}


DataReader::~DataReader()
{
}
