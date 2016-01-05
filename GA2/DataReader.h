#pragma once

#include <string>
#include <fstream>

#include "Node.h"
#include <vector>

using namespace std;

class DataReader
{
public:
	static vector<Node> GetData(string path);

	~DataReader();

private: 

	DataReader();


};

