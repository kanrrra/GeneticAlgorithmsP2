//#include "stdafx.h"
#include "Node.h"


Node::Node(string inputData)
{
	istringstream iss(inputData);

	int nofLinks, id;
	string location;

	iss >> id >> location >> nofLinks;

	_id = --id;

	for (auto badChar : _removeChars){
		location.erase(remove(location.begin(), location.end(), badChar), location.end());
	}

	auto commaPos = find(location.begin(), location.end(), ',');
	string x = string(location.begin(), commaPos);
	string y = string(commaPos + 1, location.end());

	//use of iss because of c++ bug making it impossible to use stof
	istringstream coordxss(x);
	coordxss >> _x;

	istringstream coordyss(y);
	coordyss >> _y;

	//_position = PointF(stof(x), stof(y));
	//_x = stof(x);
	//_y = stof(y);

	int link;
	for (int i = 0; i < nofLinks; i++){
		iss >> link;
		_links.push_back(--link);
	}
}


Node::~Node()
{
}
