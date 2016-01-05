#pragma once

#include <string>
#include <sstream>

#include <vector>
#include <algorithm>

//#include <Windows.h>
//#include <Gdiplus.h>


using namespace std;

//using namespace Gdiplus;
#pragma comment (lib,"Gdiplus.lib")

class Node
{
public:
	Node(string inputData);
	~Node();

	int _id;
	vector<int> _links;

	//Gdiplus::PointF _position;
	float _x, _y;

private:
	string _removeChars = "()";

};

