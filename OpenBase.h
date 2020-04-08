#pragma once
#include"stdafx.h"
class OpenBase
{
public:
	OpenBase(Mat n);
	int * foobar(Mat n);
	~OpenBase();
private:
	Mat val;
};

