#include "stdafx.h"
#include "OpenBase.h"


OpenBase::OpenBase(Mat n)
{
	val = n;
	cout << "Value is" << val << endl;
}

OpenBase::~OpenBase()
{
	cout << "Object is already destroyed" << endl;
}

int * OpenBase::foobar(Mat n)
{
	int * data = new int[2];
	data[0] = 1;
	data[1] = 2;
	return data;
}

