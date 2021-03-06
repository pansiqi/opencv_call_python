// opencv_call_cpp.cpp: 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include "OpenBase.h"

extern "C"
{
	PRO_API void TestOpen(char * filename);
	PRO_API OpenBase * OpenBase_new(int rows, int cols, unsigned char* imgdata)
	{
		Mat img(rows, cols, CV_8UC3, (void *)imgdata);
		return new OpenBase(img);
	};
	PRO_API int * OpenBase_foobar(OpenBase * OpenBase, int rows, int cols, unsigned char* imgdata) 
	{
		Mat img(rows, cols, CV_8UC3, (void *)imgdata);
		return OpenBase -> foobar(img);
	}
	PRO_API double opencv_process(char * filename, int rho, double theta, int minLinLength, int maxLineGap, int area_scale);
}


void TestOpen(char * filename)
{
	cout << filename << endl;
}

double opencv_process(char * filename, int rho = 1, double theta = CV_PI / 180, int minLinLength = 50, int maxLineGap = 3, int area_scale = 2)
{
	Mat img;
	double angle_mean;
	img = OpenImage(filename);
	threshold(img);//图像二值化
	angle_mean = line_detection(img, rho, theta, minLinLength, maxLineGap);//线段检测
	if (ispower_scale)
	{
		scale_power(img, area_scale);//力的范围检测
	}
	return angle_mean;
}






