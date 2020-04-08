#include "stdafx.h"
int  postion_x = 0;
int  postion_y = 0;
double area = 0;
bool ispower_scale = false;
// TODO: 在 STDAFX.H 中引用任何所需的附加头文件，
//而不是在此文件中引用
Mat calcGrayHist(const Mat  & image)

{
	Mat histogram = Mat::zeros(Size(256, 1), CV_32SC1);
	//图像的高和宽
	int rows = image.rows;
	int cols = image.cols;
	for (int r = 0; r < rows; r++)
	{
		for (int c = 0; c < cols; c++)
		{
			int index = int(image.at<uchar>(r, c));
			histogram.at<int>(0, index) += 1;
		}
	}
	return histogram;
}

Mat equalHist(Mat image)
{
	//CV_Assert(image.type() == CV_8UC1);
	int rows = image.rows;
	int cols = image.cols;
	//计算图像的灰度直方图
	Mat grayHist = calcGrayHist(image);
	Mat zeroCumuMoment = Mat::zeros(Size(256, 1), CV_32SC1);
	//计算累加直方图
	for (int p = 0; p < 256; p++)
	{
		if (p == 0)
			zeroCumuMoment.at<int>(0, p) = grayHist.at<int>(0, 0);
		else
			zeroCumuMoment.at<int>(0, p) = zeroCumuMoment.at<int>(0, p - 1) + grayHist.at<int>(0, p);
	}
	//根据累加直方图得到输入灰度图和输出灰度图之间的映射关系
	Mat outPut_q = Mat::zeros(Size(256, 1), CV_8UC1);
	float cofficient = 256.0 / (rows*cols);
	for (int p = 0; p < 256; p++)
	{
		float q = cofficient * zeroCumuMoment.at<int>(0, p) - 1;
		if (q >= 0)
			outPut_q.at<uchar>(0, p) = uchar(floor(q));
		else
			outPut_q.at<uchar>(0, p) = 0;
	}
	//得到直方图均衡化的图像
	Mat equalHistImage = Mat::zeros(image.size(), CV_8UC1);
	for (int r = 0; r < rows; r++)
	{
		for (int c = 0; c < cols; c++)
		{
			int p = image.at<uchar>(r, c);
			equalHistImage.at<uchar>(r, c) = outPut_q.at<uchar>(0, p);

		}
	}
	return  equalHistImage;
}

//same 卷积

void conv2D(InputArray src, InputArray kernel, OutputArray dst, int ddepth, Point anchor = Point(-1, -1), int borderType = BORDER_DEFAULT)
{
	Mat kernelFlip;
	flip(kernel, kernelFlip, -1);
	filter2D(src, dst, ddepth, kernelFlip, anchor, 0.0, borderType);
}
//可分离的same卷积  先垂直卷积，然后进行水平方向的卷积
void sepConv2D_Y_X(InputArray src, OutputArray src_kerY_kerX, int ddepth, InputArray kernelY, InputArray kernelX, Point anchor = Point(-1, -1), int borderType = BORDER_DEFAULT)
{
	Mat src_kerY;
	conv2D(src, kernelY, src_kerY, ddepth, anchor, borderType);
	conv2D(src_kerY, kernelX, src_kerY_kerX, ddepth, anchor, borderType);

}

void sepConv2D_X_Y(InputArray src, OutputArray src_kerX_kerY, int ddepth, InputArray kernelX, InputArray kernelY, Point anchor = Point(-1, -1), int borderType = BORDER_DEFAULT)
{
	//输入矩阵与水平方向卷积核的卷积
	Mat src_kerX;
	conv2D(src, kernelX, src_kerX, ddepth, anchor, borderType);
	//上面得到的卷积结果，然后接着和垂直方向的卷积核卷积，得到最终的输出
	conv2D(src_kerX, kernelY, src_kerX_kerY, ddepth, anchor, borderType);
}

void prewitt(InputArray src, OutputArray dst, int ddepth, int x, int y = 0, int borderType = BORDER_DEFAULT)
{
	CV_Assert(!(x == 0 && y == 0));
	//如果 x!=0，src 和 prewitt_x卷积核进行卷积运算
	if (x != 0)
	{
		//可分离的prewitt_x卷积核
		Mat prewitt_x_y = (Mat_<float>(3, 1) << 1, 1, 1);
		Mat prewitt_x_x = (Mat_<float>(1, 3) << 1, 0, -1);
		//可分离的离散的二维卷积
		sepConv2D_Y_X(src, dst, ddepth, prewitt_x_y, prewitt_x_x, Point(-1, -1), borderType);
	}
	if (y != 0)
	{
		//可分离的prewitt_y卷积核
		Mat prewitt_y_x = (Mat_<float>(1, 3) << 1, 1, 1);
		Mat prewitt_y_y = (Mat_<float>(3, 1) << 1, 0, -1);
		//可分离的离散二维卷积
		sepConv2D_X_Y(src, dst, ddepth, prewitt_y_x, prewitt_y_y, Point(-1, -1), borderType);
	}
}
void callback_thresh(int, void*)
{

}

//高斯平滑

Mat gaussBlur(const Mat & image, Size winSize, float sigma, int ddepth = CV_64F, Point anchor = Point(-1, 1), int borderType = BORDER_DEFAULT)
{
	CV_Assert(winSize.width % 2 == 1 && winSize.height % 2 == 1);
	//构建垂直方向的高斯卷积算子
	Mat gK_y = getGaussianKernel(sigma, winSize.height, CV_64F);
	Mat gK_x = getGaussianKernel(sigma, winSize.width, CV_64F);
	gK_x = gK_x.t();//转置
					//分离高斯卷积
	Mat blurImage;
	sepConv2D_Y_X(image, blurImage, ddepth, gK_y, gK_x, Point(-1, 1));
	return blurImage;
}

void roberts(InputArray src, OutputArray dst, int ddepth, int x, int y = 0, int borderType = BORDER_DEFAULT)
{
	CV_Assert(!(x == 0 && y == 0));
	Mat roberts_1 = (Mat_<float>(2, 2) << 1, 0, 0, -1);
	Mat roberts_2 = (Mat_<float>(2, 2) << 0, 1, -1, 0);
	//当 x 不等于零时，src 和 roberts_1 卷积
	if (x != 0)
	{
		conv2D(src, roberts_1, dst, ddepth, Point(0, 0), borderType);
	}
	//当 y 不等于零时，src 和 roberts_2 卷积
	if (y != 0)
	{
		conv2D(src, roberts_2, dst, ddepth, Point(0, 0), borderType);
	}
}

int factorial(int n)
{
	int fac = 1;
	if (n == 0)
	{
		return fac;
	}
	for (int i = 1; i <= n; i++)
		fac *= i;
	return fac;
}
//计算平滑系数
Mat getPascalSmooth(int n)
{
	Mat pascalSmooth = Mat::zeros(Size(n, 1), CV_32FC1);
	for (int i = 0; i < n; i++)
		pascalSmooth.at<float>(0, i) = factorial(n - 1) / (factorial(i) * factorial(n - 1 - i));
	return pascalSmooth;
}
//计算差分
Mat getPascalDiff(int n)
{
	Mat pascalDiff = Mat::zeros(Size(n, 1), CV_32FC1);
	Mat pascalSmooth_previous = getPascalSmooth(n - 1);
	for (int i = 0; i < n; i++)
	{
		if (i == 0)
			pascalDiff.at<float>(0, i) = 1;
		else if (i == n - 1)
			pascalDiff.at<float>(0, i) = -1;
		else
			pascalDiff.at<float>(0, i) = pascalSmooth_previous.at<float>(0, i) - pascalSmooth_previous.at<float>(0, i - 1);
	}
	return pascalDiff;
}


//求两点之间的角度
double get_point_angle(CvPoint pointO, CvPoint pointA)
{
	double angle = 0;
	CvPoint point;
	double temp;

	point = cvPoint((pointA.x - pointO.x), (pointA.y - pointO.y));

	if ((0 == point.x) && (0 == point.y))
	{
		return 0;
	}

	if (0 == point.x)
	{
		angle = 90;
		return angle;
	}

	if (0 == point.y)
	{
		angle = 0;
		return angle;
	}

	temp = fabsf(float(point.y) / float(point.x));
	temp = atan(temp);
	temp = temp * 180 / CV_PI;

	if ((0 < point.x) && (0 < point.y))
	{
		angle = 180 - temp;
		return angle;
	}

	if ((0 > point.x) && (0 < point.y))
	{
		angle = 360 - (180 - temp);
		return angle;
	}

	if ((0 < point.x) && (0 > point.y))
	{
		angle = temp;
		return angle;
	}

	if ((0 > point.x) && (0 > point.y))
	{
		angle = 180 - temp;
		return angle;
	}

	printf("sceneDrawing :: getAngle error!");
	return -1;
}

Mat OpenImage(char *  filename)
{
	cout << filename << endl;
	Mat img = imread(filename, IMREAD_GRAYSCALE);
	cout << "Open image:" << img  << endl;
	cout << "Calling Successful " << endl;
	return img;
}
//线段检测
double line_detection(Mat img,int rho,double theta,int minLinLength, int maxLineGap)
{
	/*
	img:打开的文件
	image: 边缘检测的输出图像. 它应该是个灰度图 (但事实上是个二值化图)
	rho : 　参数极径  以像素值为单位的分辨率. 我们使用 1 像素.
	theta: 参数极角  以弧度为单位的分辨率. 我们使用 1度 (即CV_PI/180)
	threshold: 要”检测” 一条直线所需最少的的曲线交点 
	minLinLength: 能组成一条直线的最少点的数量. 点数量不足的直线将被抛弃.线段的 最小长度
	maxLineGap:线段上最近两点之间的阈值
	*/
	vector<Vec4i> lines;
	HoughLinesP(img, lines, rho, theta, minLinLength, minLinLength, maxLineGap);
	cout << lines.size() << endl;
	double sum = 0;
	double angle_mean = 0;
	for (size_t i = 0; i < lines.size(); i++)
	{
		Vec4i I = lines[i];
		double angle = get_point_angle(Point(I[0], I[1]), Point(I[2], I[3]));
		cout << "angle:" << angle << endl;//输出角度
		sum = sum + angle;
	}
	if (lines.size() != 0)
	{
		angle_mean = sum / lines.size();
		cout << "the mean angle is " << angle_mean << endl;
		ispower_scale = true;
		return angle_mean;
	}
	else
		cout << "Point" << endl;
}


//计算力的范围
double  scale_power(Mat img,int area_scale)
{
	/*
	img :打开的图像
	area_scale :设定力的范畴
	*/
	vector<vector<Point>> contours;
	findContours(img, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
	int num = contours.size();
	for (int i = 0; i < num; i++)
	{
		Rect rect = boundingRect(contours[i]);
		rect.x;
		if (rect.area() > area_scale)
		{
			cout << "power scale:" <<rect.area() << endl;
			area = rect.area();
			/*
			TODO:当力的拟合范围不只是一个时，采用最大的那个
			*/
			center_point(rect);
		}
	}
	ispower_scale = false;
	return area;
}

//计算力的中心点
/*
 ____
|\   |
| \  |          
|  \ |
|___\|

*/
void center_point(Rect rect)
{
	int postion[2];
	int  x,y;
	double x0, y0;
	x0 = (rect.br().x - rect.tl().x)/2;
	y0 = (rect.tl().y - rect.br().y)/2;
	x = rect.tl().x + int(x0);
	y = rect.br().y - int(y0);
	postion_x = x;
	postion_y = y;
}

//图像二值化

void threshold(Mat &img)
{
	double Imax, Imin;
	minMaxLoc(img, &Imin,&Imax,NULL,NULL);
	int mid = (Imax - Imin) / 2;
	threshold(img, img, mid, Imax, CV_THRESH_BINARY);
}

/*
TODO:将全局变量传递给函数，导出函数直接返回
*/




