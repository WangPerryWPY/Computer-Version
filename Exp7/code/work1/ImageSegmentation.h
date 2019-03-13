#ifndef _ImageSegmentation_H_
#define _ImageSegmentation_H_
#include <iostream>
#include "CImg.h"
#include <string>
#include <math.h>
#include <vector>
#include <map>
#include <float.h>
using namespace cimg_library;
using namespace std;

#define PI 3.14159265358f
#define guassian_blur 2.5
#define Min_thres 250
#define DIFF 5

const int Red[] = {255,0,0};

class ImageSegmentation
{

private:
	//源图
	CImg<double> Img;
	//灰度图
	CImg<double> grayImg;
	//高斯平滑后的图像
	CImg<double> blurImg;
	//用于图像分割的阈值
	int threshold;
	//图像分割后的图片
	CImg<double> segImg;
	//梯度图像
	CImg<double> gradImg;
	//霍夫空间图像
	CImg<double> houghImg;
	CImg<double> resultImg;
	vector< pair< pair<int, int>, int > > peaks;
	//直线点集
	vector< pair< pair<int, int>, pair<int, int> > > lines;
	//四个角点
	vector< pair<int, int > > vertex;
	double M[9];

private:
	void rgb2gray();
	void Gauss_blur();
	//迭代法求阈值
	void get_thres_iteration();
	//OSTU法求阈值
	void get_thres_ostu();
	void get_thres(string type_c);
	void Segmentation();
	void gradDection();
	void Hough_Statistics();
	void GetLine();
	void GetVertexs();
	void orderVertexs();
	void calcMatrix();
	void warping();

public:
	ImageSegmentation();
	ImageSegmentation(const char* filename);
	~ImageSegmentation();
	void correct_process(string type_c);
};

#endif