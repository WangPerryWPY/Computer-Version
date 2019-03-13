#ifndef _HOUGH_H_
#define _HOUGH_H_
#include "canny.h"
#include <vector>
#define PI 3.14159265358f
#define Min_thres 50
#define DIFF 5
const int Red[] = {255,0,0};
class hough {
public:
	hough();
	hough(const char *filename);
	hough(CImg<float> picture);
	~hough();
	vector< pair<int, int > > process();
private:
	CImg<float> Img;
	CImg<float> cannyImg;
	//霍夫空间图像
	CImg<float> houghImg;
	//画出A4纸直线后对应的图像
	CImg<float> resultImg;
	vector< pair< pair<int, int>, int > > peaks;
	//直线点集
	vector< pair< pair<int, int>, pair<int, int> > > lines;
	//四个角点
	vector< pair<int, int > > vertex;
private:
	//统计霍夫空间点
	void Hough_Statistics();
	//获取A4纸四条直线对应的点对
	void GetLine();
	//获取四个角点
	void GetVertexs();
	float distance (float , float);
};
#endif