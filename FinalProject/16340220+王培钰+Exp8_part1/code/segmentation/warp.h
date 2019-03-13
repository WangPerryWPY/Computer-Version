#ifndef _WARP_H_
#define _WARP_H_
#include "hough.h"
#include <fstream>
#include <string>
#include <sstream>
class warp
{
public:
	warp();
	warp(const char *);
	warp(CImg<float> );
	~warp();
	CImg<float> process(const char *txtname);

private:
	CImg<float> Img;
	//四个角点
	vector< pair<int, int > > vertex;
	//用于A4纸矫正的特征矩阵
	float M[9];
	CImg<float> resultImg;
	int paper_width;
	int paper_height;

private:
	void orderVertexs();
	void calcMatrix();
	void getSize();
	void warping();
	float distance (float , float);
};
#endif