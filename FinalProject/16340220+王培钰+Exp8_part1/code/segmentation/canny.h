#ifndef _CANNY_H_
#define _CANNY_H_
#include <math.h>
#include "function.h"
class canny {
private:
	//源图
	CImg<float> Img;
	//灰度图
	CImg<float> grayImg;
	//高斯平滑后的图像
	CImg<float> blurImg;
	//sobel算子分割出的边缘(梯度的幅值)
	CImg<float> sobelImg;
	//梯度的方向
	CImg<float> angleImg;
	//非极大值抑制后图像
	CImg<float> nonImg;
	//双阈值分割后的图像
	CImg<float> thresImg; 
private:
	//高斯滤波
	void Gauss_blur(float);
	//sobel算子平滑
	void sobel();
	//非极大值抑制
  	void nonMaxSupp();
  	//双阈值检测
  	void threshold(int, int); 
public:
	canny();
	canny(const char *filename);
	canny(CImg<float> picture);
	~canny();
	CImg<float> process(float, int, int);
};
#endif