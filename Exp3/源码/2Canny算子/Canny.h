#ifndef CANNY_H_
#define CANNY_H_
#include <iostream>
#include "CImg-2.4.0_pre090618\CImg.h"
#include <cmath>
#include <vector>
using namespace std;
using namespace cimg_library;

class Canny {
private:
	CImg<int> img;
	int rows;
	int cols;
	int *smoothedim;
	int *delta_x;  //x方向的一阶导数
	int *delta_y;  //y方向的一阶导数
	float *dirim;  //梯度的方向
	int *magnitude; //梯度的幅值
	int *nms;   //非极大值抑制后得到矩阵
	int *edge;  //边缘数组
	CImg<int> thres; //canny算子边缘图像
	CImg<int> HoughImg_Circle; //霍夫空间圆图像
	//CImg<int> CircleImg; //在霍夫圆图像上生成的像素点
	CImg<int> thres_img;
	vector<pair<int, int>> Circle;
	vector<int> Circleweight;
	vector<pair<int, int>> center;
public:
	Canny();
	~Canny();
	void HoughImageDetect(string name, float sigma, float tlow, float thigh, int Nums, int min_r, int max_r);
	void RGBtoGray();
	void gaussian_smooth(float sigma);
	//计算x,y方向的一阶导数
	void derrivative_x_y();
	//计算梯度向上的方向，以正x轴为逆时针方向指定的弧度
	void radian_direction(int xdirtag, int ydirtag);
	//计算梯度的幅值
	void magnitude_x_y();
	//进行非极大值抑制
	void non_max_supp();
	//双阈值检测
	void apply_hysteresis(float tlow, float thigh);
	//整合所有获取最后的边缘图,sigma表示高斯滤波的参数，tlow和thigh为两个阈值
	CImg<int> canny_image(int sigma, float tlow, float thigh);
	//选出两个边缘点较近的距离连线
	CImg<int> CircleDetect(int Nums, int min_r, int max_r);
	CImg<int> findpixelCircle(CImg<int> thres_img1, CImg<int> thres1);
	void CircleDetectImg(float sigma, float tlow, float thigh, int Nums, int min_r, int max_r);
};

#endif