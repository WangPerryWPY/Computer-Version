/*实现canny边缘检测的原理通俗来说就是用离散化的梯度逼近函数根据二维灰度矩阵梯度向量来寻找图像灰度矩阵
的灰度跃变位置，然后再图像中将这些点连起来就形成了图像的边缘
一般对图像进行canny边缘检测要经过以下几步：
1.对图像进行灰度化处理：Gimg = 0.299R+0.587G+0.114B
2.对图像进行高斯滤波(出去图片中噪声对边缘检测的影响)
3.用一阶偏导的有限差分来计算梯度的幅值和方向
4.对梯度幅值进行非极大值抑制
5.双阈值检测和连接边缘
*/

#ifndef CANNY_H_
#define CANNY_H_
#include <iostream>
#include "CImg.h"
#include <cmath>

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
public:
	Canny();
	Canny(string name,string format);
	~Canny();
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
	//默认参数设置高斯滤波标准差为2.0，低阈值为0.25，高阈值为0.75
	CImg<int> canny_image();
	//整合所有获取最后的边缘图,sigma表示高斯滤波的参数，tlow和thigh为两个阈值
	CImg<int> canny_image(int sigma, float tlow, float thigh);
	//选出两个边缘点较近的距离连线
	CImg<int> canny_line(CImg<int> picture, int distance);
	//删掉长度小于20的边缘线
	CImg<int> delete_line(CImg<int> picture);
	//显示图像
	void display();
	//编写类过程的测试，无实际用处
	void test();
};

#endif