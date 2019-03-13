/*
网格的生成
将源图像和目标图像通过建立特征点，形成点集，然后通过三角剖分的方法形成三角形网络。
1）把源图像中形成的三角形和目标图像生成的三角形（三角形对）对应起来。（对整体的源图像和目标图像来说，他们划分出网格之后，应该是同构的。）
2）通过源三角形和目标三角形的顶点坐标值，使用仿射变换求出变换从源三角形到目标三角形的变换矩阵T。
*/
/*
中间帧生成
对其中一个三角形对来说,中间帧的生成过程是这样的：
1) 通过变换矩阵T，求得三角形对的对应像素点坐标。
2) 定位源三角形内部像素点P0的RGB值，经过线型插值运算：Pinternal=(1-1/n)P0+(1/n)P1（Pinternal是中间帧像素点RGB值，P1是目标像素点RGB值，n为变形动画的帧数）获得中间帧中点Pinternal的RGB值。
通过以上方法，求得其余三角形的中间帧点Pinternal的RGB值，并将他们写入中间帧缓存中，最终生成中间帧图像。
*/
#include <iostream>
#include "CImg.h"
#include <vector>
using namespace std;
using namespace cimg_library;

#define frame 11 //图片变换的帧数

//像素点
struct point {
	int x, y;
	point(int x0, int y0) :x(x0), y(y0) {}
};

//网格三角
struct triangle {
	point a, b, c;
	triangle(point a0, point b0, point c0) :a(a0), b(b0), c(c0) {}
};

class ImageMorphing {
private:
	CImg<float> source; //源图
	CImg<float> target; //目标图
	CImgList<float> result; //结果图片集
	vector<point> source_point; //源图网格点
	vector<point> target_point; //目标图网格点
	vector<triangle> source_triangle; //源图网格
	vector<triangle> target_triangle; //目标图网格
	vector<vector<triangle>> mid_triangle; //中间帧网格
	vector<vector<int>> triangle_grid; //需要参与变换的网格
	vector<vector<CImg<float>>> source_matrix; //源图到中间帧的变换矩阵
	vector<vector<CImg<float>>> target_matrix; //目标图到中间帧的变换矩阵
public:
	ImageMorphing();
	void get_middleGrid();
	void AffineTransform();
	void Morphing_process();
};