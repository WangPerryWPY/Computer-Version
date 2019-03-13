#ifndef _IMAGESEGMENTATION_H
#define _IMAGESEGMENTATION_H

#include <list>
#include <fstream>
#include <sstream>
#include "warp.h"
#include <direct.h>

#define border_DIFF 8        // 图像边缘设为白色的距离

struct Point {
	int x, y;
	Point() : x(-1), y(-1) {}
	Point(int posX, int posY) : x(posX), y(posY) {}
};

class ImageSegmentation {
public:
	ImageSegmentation();
	ImageSegmentation(const char *);
	ImageSegmentation(CImg<float> );
	~ImageSegmentation();
	void process(const string, const char*);

private:
	CImg<float> warpImg;
	CImg<float> binaryImg, tagImg, histogramImg, dividingImg;
	vector< CImg<float> > subImageSet;     // 一行行数字图像
	vector<Point> dividePoints;          // 直方图峰值划分线点集
	int tagAccumulate;              	// 类别tag累加值
	vector<int> classTagSet;             // 类别tag列表
	vector< list<Point> > pointPosListSet; // 装载类别tag对应的所有点的位置的list的列表
	vector< list<Point> > pointPosListSetForDisplay;
	string basePath;                     // 单个数字图片生成、预测结果文本存放路径

private:
	// 图像二值化处理
	CImg<float> AdaptiveThreshold(CImg<float> warpingResult);
	// 做y方向的直方图，找到行与行之间的分割线
	void findDividingLine();
	// 通过分割线，将图片划分为一行行
	void divideIntoBarItemImg();
	//对每一张划分的图的数字，做扩张
	void toDilate(int barItemIndex);
	// 连通区域标记算法
	void connectedRegionsTagging(int barItemIndex);
	// 存储分割后每一张数字的图以及对应的文件名称
	void saveNum(int barItemIndex);
	// 添加新的类tag
	void addNewTag(int x, int y, int barItemIndex);
	// 在正上、左上、左中、左下这四个邻点中找到最小的tag
	void findMinTag(int x, int y, int &minTag, Point &minTagPointPos, int barItemIndex);
	// 合并某个点(x,y)所属类别
	void mergeTagImageAndList(int x, int y, const int minTag, const Point minTagPointPos, int barItemIndex);
	// 获取单个数字的包围盒
	void BoundingOfSingleNum(int listIndex, int& xMin, int& xMax, int& yMin, int& yMax);
	// 根据X方向直方图判断真实的拐点
	vector<int> getInflectionPosXs(const CImg<float>& XHistogramImage);
	// 获取一行行的子图的水平分割线
	vector<int> DivideLineXofSubImage(const CImg<float>& subImg);
	// XY方向的正扩张
	int Dilate(const CImg<float>& Img, int x, int y);
	// 分割行子图，得到列子图
	vector< CImg<float> > RowItemImg(const CImg<float>& lineImg, vector<int> _dividePosXset);
};

#endif