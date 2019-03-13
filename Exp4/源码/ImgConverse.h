#include <iostream>
#include "CImg.h"
using namespace std;
using namespace cimg_library;
class ImgConverse{
private:
	CImg<int> img; //图像
	CImg<int> gray_img; //灰度图像
	CImg<int> hist_equalImg; //直方图均匀化后的图片
	CImg<float> sourceImg; //进行图像转换的原图
	CImg<float> targetImg; //进行图像转换的目标背景图
	CImg<float> sourcelab; //lab空间图
	CImg<float> targetlab;
	//用于转换计算的均值
	float mean0[3];
	float mean1[3];
	//用于转换计算的标准差
	float std0[3];
	float std1[3];
	CImg<float> result; //转换后的结果图
public:
	ImgConverse();
	void toGray(); //将图片处理为灰度图
	CImg<int> Hist(CImg<int> ); //进行直方图均衡的主体函数
	void Hist_gray(); //对灰度图进行直方图均衡处理
	//对rgb图进行直方图均衡处理，以下是两种方法
	void Hist_color();  //1.对rgb空间分别进行直方图均衡
	void RGB_HSI();  //2.将图片由rgb空间转为hsi空间并对其亮度空间进行直方图均衡再转为rgb进行显示
	CImg<float> RGB_LAB(CImg<float>, float *, float *);//将图片转为lab空间
	CImg<float> LAB_RGB(CImg<float>);
	void colorTransfer(); //对lab空间进行转换处理
	void Hist_impress(); //对灰度处理用梯度求边缘进行改进
};