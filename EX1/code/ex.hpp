#ifndef _EX_HPP_
#define _EX_HPP_
#include "../CImg.h"
#include <cmath>
#include <string>
using namespace std;
using namespace cimg_library;
const double pi(3.14159265);
class Test
{
	public:
		Test();
		~Test();
		void Todisplay();
		void change(); //把1.bmp文件的白色区域变成红色，黑色区域变成绿色
		void DrawCircle_blue1(); //不使用CImg函数在图上绘制一个圆形区域，圆心坐标(50,50)，半径为30，填充颜色为蓝色
		void DrawCircle_yellow1();//不使用CImg函数在图上绘制一个圆形区域，圆心坐标(50,50)，半径为3，填充颜色为黄色
		void DrawLine1();//不使用CImg函数 在图上绘制一条长为100的直线段，起点坐标为(0, 0)，方向角为35度，直线的颜色为蓝色。 
		//下面三个函数分别对应使用CImg函数的上述三个操作
		void DrawCircle_blue2();
		void DrawCircle_yellow2();
		void DrawLine2();
		CImg<unsigned char> getSrcImg();
	private:
		//string name; //图片的名称
		CImg<unsigned char> SrcImg; //定义一副图片
};
#endif