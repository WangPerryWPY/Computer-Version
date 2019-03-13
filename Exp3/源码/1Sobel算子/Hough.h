#include <string>
#include <vector>
#include <cmath>
#include "CImg-2.4.0_pre090618\CImg.h"

using namespace std;
using namespace cimg_library;

typedef struct parameter {
  int x; 
  int y; 
}parameter;

typedef struct lineParameter {
  double k;
  double b;
  bool operator < (const lineParameter a) const {
    return abs(k) < abs(a.k);
  }
  bool operator > (const lineParameter a) const {
    return abs(k) > abs(a.k);
  }
  /*bool operator == (const lineParameter a) const {
    return k == a.k && b == a.b;
  }*/
}lineParameter;

class Hough {
private:
  CImg<unsigned char> img; //原图像
  CImg<unsigned char> grayscaled; //灰度图
  CImg<unsigned char> gFiltered; // 高斯滤波
  CImg<unsigned char> sFiltered; // Sobel滤波
  CImg<float> angles; // 角度
  CImg<unsigned char> non; // 非极大值抑制
  CImg<unsigned char> thres; // 双阈值
  CImg<int> HoughImg; //直线的霍夫空间图像
  vector<lineParameter> HoughLine; //直线检测得出的直线的参数
  vector<lineParameter> Hough_Line; //临时存储直线参数
  vector<int> nums;  //用于筛选出最长的四条边
  vector<parameter> HoughPoint; //直线检测的角点
  CImg<unsigned char> HoughImg_Circle; //霍夫空间圆图像
  CImg<unsigned char> CircleImg; //在霍夫圆图像上生成的像素点
  CImg<unsigned char> thres_img; //边缘检测图像的圆检测图
public:
  Hough(); 
  //直线检测操作
  void LineDetect_Image(string, int, double, int, int, int, int, int);
  //圆检测操作
  void CircleDetect_Image(string, int, double, int, int, int, int, int);
  vector<vector<double> > createFilter(int, double); //高斯滤波
  vector<vector<double> > createSobelFilterX();
  vector<vector<double> > createSobelFilterY();
  void toGrayScale();
  void useFilter(int, double);
  void sobel(); 
  void nonMaxSupp(); 
  void threshold(int, int);
  //将直线变换为极坐标
  void Polar_Line();
  //霍夫直线检测
  CImg<unsigned char> LineDetect(int, int, int);
  //描绘角点
  CImg<unsigned char> drawPoint();
  //霍夫圆检测
  void CircleDetect(int, int, int);
  //查找圆对应的像素点
  void findpixelCircle();
};
