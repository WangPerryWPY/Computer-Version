#include <iostream>
#include "CImg.h"
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <dirent.h>
#define PI 3.14159265358
#define STITCH_NUM 20
using namespace cimg_library;
using namespace std;

extern "C" {
    #include "vl/generic.h"
    #include "vl/sift.h"
    #include "vl/kdtree.h"
}

//sift点对
struct point_pair {
    VlSiftKeypoint a;
    VlSiftKeypoint b;
    point_pair(VlSiftKeypoint _a, VlSiftKeypoint _b) {
        a = _a;
            b = _b;
    }
};

//单应矩阵点对
struct points {
    float x1,x2,x3,x4,x5,x6,x7,x8;
    points(float _x1, float _x2, float _x3, float _x4, float _x5, float _x6, float _x7, float _x8) {
        x1 = _x1;
        x2 = _x2;
        x3 = _x3;
        x4 = _x4;
        x5 = _x5;
        x6 = _x6;
        x7 = _x7;
        x8 = _x8;
    }
};


class ImageStitching {

private:
    //源图列表
    CImgList<float> imgs;
    //拼接后的图像
    CImg<float> resultImg;

private:
    //柱面投影
    CImg<float> CylindricalProjection(CImg<float> pic);
    //图像从rgb空间转为灰度空间
    CImg<float> convertTogray(CImg<float> pic);
    //提取图像中的特征点
    map<vector<float>, VlSiftKeypoint> SIFTFeatures(CImg<float> pic);
    //k-d树进行特征点匹配
    vector<point_pair> KDtreeMatch(map<vector<float>, VlSiftKeypoint> feature_a, map<vector<float>, VlSiftKeypoint> feature_b);
    //利用RANSAC算法求单应矩阵
    points RANSAC(vector<point_pair> pairs);
    //通过单应矩阵扭曲两幅图像内容
    void WarpTwoImg(CImg<float> src, CImg<float> &dst, points H, float offset_x, float offset_y);
    //移动两幅图像
    void MoveTwoImg(CImg<float> src, CImg<float> &dst, int offset_x, int offset_y);
    //图像拼接
    CImg<float> Blend(CImg<float> pic1, CImg<float> pic2);

public:
    ImageStitching();
    ~ImageStitching();
    void StitchProcess();
};