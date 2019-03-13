#include "ex.hpp"
using namespace std;
bool cmp(double x , double y);

Test::Test() {
	//name = "work1";
	SrcImg.load_bmp("1.bmp");
}

Test::~Test() {}

void Test::Todisplay() {
	//string t = name;
	SrcImg.display("homework");
}

CImg<unsigned char> Test::getSrcImg() {
	return SrcImg;
}

void Test::change() {
	//name = "work2";
	//CImg<unsigned char> Img = SrcImg;
	cimg_forXY(SrcImg, x, y) {
		if (SrcImg(x,y,0) == 255 && SrcImg(x,y,1) == 255 && SrcImg(x,y,2) == 255) {
			SrcImg(x,y,0) = 255;
			SrcImg(x,y,1) = 0;
			SrcImg(x,y,2) = 0;
 		}
	}
	cimg_forXY(SrcImg, x, y) {
		if (SrcImg(x,y,0) == 0 && SrcImg(x,y,1) == 0 && SrcImg(x,y,2) == 0) {
			SrcImg(x,y,0) = 0;
			SrcImg(x,y,1) = 255;
			SrcImg(x,y,2) = 0;
 		}
	}
}

void Test::DrawCircle_blue1() {
	cimg_forXY(SrcImg, x, y) {
		if (pow(pow(x-50,2)+pow(y-50,2),0.5) < 30) {
			SrcImg(x,y,0) = 0;
			SrcImg(x,y,1) = 0;
			SrcImg(x,y,2) = 255;
		}
	}
}

void Test::DrawCircle_yellow1() {
	cimg_forXY(SrcImg, x, y) {
		if (pow(pow(x-50,2)+pow(y-50,2),0.5) < 3) {
			SrcImg(x,y,0) = 200;
			SrcImg(x,y,1) = 155;
			SrcImg(x,y,2) = 0;
		}
	}
}

void Test::DrawLine1() {
	double x0 = 100*cos(35*pi/180);
	double y0 = 100*sin(35*pi/180);
	cimg_forXY(SrcImg, x, y) {
		if (x == 0) {
			if (y == 0) {
				SrcImg(x,y,0) = 0;
				SrcImg(x,y,1) = 0;
				SrcImg(x,y,2) = 255;
			}
		}
		else {
			if (cmp((double)y, (double)x*tan(35*pi/180)) && (double)x <= x0 && (double)y <= y0) {
				SrcImg(x,y,0) = 0;
				SrcImg(x,y,1) = 0;
				SrcImg(x,y,2) = 255;
			}
		}
	}
}

void Test::DrawCircle_blue2() {
	unsigned char blue[] = {0,0,255};
	SrcImg.draw_circle(50, 50, 30, blue);
}

void Test::DrawCircle_yellow2() {
	unsigned char yellow[] = {200, 155, 0};
	SrcImg.draw_circle(50, 50, 3, yellow);
}
		
void Test::DrawLine2() {
	unsigned char blue[] = {0,0,255};
	SrcImg.draw_line(0,0,100*cos(35*pi/180),100*sin(35*pi/180),blue);
}

bool cmp(double x , double y) { //compare x and y，如果差值小于一定范围则近似相等
	if (abs(x - y) <= 0.5)
		return 1;
	return 0;
}