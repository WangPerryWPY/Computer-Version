#include <algorithm>
#include "pch.h"
#include <iostream>
#include "CImg.h"
#include <string>
#include <cstring>
#include <math.h>
#include <map>
#include "ImgConverse.h"
#include "canny.h"
using namespace cimg_library;
using namespace std;

CImg<int> test(CImg<int> a) {
	cimg_forXY(a, x, y) {
		a(x, y, 0) = 255 - a(x, y, 0);
		a(x, y, 1) = 255 - a(x, y, 1);
		a(x, y, 2) = 255 - a(x, y, 2);
	}
	return a;
}

int main()
{
	ImgConverse ImgCon;
	//ImgCon.Hist_gray();
	//ImgCon.Hist_impress();
	//ImgCon.Hist_color();
	//ImgCon.RGB_HSI();
	ImgCon.colorTransfer();
	return 0;
}