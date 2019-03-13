#ifndef _FUNCTION_H_
#define _FUNCTION_H_
#include <iostream>
#include "CImg.h"
using namespace cimg_library;
using namespace std;
const int filter[3][3] = {
							{0, 1, 0},
							{1, 0 ,1},
							{0, 1, 0}
							};
CImg<float> rgb2gray(CImg<float> );
CImg<float> dilate_white(CImg<float>, int);
CImg<float> dilate_black(CImg<float>, int);
CImg<float> dilate_a(CImg<float> );
CImg<float> dilate_b(CImg<float> );
CImg<float> dilate_X(CImg<float> );
CImg<float> dilate_Y(CImg<float> );
#endif