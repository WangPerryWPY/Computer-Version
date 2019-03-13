#include "Canny.h"
#include <iostream>
#include "CImg-2.4.0_pre090618\CImg.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <algorithm>
#define pi 3.1415926
#define BOOSTBLURFACTOR 90.0
#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
using namespace std;
using namespace cimg_library;
void make_gaussian_kernel(float sigma, float **kernel, int *windowsize);
double angle_radians(double x, double y);
void follow_edges(int *edgemapptr, int *edgemagptr, int lowval, int cols);
Canny::Canny() {}
Canny::~Canny() {
	delete[] delta_x;
	delete[] delta_y;
	delete[] dirim;
	delete[] magnitude;
	delete[] nms;
	delete[] edge;
	delete[] smoothedim;
}



void Canny::HoughImageDetect(string name, float sigma, float tlow, float thigh, int Nums, int min_r, int max_r) {
	const char *a = name.c_str();
	img.load_bmp(a);
	CImg<int> thres1 = canny_image(sigma, tlow, thigh);
	//thres1.display();
	thres1.save("wow.bmp");
	CImg<int> thres_img1 = CircleDetect(Nums, min_r, max_r);
	string aa = name.substr(9,1) + "11.bmp";
	thres_img1.save(aa.c_str());
	//CImg<int> CircleImg = findpixelCircle(thres_img1, thres1);
	

	//thres_img1.display();
	//CircleImg.display();
	
	//CircleImg.save("b.bmp");
}

void Canny::RGBtoGray() {
	int r = 0, g = 0, b = 0;
	cimg_forXY(img, x, y) {
		r = img(x, y, 0);
		g = img(x, y, 1);
		b = img(x, y, 2);
		img(x, y, 0) = img(x, y, 1) = img(x, y, 2) = (r*0.2126 + g * 0.7152 + b * 0.0722);
	}
	img.resize(rows, cols, 1, 1, 5);
}

void Canny::gaussian_smooth(float sigma)
{
	int *tempim = new int[rows*cols];
	int r, c,rr,cc,
		windowsize,  //高斯核维度
		center;      //核中心
	float *kernel,
		dot,
		sum;
	make_gaussian_kernel(sigma, &kernel, &windowsize);
	center = windowsize / 2;
	for (r = 0; r < rows; r++) {
		for (c = 0; c < cols; c++) {
			dot = 0.0;
			sum = 0.0;
			for (cc = (-center); cc <= center; cc++) {
				if (((c + cc) >= 0) && ((c + cc) < cols)) {
					dot += (float)img(r,c+cc) * kernel[center + cc];
					sum += kernel[center + cc];
				}
			}
			tempim[r*cols + c] = dot / sum;
		}
	}
	for (c = 0; c < cols; c++) {
		for (r = 0; r < rows; r++) {
			sum = 0.0;
			dot = 0.0;
			for (rr = (-center); rr <= center; rr++) {
				if (((r + rr) >= 0) && ((r + rr) < rows)) {
					dot += tempim[(r + rr)*cols + c] * kernel[center + rr];
					sum += kernel[center + rr];
				}
			}
			smoothedim[r*cols + c] = (short int)(dot*BOOSTBLURFACTOR / sum + 0.5);
		}
	}
}

void Canny::derrivative_x_y() {
	int r = 0, c = 0, pos = 0;
	//计算x方向的一阶导数，判断边界避免遗失边界像素点
	for (r = 0; r < rows; r++) {
		pos = r * cols;
		delta_x[pos] = smoothedim[pos + 1] - smoothedim[pos];
		pos++;
		for (c = 1; c < (cols - 1); c++, pos++) {
			delta_x[pos] = smoothedim[pos + 1] - smoothedim[pos - 1];
		}
		delta_x[pos] = smoothedim[pos] - smoothedim[pos - 1];
	}
	//计算y方向的一阶导数，判断边界避免遗失边界像素点
	for (c = 0; c < cols; c++) {
		pos = c;
		delta_y[pos] = smoothedim[pos + cols] - smoothedim[pos];
		pos += cols;
		for (r = 1; r < (rows - 1); r++, pos += cols) {
			delta_y[pos] = smoothedim[pos + cols] - smoothedim[pos - cols];
		}
		delta_y[pos] = smoothedim[pos] - smoothedim[pos - cols];
	}
}

void Canny::radian_direction(int xdirtag, int ydirtag) {
	double dx = 0.0, dy = 0.0;
	int r = 0, c = 0, pos = 0;
	for (r = 0, pos = 0; r < rows; r++) {
		for (c = 0; c < cols; c++, pos++) {
			dx = (double)delta_x[pos];
			dy = (double)delta_y[pos];

			if (xdirtag == 1) dx = -dx;
			if (ydirtag == -1) dy = -dy;

			dirim[pos] = (float)angle_radians(dx, dy);
		}
	}
}

void Canny::magnitude_x_y() {
	int r = 0, c = 0, pos = 0, sq1 = 0, sq2 = 0;
	for (r = 0, pos = 0; r < rows; r++) {
		for (c = 0; c < cols; c++, pos++) {
			sq1 = (int)delta_x[pos] * (int)delta_x[pos];
			sq2 = (int)delta_y[pos] * (int)delta_y[pos];
			magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
		}
	}
}

void Canny::non_max_supp() {
	int rowcount = 0, colcount = 0, count = 0;
	int *magrowptr, *magptr;
	int *gxrowptr, *gxptr;
	int *gyrowptr, *gyptr, z1 = 0, z2 = 0;
	int m00, gx = 0, gy = 0;
	float mag1 = 0.0, mag2 = 0.0, xperp = 0.0, yperp = 0.0;
	int *resultrowptr, *resultptr;

	for (count = 0, resultrowptr = nms, resultptr = nms + cols * (rows - 1);
		count < cols; resultptr++, resultrowptr++, count++) {
		*resultrowptr = *resultptr = 0;
	}

	for (count = 0, resultptr = nms, resultrowptr = nms + cols - 1;
		count < rows; count++, resultptr += cols, resultrowptr += cols) {
		*resultptr = *resultrowptr = 0;
	}

	for (rowcount = 1, magrowptr = magnitude + cols + 1, gxrowptr = delta_x + cols + 1,
		gyrowptr = delta_y + cols + 1, resultrowptr = nms + cols + 1;
		rowcount < rows - 2;
		rowcount++, magrowptr += cols, gyrowptr += cols, gxrowptr += cols,
		resultrowptr += cols) {
		for (colcount = 1, magptr = magrowptr, gxptr = gxrowptr, gyptr = gyrowptr,
			resultptr = resultrowptr; colcount < cols - 2;
			colcount++, magptr++, gxptr++, gyptr++, resultptr++) {
			m00 = *magptr;
			if (m00 == 0) {
				*resultptr = NOEDGE;
			}
			else {
				xperp = -(gx = *gxptr) / ((float)m00);
				yperp = (gy = *gyptr) / ((float)m00);
			}

			if (gx >= 0) {
				if (gy >= 0) {
					if (gx >= gy)
					{
						/* 111 */
						/* Left point */
						z1 = *(magptr - 1);
						z2 = *(magptr - cols - 1);

						mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

						/* Right point */
						z1 = *(magptr + 1);
						z2 = *(magptr + cols + 1);

						mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
					}
					else
					{
						/* 110 */
						/* Left point */
						z1 = *(magptr - cols);
						z2 = *(magptr - cols - 1);

						mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

						/* Right point */
						z1 = *(magptr + cols);
						z2 = *(magptr + cols + 1);

						mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp;
					}
				}
				else
				{
					if (gx >= -gy)
					{
						/* 101 */
						/* Left point */
						z1 = *(magptr - 1);
						z2 = *(magptr + cols - 1);

						mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

						/* Right point */
						z1 = *(magptr + 1);
						z2 = *(magptr - cols + 1);

						mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
					}
					else
					{
						/* 100 */
						/* Left point */
						z1 = *(magptr + cols);
						z2 = *(magptr + cols - 1);

						mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

						/* Right point */
						z1 = *(magptr - cols);
						z2 = *(magptr - cols + 1);

						mag2 = (z1 - z2)*xperp + (m00 - z1)*yperp;
					}
				}
			}
			else
			{
				if ((gy = *gyptr) >= 0)
				{
					if (-gx >= gy)
					{
						/* 011 */
						/* Left point */
						z1 = *(magptr + 1);
						z2 = *(magptr - cols + 1);

						mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

						/* Right point */
						z1 = *(magptr - 1);
						z2 = *(magptr + cols - 1);

						mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
					}
					else
					{
						/* 010 */
						/* Left point */
						z1 = *(magptr - cols);
						z2 = *(magptr - cols + 1);

						mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

						/* Right point */
						z1 = *(magptr + cols);
						z2 = *(magptr + cols - 1);

						mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
					}
				}
				else
				{
					if (-gx > -gy)
					{
						/* 001 */
						/* Left point */
						z1 = *(magptr + 1);
						z2 = *(magptr + cols + 1);

						mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

						/* Right point */
						z1 = *(magptr - 1);
						z2 = *(magptr - cols - 1);

						mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
					}
					else
					{
						/* 000 */
						/* Left point */
						z1 = *(magptr + cols);
						z2 = *(magptr + cols + 1);

						mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

						/* Right point */
						z1 = *(magptr - cols);
						z2 = *(magptr - cols - 1);

						mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
					}
				}
			}

			/* Now determine if the current point is a maximum point */

			if ((mag1 > 0.0) || (mag2 > 0.0))
			{
				*resultptr = NOEDGE;
			}
			else
			{
				if (mag2 == 0.0)
					*resultptr = NOEDGE;
				else
					*resultptr = POSSIBLE_EDGE;
			}
		}
	}
}

void Canny::apply_hysteresis(float tlow, float thigh) {
	int r = 0, c = 0, pos = 0, numedges = 0, lowcount = 0, highcount = 0, lowthreshold = 0, highthreshold = 0,
		i = 0, *hist, rr = 0, cc = 0;
	hist = new int[32768];
	int maximum_mag = 0, sumpix = 0;
	for (r = 0, pos = 0; r < rows; r++) {
		for (c = 0; c < cols; c++, pos++) {
			if (nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
			else edge[pos] = NOEDGE;
		}
	}

	for (r = 0, pos = 0; r < rows; r++, pos += cols) {
		edge[pos] = NOEDGE;
		edge[pos + cols - 1] = NOEDGE;
	}
	pos = (rows - 1) * cols;
	for (c = 0; c < cols; c++, pos++) {
		edge[c] = NOEDGE;
		edge[pos] = NOEDGE;
	}
	for (r = 0; r < 32768; r++) hist[r] = 0;
	for (r = 0, pos = 0; r < rows; r++) {
		for (c = 0; c < cols; c++, pos++) {
			if (edge[pos] == POSSIBLE_EDGE) hist[magnitude[pos]]++;
		}
	}
	for (r = 1, numedges = 0; r < 32768; r++) {
		if (hist[r] != 0) maximum_mag = r;
		numedges += hist[r];
	}

	highcount = (int)(numedges * thigh + 0.5);

	r = 1;
	numedges = hist[1];
	while ((r < (maximum_mag - 1)) && (numedges < highcount)) {
		r++;
		numedges += hist[r];
	}
	highthreshold = r;
	lowthreshold = (int)(highthreshold * tlow + 0.5);

	for (r = 0, pos = 0; r < rows; r++) {
		for (c = 0; c < cols; c++, pos++) {
			if ((edge[pos] == POSSIBLE_EDGE) && (magnitude[pos] >= highthreshold)) {
				edge[pos] = EDGE;
				follow_edges((edge + pos), (magnitude + pos), lowthreshold, cols);
			}
		}
	}

	for (r = 0, pos = 0; r < rows; r++) {
		for (c = 0; c < cols; c++, pos++) if (edge[pos] != EDGE) edge[pos] = NOEDGE;
	}
	delete[] hist;
}

void make_gaussian_kernel(float sigma, float **kernel, int *windowsize) {
	int i = 0, center = 0;
	float x, fx, sum = 0.0;
	//根据高斯滤波核的方差计算高斯核的宽高
	*windowsize = 1 + 2 * ceil(2.5 * sigma);
	//*kernel = (float*)calloc((*windowsize), sizeof(float));
	center = (*windowsize) / 2;
	*kernel = new float[*windowsize];
	for (i = 0; i < (*windowsize); i++) {
		x = (float)(i - center);
		fx = pow(2.71828, -0.5*x*x / (sigma*sigma)) / (sigma * sqrt(6.2831853));
		(*kernel)[i] = fx;
		sum += fx;
	}
	for (i = 0; i < (*windowsize); i++) {
		(*kernel)[i] /= sum;
	}
}

double angle_radians(double x, double y) {
	double xu = 0.0, yu = 0.0, ang = 0.0;
	xu = fabs(x);
	yu = fabs(y);
	if ((xu == 0) && (yu == 0)) return(0);
	ang = atan(yu / xu);
	if (x >= 0) {
		if (y >= 0) return (ang);
		else return(2 * M_PI - ang);
	}
	else {
		if (y >= 0) return (M_PI - ang);
		else return(M_PI + ang);
	}
}

void follow_edges(int *edgemapptr, int *edgemagptr, int lowval, int cols)
{
	int *tempmagptr;
	int *tempmapptr;
	int i;
	float thethresh;
	int x[8] = { 1,1,0,-1,-1,-1,0,1 },
		y[8] = { 0,1,1,1,0,-1,-1,-1 };

	for (i = 0; i < 8; i++) {
		tempmapptr = edgemapptr - y[i] * cols + x[i];
		tempmagptr = edgemagptr - y[i] * cols + x[i];

		if ((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)) {
			*tempmapptr = EDGE;
			follow_edges(tempmapptr, tempmagptr, lowval, cols);
		}
	}
}

CImg<int> Canny::canny_image(int sigma, float tlow, float thigh) {
	rows = img.width();
	cols = img.height();
	delta_x = new int[rows*cols]; memset(delta_x, 0x0, rows*cols*sizeof(int));
	delta_y = new int[rows*cols]; memset(delta_y, 0x0, rows*cols * sizeof(int));
	dirim = new float[rows*cols]; memset(dirim, 0x0, rows*cols * sizeof(float));
	magnitude = new int[rows*cols]; memset(magnitude, 0x0, rows*cols * sizeof(int));
	nms = new int[rows*cols]; memset(nms, 0x0, rows*cols * sizeof(int));
	edge = new int[rows*cols]; memset(edge, 0x0, rows*cols * sizeof(int));
	smoothedim = new int[rows*cols];  memset(smoothedim, 0x0, rows*cols * sizeof(int));

	RGBtoGray();
	gaussian_smooth(sigma);
	derrivative_x_y();
	radian_direction(-1, -1);
	magnitude_x_y();
	non_max_supp();
	apply_hysteresis(tlow, thigh);
	//img.test();
	//CImg<int> pic(rows, cols, 1, 1, 5);
	//pic.fill(0);
	thres.resize(rows, cols, 1, 1, 5);
	thres.fill(0);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (edge[i*cols + j] < 0)
				thres(i, j) = 0;
			else if (edge[i*cols + j] > 255)
				thres(i, j) = 255;
			else
				thres(i, j) = edge[i*cols + j];
		}
	}
	return thres;
}

CImg<int> Canny::CircleDetect(int Nums, int min_r, int max_r) {
	thres_img = thres;
	thres_img.resize(thres.width(), thres.height(), 1, 3);
	vector<int> sortCircleWeight;
	int max = 0;
	int width = thres_img.width();
	int height = thres_img.height();
	CImg<int> pic = thres;
	vector<pair<int, int>> vote;
	for (int r = min_r; r < max_r; r += 5) {
		HoughImg_Circle.resize(width, height, 1, 1, 0);
		HoughImg_Circle.fill(0);
		max = 0;
		cimg_forXY(pic, x, y) {
			if (pic(x, y) == 0) {
				for (int i = 0; i < 360; i++) {
					int x0 = x - r * cos(i*M_PI / 180);
					int y0 = y - r * sin(i*M_PI / 180);
					if (x0 > 0 && x0 < width && y0 > 0 && y0 < height)
						HoughImg_Circle(x0, y0)++;
				}
			}
		}
		cimg_forXY(HoughImg_Circle, x, y) {
			if (HoughImg_Circle(x, y) > max) {
				max = HoughImg_Circle(x, y);
			}
		}
		vote.push_back(make_pair(max, r));
	}
	sort(vote.begin(), vote.end(), [](const pair<int, int>& x, const pair<int, int>& y) -> int { return x.first > y.first; });
	//int Nums = 1;
	int knums = 0;
	for (int num = 0; num < Nums; num++) {
		int i = 0;
		HoughImg_Circle.resize(width, height, 1, 1, 0);
		HoughImg_Circle.fill(0);
		int r = vote[num].second;
		cimg_forXY(pic, x, y) {
			if (pic(x, y) == 0) {
				for (int i = 0; i < 360; i++) {
					int x0 = x - r * cos(i*M_PI / 180);
					int y0 = y - r * sin(i*M_PI / 180);
					if (x0 > 0 && x0 < width && y0 > 0 && y0 < height)
						HoughImg_Circle(x0, y0)++;
				}
			}
		}
		//Circle.clear();
		//Circleweight.clear();
		cimg_forXY(HoughImg_Circle, x, y) {
			if (HoughImg_Circle(x, y) != 0) {
				Circle.push_back(make_pair(x, y));
				Circleweight.push_back(HoughImg_Circle(x, y));
			}
		}
		unsigned char blue[3] = { 0, 0, 255 };
		sortCircleWeight = Circleweight;
		sort(sortCircleWeight.begin(), sortCircleWeight.end(), greater<int>()); // 将累加矩阵从大到小进行排序
		//同个半径圆的检测
		int count = 0;
		while (1) {
			int weight = sortCircleWeight[count], index;
			vector<int>::iterator iter = find(Circleweight.begin(), Circleweight.end(), weight);
			index = iter - Circleweight.begin();
			int a = Circle[index].first, b = Circle[index].second;
			count++;

			int ii;
			for (ii = 0; ii < center.size(); ii++) {
				// 判断检测出来的圆心坐标是否跟已检测的圆心坐标的距离，如果距离过小，默认是同个圆
				if (sqrt(pow((center[ii].first - a), 2) + pow((center[ii].second - b), 2)) < 100) {
					break;
				}
			}
			if (ii == center.size()) {
				center.push_back(make_pair(a, b));
				thres_img.draw_circle(a, b, r, blue, 1, 1);
				knums++;
				break;
			}
		}
	}
	cout << "the number of the coins is: " << knums << endl;
	//thres_img.display();
	return thres_img;
}

CImg<int> Canny::findpixelCircle(CImg<int> thres_img1, CImg<int> thres1) {
	CImg<int> CircleImg = img;
	//CircleImg.resize(thres.width(), thres.height());
	cimg_forXY(CircleImg, x, y) {
		if (thres_img(x, y, 0) == 0 && thres_img(x, y, 1) == 0 && thres_img(x, y, 2) == 255 && thres(x, y) == 0) {
			CircleImg(x, y, 0) = 255;
			CircleImg(x, y, 1) = 0;
			CircleImg(x, y, 2) = 0;
		}
	}
	return CircleImg;
}

