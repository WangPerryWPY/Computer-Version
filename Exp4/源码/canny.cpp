#include "pch.h"
#include "canny.h"
#include <stdio.h>
#include <math.h>
#define M_PI 3.141592653589793

canny::canny()
{
}


canny::~canny()
{
}

CImg<int> canny::process(const char *filepath, int scale, double sigma, int thresh_low, int thresh_high) {
	img.load_bmp(filepath);

	toGrayScale();

	useFilter(scale, sigma);

	sobel();

	nonMaxSupp();

	threshold(thresh_low, thresh_high);
	return thres;
}

vector<vector<double> > canny::createFilter(int scale, double sigma) {
	vector<vector<double> > filter;

	// initialize
	int s = 2 * scale + 1;
	for (int i = 0; i < s; ++i) {
		vector<double> temp;
		for (int j = 0; j < s; j++) {
			temp.push_back(-1);
		}
		filter.push_back(temp);
	}

	float sum = 0;
	for (int x = -s / 2; x <= s / 2; ++x) {
		for (int y = -s / 2; y <= s / 2; ++y) {
			filter[x + s / 2][y + s / 2] = (exp(-(x * x + y * y) / (2.0 * sigma * sigma))) * (1 / (M_PI * (2.0 * sigma * sigma)));
			sum += filter[x + s / 2][y + s / 2];
		}
	}

	// normalize the Filter
	for (int i = 0; i < s; ++i) {
		for (int j = 0; j < s; ++j) {
			filter[i][j] /= sum;
		}
	}

	return filter;
}

vector<vector<double> > canny::createSobelFilterX() {
	// Sobel x filter
	vector<vector<double> > x_filter(3);
	x_filter[0] = { -1.0, 0, 1.0 };
	x_filter[1] = { -2.0, 0, 2.0 };
	x_filter[2] = { -1.0, 0, 1.0 };
	return x_filter;
}

vector<vector<double> > canny::createSobelFilterY() {
	// Sobel y filter
	vector<vector<double> > y_filter(3);
	y_filter[0] = { 1.0, 2.0, 1.0 };
	y_filter[1] = { 0, 0, 0 };
	y_filter[2] = { -1.0, -2.0, -1.0 };
	return y_filter;
}

void canny::toGrayScale()
{
	grayscaled = CImg<unsigned char>(img.width(), img.height(), 1);
	cimg_forXY(img, x, y) {
		unsigned char r = img(x, y, 0);
		unsigned char g = img(x, y, 1);
		unsigned char b = img(x, y, 2);
		double gray = (r * 0.2126 + g * 0.7152 + b * 0.0722);
		grayscaled.atXY(x, y) = gray;
	}
}

void canny::useFilter(int scale, double sigma) {
	vector<vector<double> > filter = createFilter(scale, sigma);
	int size = (int)filter.size() / 2;
	CImg<unsigned char> grayImg(grayscaled.width() + 2 * size, grayscaled.height() + 2 * size, 1, 1,0);
	cimg_forXY(grayImg, x, y) {
		grayImg(x, y) = 0;
	}
	for (int i = size; i < grayImg.height() - size; ++i) {
		for (int j = size; j < grayImg.width() - size; ++j) {
			grayImg(j, i) = grayscaled(j - size, i - size);
		}
	}
	grayscaled = grayImg;
	gFiltered = CImg<unsigned char>(grayscaled.width() - 2 * size, grayscaled.height() - 2 * size, 1);
	for (int i = size; i < grayscaled.height() - size; ++i) {
		for (int j = size; j < grayscaled.width() - size; ++j) {
			double sum = 0;
			for (int x = 0; x < filter.size(); ++x) {
				for (int y = 0; y < filter.size(); ++y) {
					sum += filter[x][y] * (double)(grayscaled.atXY(j + y - size, i + x - size));
				}
			}
			gFiltered.atXY(j - size, i - size) = sum;
		}
	}
}

void canny::sobel() {
	// Sobel filter
	vector<vector<double> > x_filter = createSobelFilterX();
	vector<vector<double> > y_filter = createSobelFilterY();

	int size = (int)x_filter.size() / 2;  // limit size
	const int width = gFiltered.width() - 2 * size;
	const int height = gFiltered.height() - 2 * size;
	sFiltered = CImg<unsigned char>(width, height, 1);
	angles = CImg<float>(width, height, 1);

	for (int i = size; i < gFiltered.height() - size; ++i) {
		for (int j = size; j < gFiltered.width() - size; ++j) {
			double sumx = 0;
			double sumy = 0;
			for (int x = 0; x < x_filter.size(); ++x) {
				for (int y = 0; y < y_filter.size(); ++y) {
					sumx += x_filter[x][y] * (double)(gFiltered.atXY(j + y - size, i + x - size));
					sumy += y_filter[x][y] * (double)(gFiltered.atXY(j + y - size, i + x - size));
				}
			}
			double sumxsq = sumx * sumx;
			double sumysq = sumy * sumy;
			double sq2 = sqrt(sumxsq + sumysq);
			if (sq2 > 255) sq2 = 255;
			sFiltered.atXY(j - size, i - size) = sq2;
			if (sumx == 0) {
				angles.atXY(j - size, i - size) = 90;
			}
			else {
				angles.atXY(j - size, i - size) = atan(sumy / sumx);
			}
		}
	}
}

void canny::nonMaxSupp() {
	non = CImg<unsigned char>(sFiltered.width() - 2, sFiltered.height() - 2, 1);
	for (int i = 1; i < sFiltered.height() - 1; ++i) {
		for (int j = 1; j < sFiltered.width() - 1; ++j) {
			float tangent = angles.atXY(j, i);
			non.atXY(j - 1, i - 1) = sFiltered.atXY(j, i);
			// horizontal edge
			if (((-22.5 < tangent) && (tangent <= 22.5)) || ((157.5 < tangent) && (tangent <= -157.5))) {
				if ((sFiltered.atXY(j, i) < sFiltered.atXY(j + 1, i)) || (sFiltered.atXY(j, i) < sFiltered.atXY(j - 1, i)))
					non.atXY(j - 1, i - 1) = 0;
			}
			// vertical edge
			if (((-112.5 < tangent) && (tangent <= -67.5)) || ((67.5 < tangent) && (tangent <= 112.5))) {
				if ((sFiltered.atXY(j, i) < sFiltered.atXY(j, i + 1)) || (sFiltered.atXY(j, i) < sFiltered.atXY(j, i - 1)))
					non.atXY(j - 1, i - 1) = 0;
			}
			// -45 degree edge
			if (((-67.5 < tangent) && (tangent <= -22.5)) || ((112.5 < tangent) && (tangent <= 157.5))) {
				if ((sFiltered.atXY(j, i) < sFiltered.atXY(j + 1, i - 1)) || (sFiltered.atXY(j, i) < sFiltered.atXY(j - 1, i + 1)))
					non.atXY(j - 1, i - 1) = 0;
			}
			// 45 degree edge
			if (((-157.5 < tangent) && (tangent <= -112.5)) || ((22.5 < tangent) && (tangent <= 67.5))) {
				if ((sFiltered.atXY(j, i) < sFiltered.atXY(j + 1, i + 1)) || (sFiltered.atXY(j, i) < sFiltered.atXY(j - 1, i - 1)))
					non.atXY(j - 1, i - 1) = 0;
			}
		}
	}
}

void canny::threshold(int low, int high) {
	thres = CImg<unsigned char>(non);
	for (int i = 0; i < non.height(); ++i) {
		for (int j = 0; j < non.width(); ++j) {
			if (thres.atXY(j, i) > high) {
				thres.atXY(j, i) = 255;
			}
			else if (thres.atXY(j, i) < low) {
				thres.atXY(j, i) = 0;
			}
			else {
				bool anyHigh = false;
				bool anyBetween = false;
				for (int x = i - 1; x < i + 2; ++x) {
					for (int y = j - 1; y < j + 2; ++y) {
						if (x <= 0 || y <= 0 || x > thres.height() || y > thres.width()) {
							continue;
						}
						else {
							if (thres.atXY(y, x) > high) {
								thres.atXY(j, i) = 255;
								anyHigh = true;
								break;
							}
							else if (thres.atXY(y, x) <= high && thres.atXY(y, x) >= low) {
								anyBetween = true;
							}
						}
					}
					if (anyHigh) break;
				}
				if (!anyHigh && anyBetween) {
					for (int x = i - 2; x < i + 3; ++x) {
						for (int y = j - 1; y < j + 3; ++y) {
							if (x < 0 || y < 0 || x > thres.height() || y > thres.width()) {
								continue;
							}
							else {
								if (thres.atXY(y, x) > high) {
									thres.atXY(j, i) = 255;
									anyHigh = true;
									break;
								}
							}
						}
						if (anyHigh) {
							break;
						}
					}
				}
				if (!anyHigh) {
					thres.atXY(j, i) = 0;
				}
			}
		}
	}
}