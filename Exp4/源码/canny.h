#pragma once
#include <string>
#include <vector>
#include "CImg.h"
#include <math.h>
using namespace std;
using namespace cimg_library;

class canny
{
private:

	CImg<unsigned char> img; //Original Image
	CImg<unsigned char> grayscaled; // Grayscale
	CImg<unsigned char> gFiltered; // Gauss Filtered
	CImg<unsigned char> sFiltered; // Sobel Filtered
	CImg<float> angles; // Angle Map
	CImg<unsigned char> non; // Non-maxima supp.
	CImg<unsigned char> thres; // Double threshold and final
public:
	canny();
	~canny();
	CImg<int> process(const char *, int, double, int, int); // Process image
  vector<vector<double> > createFilter(int, double); // Creates a gaussian filter
  vector<vector<double> > createSobelFilterX();
  vector<vector<double> > createSobelFilterY();
  void toGrayScale(); // Gray image
  void useFilter(int, double); // Gauss filtering
  void sobel(); // Sobel filtering
  void nonMaxSupp(); // Non-maxima supp.
  void threshold(int, int); // Double threshold and finalize picture
};

