#include "Hough.h"
#include <stdio.h>
#include <string>
#include <vector>
#include "CImg-2.4.0_pre090618\CImg.h"

using namespace std;
using namespace cimg_library;
int main() {
  int scale = 1, thresh_low = 60, thresh_high = 100;
  double sigma = 6.0;
  Hough hough1;
  hough1.LineDetect_Image("Dataset1/1.bmp",scale, sigma, thresh_low, thresh_high, 4, 0, 50);
  Hough hough2;
  hough2.LineDetect_Image("Dataset1/2.bmp",scale, sigma, thresh_low, thresh_high, 4, 2, 50);
  Hough hough3;
  hough3.LineDetect_Image("Dataset1/3.bmp",scale, sigma, thresh_low, thresh_high, 4, 0, 50);
  Hough hough4;
  hough4.LineDetect_Image("Dataset1/4.bmp",scale, sigma, thresh_low, thresh_high, 4, 2, 50);
  Hough hough5;
  hough5.LineDetect_Image("Dataset1/5.bmp",scale, sigma, thresh_low, thresh_high, 4, 0, 50);
  Hough hough6;
  hough6.LineDetect_Image("Dataset1/6.bmp",scale, sigma, thresh_low, thresh_high, 4, 0, 50);
  
  Hough hough12;
  hough12.CircleDetect_Image("Dataset2/1.bmp",scale, sigma, thresh_low, thresh_high, 1, 150, 170);
  Hough hough22;
  hough22.CircleDetect_Image("Dataset2/2.bmp",scale, sigma, thresh_low, thresh_high, 4, 190, 230);
  Hough hough32;
  hough32.CircleDetect_Image("Dataset2/3.bmp",scale, sigma, thresh_low, thresh_high, 7, 140, 200);
  Hough hough42;
  hough42.CircleDetect_Image("Dataset2/4.bmp",scale, sigma, thresh_low, thresh_high, 3, 140, 210);
  Hough hough52;
  hough52.CircleDetect_Image("Dataset2/5.bmp",scale, sigma, thresh_low, thresh_high, 2, 460, 530);
  Hough hough62;
  hough62.CircleDetect_Image("Dataset2/6.bmp",scale, sigma, thresh_low, thresh_high, 5, 40, 60);
  return 0;
}