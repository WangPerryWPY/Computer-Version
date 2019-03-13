#include "canny.h"

canny::canny() {}

canny::canny(const char *filename) {
	Img.load(filename);
}

canny::canny(CImg<float> picture) {
	Img = picture;
}

canny::~canny() {}

void canny::Gauss_blur(float BLUR) {
	blurImg = grayImg.get_blur(BLUR);
	//blurImg.display();
}

void canny::sobel() {
	sobelImg.resize(blurImg._width, blurImg._height, 1, 1, 0);
	angleImg = sobelImg;
	CImg_3x3(I, float);
	cimg_for3x3(blurImg, x, y, 0, 0, I, float) {
		float ix = (Inp + 2 * Inc + Inn) - (Ipp + 2 * Ipc + Ipn);
		float iy = (Ipp + 2 * Icp + Inp) - (Ipn + 2 * Icn + Inn);
		float grad = sqrt(ix * ix + iy * iy);
		angleImg(x, y) = (ix == 0) ? 90 : atan(iy / ix);
		if (grad > 255) {
			grad = 255;
		}
		if (grad < 0) {
			grad = 0;
		}
		sobelImg(x, y) = grad;
	}
	//sobelImg.display();
	//angleImg.display();
}

void canny::nonMaxSupp() {
	nonImg = sobelImg;
	//sobelImg.display();
	cimg_for_insideXY(nonImg, x, y, 1) {
		// horizontal edge
		if (((-22.5 < angleImg(x, y)) && (angleImg(x, y) <= 22.5)) || ((157.5 < angleImg(x, y)) && (angleImg(x, y) <= -157.5))) {
	     	if (sobelImg(x, y) < sobelImg(x + 1, y) || sobelImg(x, y) < sobelImg(x - 1, y))
	        	nonImg(x, y) = 0;
	    }
	    // vertical edge
	    if (((-112.5 < angleImg(x, y)) && (angleImg(x, y) <= -67.5)) || ((67.5 < angleImg(x, y)) && (angleImg(x, y) <= 112.5))) {
	     	if (sobelImg(x, y) < sobelImg(x, y + 1) || sobelImg(x, y) < sobelImg(x, y - 1))
	       		nonImg(x, y) = 0;
	    }
	    // -45 degree edge
	    if (((-67.5 < angleImg(x, y)) && (angleImg(x, y) <= -22.5)) || ((112.5 < angleImg(x, y)) && (angleImg(x, y) <= 157.5))) {
	     	if (sobelImg(x, y) < sobelImg(x + 1, y - 1) || sobelImg(x, y) < sobelImg(x - 1, y + 1))
	          	nonImg(x, y) = 0;
	    }
	    // 45 degree edge
	    if (((-157.5 < angleImg(x, y)) && (angleImg(x, y) <= -112.5)) || ((22.5 < angleImg(x, y)) && (angleImg(x, y) <= 67.5))) {
	     	if (sobelImg(x, y) < sobelImg(x + 1, y + 1) || sobelImg(x, y) < sobelImg(x - 1, y - 1))
	          	nonImg(x, y) = 0;
	    } 
	}
	//nonImg.display();
}

void canny::threshold(int low_thres, int high_thres) {
	thresImg = nonImg;
	cimg_forXY(thresImg, x, y) {
		if (thresImg(x, y) > high_thres) {
			thresImg(x, y) = 255;
		}
		else if (thresImg(x, y) < low_thres) {
			thresImg(x, y) = 0;
		}
		else {
			bool anyHigh = false;
			bool anyBetween = false;
			for (int i = x - 1; i <= x + 1; ++i) {
				for (int j = y - 1; j <= y + 1; ++j) {
					if (i < 0 || i >= thresImg._width || j < 0 || j >= thresImg._height) {
						continue;
					}
					else {
						if (thresImg(i, j) > high_thres) {
							thresImg(x, y) = 255;
							anyHigh = true;
							break;
						}
						else if (thresImg(i, j) <= high_thres && thresImg(i, j) >= low_thres) {
							anyBetween = true;
						}
					}
				}
				if (anyHigh) break;
			}
			if (!anyHigh && anyBetween) {
				for (int i = x - 2; i <= x + 2; ++i) {
					for (int j = y - 2; j <= y + 2; ++j) {
						if (i < 0 || i >= thresImg._width || j < 0 || j >= thresImg._height) {
							continue;
						}
						else {
							if (thresImg(i, j) > high_thres) {
								thresImg(x, y) = 255;
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
				thresImg(x, y) = 0;
			}
		}
	}
	//thresImg.display("cannyImg");
}

CImg<float> canny::process(float BLUR, int low_thres, int high_thres) {
	if (Img.spectrum() == 3) {
		//转灰度图
		grayImg = rgb2gray(Img);
	}
	else {
		grayImg = Img;
	}
	//高斯滤波
	Gauss_blur(BLUR);
	//sobel算子平滑
	sobel();
	//非极大值抑制
  	nonMaxSupp();
  	//双阈值检测
  	threshold(low_thres, high_thres);
  	return thresImg;
}