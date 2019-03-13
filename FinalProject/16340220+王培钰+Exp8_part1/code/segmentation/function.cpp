#include "function.h"
CImg<float> rgb2gray(CImg<float> picture) {
	CImg<float> graypic;
	graypic.resize(picture._width, picture._height, 1, 1, 0);
	cimg_forXY(graypic, x, y) {
		float R = picture(x,y,0);
		float G = picture(x,y,1);
		float B = picture(x,y,2);
		float Gray = (R * 299 + G * 587 + B * 114 + 500) / 1000;
		graypic(x,y) = Gray;
	}
	return graypic;
}

CImg<float> dilate_black(CImg<float> Pic, int flag) {
	//对图像进行膨胀操作
	CImg<float> picture = Pic;
	int k0 = flag / 2;
	cimg_for_insideXY(picture, x, y, k0) {
		int Count = 0;
		for (int i = -k0; i <= k0; i++) {
			for (int j = -k0; j <= k0; j++) {
				if (Pic(x+i, y+j) != 255) {
					Count++;
				}
			}
		}
		if (Count >= 2) {
			for (int i = -k0; i <= k0; i++) {
				for (int j = -k0; j <= k0; j++) {
					picture(x+i, y+j) = 0;
				}
			}
		}
		else {
			for (int i = -k0; i <= k0; i++) {
				for (int j = -k0; j <= k0; j++) {
					picture(x+i, y+j) = 255;
				}
			}
		}
	}
	return picture;
}

CImg<float> dilate_white(CImg<float> Pic, int flag) {
	//对图像进行膨胀操作
	CImg<float> picture = Pic;
	int k0 = flag / 2;
	cimg_for_insideXY(picture, x, y, k0) {
		int Count = 0;
		for (int i = -k0; i <= k0; i++) {
			for (int j = -k0; j <= k0; j++) {
				if (Pic(x+i, y+j) != 0) {
					Count++;
				}
			}
		}
		if (Count >= 2) {
			for (int i = -k0; i <= k0; i++) {
				for (int j = -k0; j <= k0; j++) {
					picture(x+i, y+j) = 255;
				}
			}
		}
		else {
			for (int i = -k0; i <= k0; i++) {
				for (int j = -k0; j <= k0; j++) {
					picture(x+i, y+j) = 0;
				}
			}
		}
	}
	return picture;
}

CImg<float> dilate_a(CImg<float> Pic) {
	CImg<float> picture = Pic;
	cimg_forXY(picture, x, y) {
		int sum = 0;
		if (x < 1 || x >= Pic._width - 1 || y < 1 || y >= Pic._height - 1) {
			continue;
		}
		if (Pic(x, y) == 0) {
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					if (((int)(Pic(x+i, y+j)) ^ 255) * filter[j+1][i+1] != 0)
						sum++;
				}
			}
			if (sum == 0) {
				picture(x, y) = 255;
			}
		}
	}
	return picture;
}

CImg<float> dilate_b(CImg<float> Pic) {
	CImg<float> picture = Pic;
	cimg_forXY(picture, x, y) {
		int sum = 0;
		if (x < 1 || x >= Pic._width - 1 || y < 1 || y >= Pic._height - 1) {
			continue;
		}
		if (Pic(x, y) != 0) {
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					if (((int)(Pic(x+i, y+j)) ^ 255) * filter[j+1][i+1] != 0)
						sum++;
				}
			}
			if (sum != 0) {
				picture(x, y) = 0;
			}
		}
	}
	return picture;
}

CImg<float> dilate_X(CImg<float> Pic) {
	CImg<float> picture = Pic;
	cimg_forXY(picture, x, y) {
		if (x < 1 || x >= Pic._width - 1) {
			continue;
		}
		if (Pic(x, y) != 0) {
			if ((Pic.atXY(x - 1, y) == 0) && (Pic.atXY(x + 1, y) == 0))
				picture(x, y) = 0;
		}
	}
	return picture;
}

CImg<float> dilate_Y(CImg<float> Pic) {
	CImg<float> picture = Pic;
	cimg_forXY(picture, x, y) {
		if (y < 2 || x >= Pic._height - 2) {
			continue;
		}
		if (Pic(x, y) != 0) {
			if ((Pic.atXY(x, y - 1) == 0 || Pic.atXY(x, y - 2) == 0) && (Pic.atXY(x, y + 1) == 0 || Pic.atXY(x, y + 2) == 0))
				picture(x, y) = 0;
		}
	}
	return picture;
}