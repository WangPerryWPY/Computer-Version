#include "ImageSegmentation.h"

inline float distance (float x, float y) {
	return sqrt(x * x + y * y);
}

//Bilinear interpolation
/*inline void Interpolation(CImg<double> image, double x, double y, double P[]) {
	int x0 = floor(x);
	int x1 = (x0 < image.width() - 1) ? x0 + 1 : x0;
	int y0 = floor(y);
	int y1 = (y0 < image.height() - 1) ? y0 + 1 : y0;
	for (int c = 0; c < 3; c++) {
		double P1 = image(x0, y0, c) * (x1-x) + image(x1, y0, c) * (x-x0);
		double P2 = image(x0, y1, c) * (x1-x) + image(x1, y1, c) * (x-x0);	
		P[c] = P1 * (y1-y) + P2 * (y-y0);
	}
}*/


ImageSegmentation::ImageSegmentation() {
	threshold = 0;
}

ImageSegmentation::~ImageSegmentation() {

}

ImageSegmentation::ImageSegmentation(const char *filename) {
	threshold = 0;
	Img.load(filename);
}

void ImageSegmentation::rgb2gray() {
	grayImg.resize(Img._width, Img._height, 1, 1, 0);
	cimg_forXY(grayImg, x, y) {
		double R = Img(x,y,0);
		double G = Img(x,y,1);
		double B = Img(x,y,2);
		double Gray = (R * 299 + G * 587 + B * 114 + 500) / 1000;
		grayImg(x,y) = Gray;
	}
}

void ImageSegmentation::Gauss_blur() {
	blurImg = grayImg.get_blur(guassian_blur);
}

//迭代法求阈值
void ImageSegmentation::get_thres_iteration() {
	CImg<int> image = blurImg;
	CImg<int> hist = image.histogram(256, 0, 255);
	int size = blurImg.size();
	cimg_forX(hist, i) {
		threshold += i*hist(i);
	}
	threshold /= size;
	int threshold_new;
	while (true) {
		int t1 = 0, t2 = 0;
		int num1 = 0, num2 = 0;
		// 计算小于等于阈值threshold的灰度平均值t1以及大于阈值的t2
		cimg_forX(hist, i) {
			if (i <= threshold) {
				t1 += i * hist(i);
				num1 += hist(i);
			}
			else {
				t2 += i * hist(i);
				num2 += hist(i);
			}
		}
		if (num1 == 0 || num2 == 0)
			continue;
		t1 /= num1;
		t2 /= num2;
		threshold_new = (t1 + t2) / 2;
		// 若两个阈值相等，则返回阈值threshold，否则更新阈值继续循环
		if (threshold == threshold_new) break;
		else threshold = threshold_new;
	}
	cout << "threshold = " << threshold << endl;
}
//OSTU法求阈值
void ImageSegmentation::get_thres_ostu() {
	//定义类间方差
	double variance = 0.0;
	CImg<int> image = blurImg;
	CImg<int> hist = image.histogram(256, 0, 255);
	int size = blurImg.size();
	for (int i = 0; i < 256; i++) {
		//定义前景图，背景图的像素点所占比例以及平均灰度
		double p1 = 0.0, p2 = 0.0, g1 = 0.0, g2 = 0.0;
		cimg_forX(hist, j) {
			if (j <= i) {
				p1 += hist(j);
				g1 += j*hist(j);
			}
			else {
				p2 += hist(j);
				g2 += j*hist(j);
			}
		}
		if (p1 == 0 || p2 == 0)
			continue;
		g1 /= p1;
		p1 /= size;
		g2 /= p2;
		p2 /= size;
		double temp_variance = p1 * p2 * (g1 - g2) * (g1 - g2);
		if (variance < temp_variance) {
			variance = temp_variance;
			threshold = i;
		}
	}
	cout << "threshold = " << threshold << endl;
}

void ImageSegmentation::get_thres(string type_c) {
	if (type_c == "Iteration") {
		get_thres_iteration();
	}
	else {
		get_thres_ostu();
	}
}

void ImageSegmentation::Segmentation() {
	segImg.resize(Img._width, Img._height, 1, 1, 0);
	cimg_forXY(blurImg, x, y) {
		if (blurImg(x,y) > threshold) {
			segImg(x,y) = 0;
		}
		else
			segImg(x,y) = 255;
	}
	//segImg.display();
}

void ImageSegmentation::gradDection() {
	gradImg.resize(segImg._width, segImg._height, 1, 1, 0);
	CImg_3x3(I, double);
	cimg_for3x3(segImg, x, y, 0, 0, I, double) {
		const double ix = (Inn + 2 * Icn + Ipn) - (Ipp + 2 * Icp + Inp);
		const double iy = (Inp + 2 * Inc + Inn) - (Ipp + 2 * Ipc + Ipn);
		double grad = sqrt(ix * ix + iy * iy);
		if (grad > 255) grad = 255;
		if (grad < 0 ) grad = 0;
		gradImg(x, y) = grad;
	}
	gradImg.display();
}

void ImageSegmentation::Hough_Statistics() {
	//double maxDistance = sqrt(Img._width*Img._width + Img._height*Img._height);
	double w = Img._width;
	double h = Img._height;
	double center_x = w/2;
	double center_y = h/2;
	double hough_h = ((sqrt(2.0) * (double)(h>w?h:w)) / 2.0);
	houghImg.resize(180, hough_h * 2, 1, 1, 0);
	cimg_forXY(gradImg, x, y) {
		if (gradImg(x,y) != 0) {
			cimg_forX(houghImg, angle) {
				double _angle = (double)PI*angle / 180.0f;
				int polar = (int)((((double)x - center_x)*cos(_angle) + ((double)y - center_y)*sin(_angle)) + hough_h);
				//cout << polar << endl;
				houghImg(angle, polar) += 1;
			}
		}
	}
	//houghImg.display();
}

void ImageSegmentation::GetLine() {
	resultImg = Img;
	//剔除掉可能出现的重合线，方法是取9x9空间内的霍夫最大值
	int hough_h = houghImg._height;
	//int hough_w = houghImg._width;
	int img_h = Img._height;
	int img_w = Img._width;
	const int y_min = 0;
	const int y_max = Img._height - 1;
	const int x_min = 0;
	const int x_max = Img._width - 1;
	cimg_forXY(houghImg, angle, polar) {
		if (houghImg(angle, polar) >= Min_thres) {
			int max = houghImg(angle, polar);
			for(int ly=-DIFF;ly<=DIFF;ly++) {
				for(int lx=-DIFF;lx<=DIFF;lx++) {
					if( (ly+polar>=0 && ly+polar<houghImg._height) && (lx+angle>=0 && lx+angle<houghImg._width) ) {
						if( (int)houghImg(angle + lx, polar + ly ) > max ) {
							max = houghImg(angle + lx, polar + ly );
							ly = lx = DIFF + 1;
						}
					}
				}
			}
			if (max > (int)houghImg(angle, polar) )
				continue;
			peaks.push_back(pair< pair<int, int>, int >(pair<int, int>(angle, polar), houghImg(angle, polar)));
		}
	}
	sort(peaks.begin(), peaks.end(), [](const pair< pair<int, int>, int > &a, const pair< pair<int, int>, int > &b) -> int {return a.second > b.second ;});
	for (int i = 0; lines.size() != 4; i++) {
		int angle = peaks[i].first.first;
		int polar = peaks[i].first.second;
		//cout << angle << endl << polar << endl;
		int x1, y1, x2, y2;
		x1 = y1 = x2 = y2 = 0;
		double _angle = (double)PI*angle / 180.0f;
		if(angle >= 45 && angle <= 135) {
			x1 = 0;
			y1 = ((double)(polar-(hough_h/2)) - ((x1 - (img_w/2) ) * cos(_angle))) / sin(_angle) + (img_h / 2);
			x2 = img_w;
			y2 = ((double)(polar-(hough_h/2)) - ((x2 - (img_w/2) ) * cos(_angle))) / sin(_angle) + (img_h / 2);
		}
		else {
			y1 = 0;
			x1 = ((double)(polar-(hough_h/2)) - ((y1 - (img_h/2) ) * sin(_angle))) / cos(_angle) + (img_w / 2);
			y2 = img_h;
			x2 = ((double)(polar-(hough_h/2)) - ((y2 - (img_h/2) ) * sin(_angle))) / cos(_angle) + (img_w / 2);
		}
		//if
		bool flag = true;
		for (int k = 0; k < lines.size(); k++) {
			if (distance(lines[k].first.first - x1, lines[k].first.second - y1) < 100 && distance(lines[k].second.first - x2, lines[k].second.second - y2) < 100) {
				flag = false;
				break;
			}
		}
		if (flag == true) {
			lines.push_back(pair< pair<int, int>, pair<int, int> >(pair<int, int>(x1, y1), pair<int, int>(x2, y2)));
		}
	}
	for (int i = 0; i < lines.size(); i++) {
		cout << lines[i].first.first << ", " << lines[i].first.second << "  ..  " << lines[i].second.first << ", " << lines[i].second.second << endl;
		resultImg.draw_line(lines[i].first.first, lines[i].first.second, lines[i].second.first, lines[i].second.second, Red);
	}
	//resultImg.draw_line(200, 3458, 2500, 3459, Red);
	//resultImg.display();
}

void ImageSegmentation::GetVertexs() {
	for (int i = 0; i < lines.size(); i++) {
		double k0, b0;
		if (lines[i].first.first == lines[i].second.first) {
			k0 = DBL_MAX;
			b0 = lines[i].first.first;
		}
		else {
			k0 = (double) (lines[i].first.second - lines[i].second.second) / (lines[i].first.first - lines[i].second.first);
			b0 = (double) (lines[i].first.second * lines[i].second.first - lines[i].second.second * lines[i].first.first) / (lines[i].second.first - lines[i].first.first);
		}
		for (int j = i + 1; j < lines.size(); j++) {
			double k1, b1;
			if (lines[j].first.first == lines[j].second.first) {
				k1 = DBL_MAX;
				b1 = lines[j].first.first;
			}
			else {
				k1 = (double) (lines[j].first.second - lines[j].second.second) / (lines[j].first.first - lines[j].second.first);
				b1 = (double) (lines[j].first.second * lines[j].second.first - lines[j].second.second * lines[j].first.first) / (lines[j].second.first - lines[j].first.first);
			}
			if (k0 == k1)
				continue;
			if (k0 == DBL_MAX) {
				int _x = b0, _y = k1 * b0 + b1;
				if (_x >= 0 && _x < Img._width && _y >= 0 && _y < Img._height)
					vertex.push_back(make_pair(_x, _y));
				continue;
			}
			if (k1 == DBL_MAX) {
				int _x = b1, _y = k0 * b1 + b0;
				if (_x >= 0 && _x < Img._width && _y >= 0 && _y < Img._height)
					vertex.push_back(make_pair(_x, _y));
				continue;
			}
			int _x = (b0 - b1) / (k1 - k0);
			int _y = (k0 * b1 - k1 * b0) / (k0 - k1);
			if (_x >= 0 && _x < Img._width && _y >= 0 && _y < Img._height)
				vertex.push_back(make_pair(_x, _y));
		}
	}
	for (int i = 0; i < vertex.size(); i++) {
		cout << vertex[i].first << "  ...  " << vertex[i].second << endl;
		resultImg.draw_circle(vertex[i].first, vertex[i].second, 50, Red);
	}
	resultImg.display();
	resultImg.save("result_a.bmp");
}

void ImageSegmentation::orderVertexs() {
	sort(vertex.begin(), vertex.end(), [](const pair<int, int> &a, const pair<int, int> &b)-> int {return distance(a.first, a.second) < distance(b.first, b.second);});
	double w = distance(vertex[0].first - vertex[1].first, vertex[0].second - vertex[1].second);
	double h = distance(vertex[0].first - vertex[2].first, vertex[0].second - vertex[2].second);
	//纸张是横向的
	if (vertex[1].first < vertex[2].first && h > w) {
		swap(vertex[1], vertex[2]);
		swap(vertex[2], vertex[3]);
		vertex.push_back(vertex[0]);
		vertex.erase(vertex.begin());
	}
	//纸张是竖向的
	else {
		swap(vertex[2], vertex[3]);
	}
}

void ImageSegmentation::calcMatrix() {
  double x0 = vertex[0].first, x1 = vertex[1].first, x2 = vertex[2].first, x3 = vertex[3].first;
  double y0 = vertex[0].second, y1 = vertex[1].second, y2 = vertex[2].second, y3 = vertex[3].second;
  double dx3 = x0 - x1 + x2 - x3;
  double dy3 = y0 - y1 + y2 - y3;
  if (fabs(dx3) < 10e-5 && fabs(dy3) < 10e-5) {
    M[0] = x1 - x0, M[1] = y1 - y0, M[2] = 0;
    M[3] = x2 - x1, M[4] = y2 - y1, M[5] = 0;
    M[6] = x0, M[7] = y0, M[8] = 1;
  }
  else {
    double dx1 = x1 - x2, dx2 = x3 - x2, dy1 = y1 - y2, dy2 = y3 - y2;
    double det = dx1 * dy2 - dx2 * dy1;
    double a13 = (dx3 * dy2 - dx2 * dy3) / det;
    double a23 = (dx1 * dy3 - dx3 * dy1) / det;
    M[0] = x1 - x0 + a13 * x1, M[1] = y1 - y0 + a13 * y1, M[2] = a13;
    M[3] = x3 - x0 + a23 * x3, M[4] = y3 - y0 + a23 * y3, M[5] = a23;
    M[6] = x0, M[7] = y0, M[8] = 1;
  }
}

void ImageSegmentation::warping() {
	resultImg = Img;
	double P[3];
  	resultImg.resize(1050, 1485);  // 标准A4纸比例
  	double width = resultImg.width(), height = resultImg.height();
	cimg_forXY(resultImg, x, y) {
		double _x = x / width, _y = y / height;
		double denominator = M[2] * _x + M[5] * _y + M[8];
		double tx = (M[0] * _x + M[3] * _y + M[6]) / denominator;
		double ty = (M[1] * _x + M[4] * _y + M[7]) / denominator;
		/*Interpolation(Img, tx, ty, P);
		resultImg(x,y,0) = P[0];
		resultImg(x,y,1) = P[1];
		resultImg(x,y,2) = P[2];
		//cout << P[0] << endl;	*/
		cimg_forC(resultImg, c) {
			resultImg(x,y,c) = Img((int)tx, (int)ty, c);
		}
  	}
  	resultImg.display();
  	resultImg.save("result_b.bmp");
}

void ImageSegmentation::correct_process(string type_c) {
	rgb2gray();
	Gauss_blur();
	get_thres(type_c);
	Segmentation();
	gradDection();
	Hough_Statistics();
	GetLine();
	GetVertexs();
	orderVertexs();
	calcMatrix();
	warping();
}