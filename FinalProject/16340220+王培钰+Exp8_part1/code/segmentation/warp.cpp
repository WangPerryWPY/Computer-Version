#include "warp.h"

warp::warp(){}

warp::warp(const char *filename) {
	Img.load(filename);
}

warp::warp(CImg<float> picture) {
	Img = picture;
}
warp::~warp(){}

float warp::distance (float x, float y) {
	return sqrt(x * x + y * y);
}

bool sortCmp(const pair<int, int> &a, const pair<int, int> &b) {
	return sqrt(a.first * a.first + a.second * a.second) < sqrt(b.first * b.first + b.second * b.second);
}

void warp::orderVertexs() {
	sort(vertex.begin(), vertex.end(), sortCmp);
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

void warp::getSize() {
	double w1 = distance(vertex[0].first - vertex[1].first, vertex[0].second - vertex[1].second);
	double h1 = distance(vertex[1].first - vertex[2].first, vertex[1].second - vertex[2].second);
	double w2 = distance(vertex[2].first - vertex[3].first, vertex[2].second - vertex[3].second);
	double h2 = distance(vertex[0].first - vertex[3].first, vertex[0].second - vertex[3].second);
	paper_width = w1 > w2 ? w2 : w1;
	paper_height = h1 > h2 ? h2 : h1;
	if (paper_width > 595) {
		paper_width = 595;
	}
	if (paper_height > 842) {
		paper_height = 842;
	}
}

void warp::calcMatrix() {
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

void warp::warping() {
	double P[3];
  	resultImg.resize(paper_width, paper_height); 
  	double width = resultImg.width(), height = resultImg.height();
	cimg_forXY(resultImg, x, y) {
		double _x = x / width, _y = y / height;
		double denominator = M[2] * _x + M[5] * _y + M[8];
		double tx = (M[0] * _x + M[3] * _y + M[6]) / denominator;
		double ty = (M[1] * _x + M[4] * _y + M[7]) / denominator;
		cimg_forC(resultImg, c) {
			resultImg(x,y,c) = Img((int)tx, (int)ty, c);
		}
  	}
  	resultImg.display("A4Img");
}

CImg<float> warp::process(const char *txtname) {
	hough toHough(Img);
	vertex = toHough.process();
	orderVertexs();
	ofstream oFile;
    oFile.open(txtname, ios::app);
    int i = 0;
    for (i = 0; i < vertex.size() - 1; i++) {
    	stringstream s0, s1;
    	s0 << vertex[i].first;
    	s1 << vertex[i].second;
    	oFile << "(" << s0.str() << ", " << s1.str() << ")" << endl;
    }
    stringstream s0, s1;
    s0 << vertex[i].first;
    s1 << vertex[i].second;
    oFile << "(" << s0.str() << ", " << s1.str() << ")";
	oFile.close();
	getSize();
	calcMatrix();
	warping();
	return resultImg;
}