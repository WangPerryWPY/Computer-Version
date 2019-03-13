#include "hough.h"

hough::hough() {

}

hough::~hough() {

}

hough::hough(CImg<float> picture) {
	Img = picture;
}

hough::hough(const char* filename) {
	Img.load(filename);
}

float hough::distance (float x, float y) {
	return sqrt(x * x + y * y);
}


bool sortCmp(const pair< pair<int, int>, int > &a, const pair< pair<int, int>, int > &b) {
	return a.second > b.second;
}

void hough::Hough_Statistics() {
	//double maxDistance = sqrt(Img._width*Img._width + Img._height*Img._height);
	double w = Img._width;
	double h = Img._height;
	double center_x = w/2;
	double center_y = h/2;
	double hough_h = ((sqrt(2.0) * (double)(h>w?h:w)) / 2.0);
	houghImg.resize(180, hough_h * 2, 1, 1, 0);
	cimg_forXY(cannyImg, x, y) {
		if (cannyImg(x,y) != 0) {
			cimg_forX(houghImg, angle) {
				double _angle = (double)PI*angle / 180.0f;
				int polar = (int)((((double)x - center_x)*cos(_angle) + ((double)y - center_y)*sin(_angle)) + hough_h);
				houghImg(angle, polar) += 1;
			}
		}
	}
}

void hough::GetLine() {
	resultImg = Img;
	//剔除掉可能出现的重合线，方法是取9x9空间内的霍夫最大值
	int hough_h = houghImg._height;
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
	sort(peaks.begin(), peaks.end(), sortCmp);
	for (int i = 0; lines.size() != 4; i++) {
		int angle = peaks[i].first.first;
		int polar = peaks[i].first.second;
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
		bool flag = true;
		for (int k = 0; k < lines.size(); k++) {
			if (distance((float)(lines[k].first.first - x1), (float)(lines[k].first.second - y1)) < 100 && distance((float)(lines[k].second.first - x2), (float)(lines[k].second.second - y2)) < 100) {
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
	//resultImg.display();
}

void hough::GetVertexs() {
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
	resultImg.display("houghImg");
}

vector< pair<int, int > > hough::process() {
	canny toCanny(Img);
	cannyImg = toCanny.process(1.0, 40, 60);
	Hough_Statistics();
	GetLine();
	GetVertexs();
	return vertex;
}