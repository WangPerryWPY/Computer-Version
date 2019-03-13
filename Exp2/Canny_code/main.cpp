#include "pch.h"
#include "Canny.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

int main() {
	Canny img("lena.jpg","jpg");
	CImg<int> picture = img.canny_image(2.0,0.25,0.75);
	int count = 0;
	cimg_forXY(picture, x, y) {
		if (picture(x, y) == 0) {
			count++;
		}
	}
	cout << count << endl;
	CImg<int> picture_a = img.canny_line(picture, 10);
	int count1 = 0;
	cimg_forXY(picture_a, x, y) {
		if (picture_a(x, y) == 0) {
			count1++;
		}
	}
	cout << count1 << endl;
	CImg<int> picture_b = img.delete_line(picture_a);
	int count2 = 0;
	cimg_forXY(picture_b, x, y) {
		if (picture_b(x, y) == 0) {
			count2++;
		}
	}
	cout << count2 << endl;
	picture.display();
	picture.save("result.jpg");
	//CImg<int> picture_a = img.canny_line(picture, 10);
	//CImg<int> picture_b = img.delete_line(picture_a);
	picture_a.save("result_a.jpg");
	picture_a.display();
	picture_b.save("result_b.jpg");
	picture_b.display();
	
}