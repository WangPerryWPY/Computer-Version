#include "Canny.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

int main() {
	//Canny img1;
	//img1.HoughImageDetect("Dataset2/1.bmp",11.0,0.5,0.85, 1, 150, 170);
	//Canny img2;
	//img2.HoughImageDetect("Dataset2/2.bmp",8.0,0.4,0.95, 4, 190, 230);
	//Canny img3;
	//img3.HoughImageDetect("Dataset2/3.bmp",8.0,0.4,0.95, 7, 140, 200);
	//Canny img4;
	//img4.HoughImageDetect("Dataset2/4.bmp",11.0,0.5,0.85, 3, 140, 210);
	//Canny img5;
	//img5.HoughImageDetect("Dataset2/5.bmp",11.0,0.5,0.85, 2, 460, 530);
	Canny img6;
	img6.HoughImageDetect("Dataset2/6.bmp", 5.0,0.25,0.75, 9, 40, 60);
}