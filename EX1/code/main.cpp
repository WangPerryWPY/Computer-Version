#include "ex.hpp"
using namespace std;
int main(int argc, char const *argv[])
{
	Test pic;
	Test pic1;
	pic.change();
	pic.DrawCircle_blue1();
	pic.DrawCircle_yellow1();
	pic.DrawLine1();
	pic.Todisplay();
	CImg<unsigned char> temp = pic.getSrcImg();
	temp.save("2.bmp");
	pic1.change();
	pic1.DrawCircle_blue2();
	pic1.DrawCircle_yellow2();
	pic1.DrawLine2();
	pic1.Todisplay();
	CImg<unsigned char> temp1 = pic1.getSrcImg();
	temp1.save("2.1.bmp");
	return 0;
}