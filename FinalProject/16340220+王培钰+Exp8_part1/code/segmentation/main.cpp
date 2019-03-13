#include "ImageSegmentation.h"
int main () {
	for (int i = 9; i <= 9; i++) {
		CImg<float> img;
		string inpath = "./ImageData/";
		string txtpath = "./imagetxt/";
		string outpath = "./imageoutput/";
		stringstream ss;
    	ss << i;
    	string picname = ss.str() + ".jpg";
    	string txtname = txtpath + ss.str() + ".txt";
		string pathpic = inpath + picname;
		string outfile = outpath + "image" + ss.str();
		ofstream oFile;
	    oFile.open(txtname.c_str(), ios::out|ios::trunc);
	    oFile << picname << endl;
		oFile.close();
		img.load(pathpic.c_str());
		ImageSegmentation picture(img);
		picture.process(outfile, txtname.c_str());
	}
}