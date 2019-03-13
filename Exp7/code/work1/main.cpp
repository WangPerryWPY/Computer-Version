#include "ImageSegmentation.h"
int main () {
	ImageSegmentation seg("Dataset/1.jpg");
	seg.correct_process("OSTU");
}