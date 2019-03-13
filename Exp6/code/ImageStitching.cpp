#include "ImageStitching.h"

float Interpolation(CImg<float> image, float x, float y, int channel);
bool IsBlack(CImg<float> img, int x, int y);
points HomographyMatrix(vector<point_pair> pair);
float MinX(CImg<float> src, points H);
float MinY(CImg<float> src, points H);
float MaxX(CImg<float> src, points H);
float MaxY(CImg<float> src, points H);
float getX(float x, float y, points H);
float getY(float x, float y, points H);

ImageStitching::ImageStitching() {
	struct dirent *ptr;
	DIR *dir;
	dir = opendir("TEST-ImageData2");
	while ((ptr = readdir(dir)) != NULL) {
		if(ptr->d_name[0] == '.')
            		continue;
		string file = string("TEST-ImageData2/") + string(ptr->d_name);
		cout << file << endl;
		const char *Ff = file.c_str();
		CImg<float> picture;
		picture.load(Ff);
		imgs.push_back(picture);
	}
}

ImageStitching::~ImageStitching() {
	imgs.clear();
}

CImg<float> ImageStitching::CylindricalProjection(CImg<float> pic) {
    int width = pic._width;
    int height = pic._height;
    int channel = pic._spectrum;
    CImg<float> result(width, height, 1, channel, 0);
    float  R = width / (2*tan(28.0f/2.0f*PI/180.0f));
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            float x0 = x - (float)(width / 2);
            float y0 = y - (float) (height/2);
            float _x = x0*sqrt(R*R + x0*x0) / R + (float) width / 2;
            float _y = y0*sqrt(R*R + x0*x0) / R + (float) height / 2;
            if (_x >= 0 && _x < width && _y >=0 && _y <= height) {
                for (int c = 0;  c < channel; c++) {
                    result(x,y,c) = Interpolation(pic,_x,_y,c);
                }
            }
        }
    }
    return result;
}

CImg<float> ImageStitching::convertTogray(CImg<float> pic) {
	CImg<float> picture(pic._width, pic._height, 1,1,0);
	cimg_forXY(pic,x,y) {
		float R = pic(x,y,0);
		float G = pic(x,y,1);
		float B = pic(x,y,2);		
		picture(x,y) =  R*0.299 + G*0.587 + B*0.114;
	}
	return picture;
}

map<vector<float>, VlSiftKeypoint> ImageStitching::SIFTFeatures(CImg<float> pic) {
    int noctaves = 4, nlevels = 2, o_min = 0;
    vl_sift_pix *ImageData = new vl_sift_pix[pic._width*pic._height];
    for (int i = 0; i < pic._width; i++) {
        for (int j = 0; j < pic._height; j++) {
            ImageData[j*pic._width + i] = pic(i,j,0);
        }
    }
    // 定义VlSiftFilt结构体指针
	VlSiftFilt *SiftFilt = NULL;
	// 创建一个新的sift滤波器
	SiftFilt = vl_sift_new(pic._width, pic._height, noctaves, nlevels, o_min);
    //int KeyPoint = 0;
    map<vector<float>, VlSiftKeypoint> Feature;
    if (vl_sift_process_first_octave(SiftFilt, ImageData) != VL_ERR_EOF) {
        while(true) {
            //计算每组中的关键点
            vl_sift_detect(SiftFilt);
			//遍历并绘制每个点  
			//KeyPoint += SiftFilt->nkeys;//检测到的关键点的数目
            VlSiftKeypoint *pKeyPoint = SiftFilt->keys;//检测到的关键点
			for (int i = 0; i<SiftFilt->nkeys; i++) {
                VlSiftKeypoint tempKeyPoint = *pKeyPoint;
                pKeyPoint++;
                double angles[4];
				int angleCount = vl_sift_calc_keypoint_orientations(SiftFilt, angles, &tempKeyPoint);//计算关键点的方向
				for (int j = 0; j<angleCount; j++) {
                    vector<float> Descriptor;
                    double tempAngle = angles[j];
					vl_sift_pix descriptors[128];
					// 计算每个方向的描述
					vl_sift_calc_keypoint_descriptor(SiftFilt, descriptors, &tempKeyPoint, tempAngle);
                    for (int k = 0; k < 128 ;k++) {
                        Descriptor.push_back(descriptors[k]);
                    }
                    Feature.insert(pair<vector<float>, VlSiftKeypoint>(Descriptor, tempKeyPoint));
                }
                
            }
            //下一阶  
			if (vl_sift_process_next_octave(SiftFilt) == VL_ERR_EOF)
			{
				break;
			}
			//free(pKeyPoint);  
			//KeyPoint = NULL;
        }
    }
    vl_sift_delete(SiftFilt);
	delete[]ImageData;
	ImageData = NULL;
    return Feature;
}

vector<point_pair> ImageStitching::KDtreeMatch(map<vector<float>, VlSiftKeypoint> feature_a, map<vector<float>, VlSiftKeypoint> feature_b) {
	vector<point_pair> result;

    VlKDForest* forest = vl_kdforest_new(VL_TYPE_FLOAT, 128, 1, VlDistanceL1);

	float *data = new float[128 * feature_a.size()];
	int k = 0;
	for (auto it = feature_a.begin(); it != feature_a.end(); it++) {
		const vector<float> &descriptors = it->first;
		for (int i = 0; i < 128; i++) {
			data[i + 128 * k] = descriptors[i];
		}
		k++;
	}

	vl_kdforest_build(forest, feature_a.size(), data);
	VlKDForestSearcher* searcher = vl_kdforest_new_searcher(forest);
	VlKDForestNeighbor neighbours[2];

	for (auto it = feature_b.begin(); it != feature_b.end(); it++){
		float *temp_data = new float[128];

		for (int i = 0; i < 128; i++) {
			temp_data[i] = (it->first)[i];
		}

		int nvisited = vl_kdforestsearcher_query(searcher, neighbours, 2, temp_data);

		float ratio = neighbours[0].distance / neighbours[1].distance;
		if (ratio < 0.5) {
			vector<float> des(128);
			for (int j = 0; j < 128; j++) {
				des[j] = data[j + neighbours[0].index * 128];
			}

			VlSiftKeypoint left = feature_a.find(des)->second;
			VlSiftKeypoint right = it->second;
			result.push_back(point_pair(left, right));
		}

		delete[] temp_data;
		temp_data = NULL;
	}

	vl_kdforestsearcher_delete(searcher);
	vl_kdforest_delete(forest);

	delete[] data;
	data = NULL;

	return result;
}

points ImageStitching::RANSAC(vector<point_pair> pairs) {
	srand(time(0));

	int iterations = ceil(log(1 - 0.99) / log(1 - pow(0.5, 4)));

	vector<int> max_inliner_indexs;

	while (iterations--) {
		vector<point_pair> random_pairs;
		set<int> seleted_indexs;

		for (int i = 0; i < 4; i++) {
			int index = rand() % pairs.size();
			while (seleted_indexs.find(index) != seleted_indexs.end()) {
				index = rand() % pairs.size();
			}
			seleted_indexs.insert(index);

			random_pairs.push_back(pairs[index]);
		}

		points H = HomographyMatrix(random_pairs);
		//cout << H.x1 << endl;
		vector<int> cur_inliner_indexs;
		for (int i = 0; i < pairs.size(); i++) {
			if (seleted_indexs.find(i) != seleted_indexs.end()) {
				continue;
			}

			float real_x = pairs[i].b.x;
			float real_y = pairs[i].b.y;

			float x = H.x1 * pairs[i].a.x + H.x2 * pairs[i].a.y + H.x3 * pairs[i].a.x * pairs[i].a.y + H.x4;
			float y = H.x5 * pairs[i].a.x + H.x6 * pairs[i].a.y + H.x7 * pairs[i].a.x * pairs[i].a.y + H.x8;

			float distance = sqrt((x - real_x) * (x - real_x) + (y - real_y) * (y - real_y));
			if (distance < 4) {
				cur_inliner_indexs.push_back(i);
			}
		}
		if (cur_inliner_indexs.size() > max_inliner_indexs.size()) {
			max_inliner_indexs = cur_inliner_indexs;
		}
	}

	int calc_size = max_inliner_indexs.size();

	CImg<double> A(4, calc_size, 1, 1, 0);
	CImg<double> b(1, calc_size, 1, 1, 0);

	for (int i = 0; i < calc_size; i++) {
		int cur_index = max_inliner_indexs[i];

		A(0, i) = pairs[cur_index].a.x;
		A(1, i) = pairs[cur_index].a.y;
		A(2, i) = pairs[cur_index].a.x * pairs[cur_index].a.y;
		A(3, i) = 1;

		b(0, i) = pairs[cur_index].b.x;
	}

	CImg<double> x1 = b.get_solve(A);

	for (int i = 0; i < calc_size; i++) {
		int cur_index = max_inliner_indexs[i];

		b(0, i) = pairs[cur_index].b.y;
	}

	CImg<double> x2 = b.get_solve(A);

	return points(x1(0, 0), x1(0, 1), x1(0, 2), x1(0, 3), x2(0, 0), x2(0, 1), x2(0, 2), x2(0, 3));
}

void ImageStitching::WarpTwoImg(CImg<float> src, CImg<float> &dst, points H, float offset_x, float offset_y) {
	for (int dst_x = 0; dst_x < dst.width(); dst_x++) {
		for (int dst_y = 0; dst_y < dst.height(); dst_y++) {
			int src_x = H.x1 * (dst_x + offset_x) + H.x2 * (dst_y + offset_y) + H.x3 * (dst_x + offset_x) * (dst_y + offset_y) + H.x4;
			int src_y = H.x5 * (dst_x + offset_x) + H.x6 * (dst_y + offset_y) + H.x7 * (dst_x + offset_x) * (dst_y + offset_y) + H.x8;

			if (src_x >= 0 && src_x < src.width() && src_y >= 0 && src_y < src.height()) {
				for (int k = 0; k < src.spectrum(); k++) {
					dst(dst_x, dst_y, k) = Interpolation(src, src_x, src_y, k);
				}
			}
		}
	}
}
void ImageStitching::MoveTwoImg(CImg<float> src, CImg<float> &dst, int offset_x, int offset_y) {
	for (int dst_x = 0; dst_x < dst.width(); dst_x++) {
		for (int dst_y = 0; dst_y < dst.height(); dst_y++) {
			int src_x = dst_x + offset_x;
			int src_y = dst_y + offset_y;

			if (src_x >= 0 && src_x < src.width() && src_y >= 0 && src_y < src.height()) {
				for (int k = 0; k < src.spectrum(); k++) {
					dst(dst_x, dst_y, k) = src(src_x, src_y, k);
				}
			}
		}
	}
}

CImg<float> ImageStitching::Blend(CImg<float> pic1, CImg<float> pic2) {
	double sum_a_x = 0;
	double sum_a_y = 0;
	int a_n = 0;

	double sum_overlap_x = 0;
	double sum_overlap_y = 0;
	int overlap_n = 0;
	if (pic1.width() > pic1.height()) {
		for (int x = 0; x < pic1.width(); x++) {
			if (!IsBlack(pic1, x, pic1.height() / 2)) {
				sum_a_x += x;
				a_n++;
			}

			if (!IsBlack(pic1, x, pic1.height() / 2) && !IsBlack(pic2, x, pic1.height() / 2)) {
				sum_overlap_x += x;
				overlap_n++;
			}
		}
	}
	else {
		for (int y = 0; y < pic1.height(); y++) {
			if (!IsBlack(pic1, pic1.width() / 2, y)) {
				sum_a_y += y;
				a_n++;
			}

			if (!IsBlack(pic1, pic1.width() / 2, y) && !IsBlack(pic2, pic2.width() / 2, y)) {
				sum_overlap_y += y;
				overlap_n++;
			}
		}
	}

	int min_len = (pic1.width() < pic1.height()) ? pic1.width() : pic1.height();

	int n_level = floor(log2(min_len));

	vector<CImg<float> > a_pyramid(n_level);
	vector<CImg<float> > b_pyramid(n_level);
	vector<CImg<float> > mask(n_level);

	// Initialize the base.
	a_pyramid[0] = pic1;
	b_pyramid[0] = pic2;
	mask[0] = CImg<float>(pic1.width(), pic1.height(), 1, 1, 0);

	if (pic1.width() > pic1.height()) {
		if (sum_a_x / a_n < sum_overlap_x / overlap_n) {
			for (int x = 0; x < sum_overlap_x / overlap_n; x++) {
				for (int y = 0; y < pic1.height(); y++) {
					mask[0](x, y) = 1;
				}
			}
		}
		else {
			for (int x = sum_overlap_x / overlap_n + 1; x < pic1.width(); x++) {
				for (int y = 0; y < pic1.height(); y++) {
					mask[0](x, y) = 1;
				}
			}
		}
	}
	else {
		if (sum_a_y / a_n < sum_overlap_y / overlap_n) {
			for (int x = 0; x < pic1.width(); x++) {
				for (int y = 0; y < sum_overlap_y / overlap_n; y++) {
					mask[0](x, y) = 1;
				}
			}
		}
		else {
			for (int x = 0; x < pic1.width(); x++) {
				for (int y = sum_overlap_y / overlap_n; y < pic1.height(); y++) {
					mask[0](x, y) = 1;
				}
			}
		}
	}

	// Down sampling a and b, building Gaussian pyramids.
	for (int i = 1; i < n_level; i++) {
		a_pyramid[i] = a_pyramid[i - 1].get_blur(2).get_resize(a_pyramid[i - 1].width() / 2, a_pyramid[i - 1].height() / 2, 1, a_pyramid[i - 1].spectrum(), 3);
		b_pyramid[i] = b_pyramid[i - 1].get_blur(2).get_resize(b_pyramid[i - 1].width() / 2, b_pyramid[i - 1].height() / 2, 1, b_pyramid[i - 1].spectrum(), 3);

		mask[i] = mask[i - 1].get_blur(2).get_resize(mask[i - 1].width() / 2, mask[i - 1].height() / 2, 1, mask[i - 1].spectrum(), 3);
	}

	// Building Laplacian pyramids.
	for (int i = 0; i < n_level - 1; i++) {
		a_pyramid[i] = a_pyramid[i] - a_pyramid[i + 1].get_resize(a_pyramid[i].width(), a_pyramid[i].height(), 1, a_pyramid[i].spectrum(), 3);
		b_pyramid[i] = b_pyramid[i] - b_pyramid[i + 1].get_resize(b_pyramid[i].width(), b_pyramid[i].height(), 1, b_pyramid[i].spectrum(), 3);
	}

	vector<CImg<float> > blend_pyramid(n_level);

	for (int i = 0; i < n_level; i++) {
		blend_pyramid[i] = CImg<float>(a_pyramid[i].width(), a_pyramid[i].height(), 1, a_pyramid[i].spectrum(), 0);
		cimg_forXYC(blend_pyramid[i], x, y, c) {
			blend_pyramid[i](x, y, c) = a_pyramid[i](x, y, c) * mask[i](x, y) + b_pyramid[i](x, y, c) * (1.0 - mask[i](x, y));
		}
	}

	CImg<float> res = blend_pyramid[n_level - 1];
	for (int i = n_level - 2; i >= 0; i--) {
		res.resize(blend_pyramid[i].width(), blend_pyramid[i].height(), 1, blend_pyramid[i].spectrum(), 3);
		cimg_forXYC(blend_pyramid[i], x, y, c) {
			res(x, y, c) = blend_pyramid[i](x, y, c) + res(x, y, c);
			if (res(x, y, c) > 255) res(x, y, c) = 255;
			else if (res(x, y, c) < 0) res(x, y, c) = 0;
		}
	}

	return res;
}

void ImageStitching::StitchProcess() {
	int size = imgs.size();
	// Save the features and corresponding coordinates of each image
	vector<map<vector<float>, VlSiftKeypoint>> features(size);
	for (int i = 0; i < size; i++) {
		imgs[i] = CylindricalProjection(imgs[i]);	
		CImg<float> gray = convertTogray(imgs[i]);
		features[i] = SIFTFeatures(gray);
		/*if (i == 0) {
			auto iter = features[i].begin();
			while (iter != features[i].end()) {
				cout << (iter->second).x << endl;
				iter++;	
			}
		}*/
		
	}
	bool need_stitching[STITCH_NUM][STITCH_NUM] = { false };
	//Record adjacent images for each image
	vector< vector<int> > matching_index(size);

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i == j)
				continue;

			vector<point_pair> pairs = KDtreeMatch(features[i], features[j]);
			if (pairs.size() >= 20) {
				need_stitching[i][j] = true;
				matching_index[i].push_back(j);
			}

		}
	}
	int one_side = 0;

	for (int i = 0; i < matching_index.size(); i++) {
		if (matching_index[i].size() == 1) {
			one_side = i;
			break;
		}
	}

	int start_index = one_side;
	int pre_start_index = -1;
	int n = matching_index.size() / 2;	
	while (n--) {
		for (int i = 0; i < matching_index[start_index].size(); i++) {
			if (matching_index[start_index][i] != pre_start_index) {
				pre_start_index = start_index;
				start_index = matching_index[start_index][i];

				break;
			}
		}
	}
	// Record the previous stitched picture
	int prev_dst_index = start_index;

	queue<int> unstitched_index;
	unstitched_index.push(start_index);

	resultImg = imgs[start_index];

	while (!unstitched_index.empty()) {
		int src_index = unstitched_index.front();
		unstitched_index.pop();

		for (int j = matching_index[src_index].size() - 1; j >= 0; j--) {
			int dst_index = matching_index[src_index][j];

			if (!need_stitching[src_index][dst_index]) {
				continue;
			}
			else {
				need_stitching[src_index][dst_index] = false;
				need_stitching[dst_index][src_index] = false;
				unstitched_index.push(dst_index);
			}

			// Feature matching using kd tree
			vector<point_pair> src_to_dst_pairs = KDtreeMatch(features[src_index], features[dst_index]);
			vector<point_pair> dst_to_src_pairs = KDtreeMatch(features[dst_index], features[src_index]);
			/*for (int aa = 0; aa < src_to_dst_pairs.size(); aa++) {
				cout << src_to_dst_pairs[aa].a.x << endl;
			}*/
			if (src_to_dst_pairs.size() > dst_to_src_pairs.size()) {
				dst_to_src_pairs.clear();
				for (int i = 0; i < src_to_dst_pairs.size(); i++) {
					point_pair temp(src_to_dst_pairs[i].b, src_to_dst_pairs[i].a);
					dst_to_src_pairs.push_back(temp);
				}
			}
			else {
				src_to_dst_pairs.clear();
				for (int i = 0; i < dst_to_src_pairs.size(); i++) {
					point_pair temp(dst_to_src_pairs[i].b, dst_to_src_pairs[i].a);
					src_to_dst_pairs.push_back(temp);
				}
			}
			/*for (int aa = 0; aa < src_to_dst_pairs.size(); aa++) {
				cout << src_to_dst_pairs[aa].a.x << endl;
			}*/
			// Using RANSAC algorithm to find homography matrix
			points forward_H = RANSAC(dst_to_src_pairs);
			points backward_H = RANSAC(src_to_dst_pairs);

			// Find the image size after stitching
			float min_x = MinX(imgs[dst_index], forward_H);
			min_x = (min_x < 0) ? min_x : 0;
			float min_y = MinY(imgs[dst_index], forward_H);
			min_y = (min_y < 0) ? min_y : 0;
			float max_x = MaxX(imgs[dst_index], forward_H);
			max_x = (max_x >= resultImg.width()) ? max_x : resultImg.width();
			float max_y = MaxY(imgs[dst_index], forward_H);
			max_y = (max_y >= resultImg.height()) ? max_y : resultImg.height();

			int new_width = ceil(max_x - min_x);
			int new_height = ceil(max_y - min_y);

			CImg<float> a(new_width, new_height, 1, imgs[dst_index].spectrum(), 0);
			CImg<float> b(new_width, new_height, 1, imgs[dst_index].spectrum(), 0);

			WarpTwoImg(imgs[dst_index], a, backward_H, min_x, min_y);
			MoveTwoImg(resultImg, b, min_x, min_y);

			// Update features based on homography matrix
			for (auto iter = features[dst_index].begin(); iter != features[dst_index].end(); iter++) {
				float cur_x = iter->second.x;
				float cur_y = iter->second.y;
				iter->second.x = getX(cur_x, cur_y, forward_H) - min_x;
				iter->second.y = getY(cur_x, cur_y, forward_H) - min_y;
				iter->second.ix = int(iter->second.x);
				iter->second.iy = int(iter->second.y);
			}
			// Update features based on offset
			for (auto iter = features[prev_dst_index].begin(); iter != features[prev_dst_index].end(); iter++) {
				iter->second.x -= min_x;
				iter->second.y -= min_y;
				iter->second.ix = int(iter->second.x);
				iter->second.iy = int(iter->second.y);
			}
			// Image Blend
			resultImg = Blend(a, b);
			prev_dst_index = dst_index;
		}
	}
	resultImg.save("result/result2.jpg");
}

//Bilinear interpolation
float Interpolation(CImg<float> image, float x, float y, int channel) {
	 int x0 = floor(x);
	  int x1 = (x0 < image.width() - 1) ? x0 + 1 : x0;

	  int y0 = floor(y);
	  int y1 = (y0 < image.height() - 1) ? y0 + 1 : y0;

	  float P1 = image(x0, y0, channel) * (x1-x) + image(x1, y0, channel) * (x-x0);
	  float P2 = image(x0, y1, channel) * (x1-x) + image(x1, y1, channel) * (x-x0);	
	  float P = P1 * (y1-y) + P2 * (y-y0);
	  return P;
}

points HomographyMatrix(vector<point_pair> pair) {
	
	float x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8;
	
	float u0 = pair[0].a.x, v0 = pair[0].a.y;
	float u1 = pair[1].a.x, v1 = pair[1].a.y;
	float u2 = pair[2].a.x, v2 = pair[2].a.y;
	float u3 = pair[3].a.x, v3 = pair[3].a.y;

	float x0 = pair[0].b.x, y0 = pair[0].b.y;
	float x1 = pair[1].b.x, y1 = pair[1].b.y;
	float x2 = pair[2].b.x, y2 = pair[2].b.y;
	float x3 = pair[3].b.x, y3 = pair[3].b.y;

	x_1 = -(u0*v0*v1*x2 - u0*v0*v2*x1 - u0*v0*v1*x3 + u0*v0*v3*x1 - u1*v0*v1*x2 + u1*v1*v2*x0 + u0*v0*v2*x3 - u0*v0*v3*x2 + u1*v0*v1*x3 - u1*v1*v3*x0 + u2*v0*v2*x1 - u2*v1*v2*x0
		- u1*v1*v2*x3 + u1*v1*v3*x2 - u2*v0*v2*x3 + u2*v2*v3*x0 - u3*v0*v3*x1 + u3*v1*v3*x0 + u2*v1*v2*x3 - u2*v2*v3*x1 + u3*v0*v3*x2 - u3*v2*v3*x0 - u3*v1*v3*x2 + u3*v2*v3*x1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	x_2 = (u0*u1*v0*x2 - u0*u2*v0*x1 - u0*u1*v0*x3 - u0*u1*v1*x2 + u0*u3*v0*x1 + u1*u2*v1*x0 + u0*u1*v1*x3 + u0*u2*v0*x3 + u0*u2*v2*x1 - u0*u3*v0*x2 - u1*u2*v2*x0 - u1*u3*v1*x0
		- u0*u2*v2*x3 - u0*u3*v3*x1 - u1*u2*v1*x3 + u1*u3*v1*x2 + u1*u3*v3*x0 + u2*u3*v2*x0 + u0*u3*v3*x2 + u1*u2*v2*x3 - u2*u3*v2*x1 - u2*u3*v3*x0 - u1*u3*v3*x2 + u2*u3*v3*x1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	x_3 = (u0*v1*x2 - u0*v2*x1 - u1*v0*x2 + u1*v2*x0 + u2*v0*x1 - u2*v1*x0 - u0*v1*x3 + u0*v3*x1 + u1*v0*x3 - u1*v3*x0 - u3*v0*x1 + u3*v1*x0
		+ u0*v2*x3 - u0*v3*x2 - u2*v0*x3 + u2*v3*x0 + u3*v0*x2 - u3*v2*x0 - u1*v2*x3 + u1*v3*x2 + u2*v1*x3 - u2*v3*x1 - u3*v1*x2 + u3*v2*x1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	x_4 = (u0*u1*v0*v2*x3 - u0*u1*v0*v3*x2 - u0*u2*v0*v1*x3 + u0*u2*v0*v3*x1 + u0*u3*v0*v1*x2 - u0*u3*v0*v2*x1 - u0*u1*v1*v2*x3 + u0*u1*v1*v3*x2 + u1*u2*v0*v1*x3 - u1*u2*v1*v3*x0 - u1*u3*v0*v1*x2 + u1*u3*v1*v2*x0
		+ u0*u2*v1*v2*x3 - u0*u2*v2*v3*x1 - u1*u2*v0*v2*x3 + u1*u2*v2*v3*x0 + u2*u3*v0*v2*x1 - u2*u3*v1*v2*x0 - u0*u3*v1*v3*x2 + u0*u3*v2*v3*x1 + u1*u3*v0*v3*x2 - u1*u3*v2*v3*x0 - u2*u3*v0*v3*x1 + u2*u3*v1*v3*x0)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	x_5 = -(u0*v0*v1*y2 - u0*v0*v2*y1 - u0*v0*v1*y3 + u0*v0*v3*y1 - u1*v0*v1*y2 + u1*v1*v2*y0 + u0*v0*v2*y3 - u0*v0*v3*y2 + u1*v0*v1*y3 - u1*v1*v3*y0 + u2*v0*v2*y1 - u2*v1*v2*y0
		- u1*v1*v2*y3 + u1*v1*v3*y2 - u2*v0*v2*y3 + u2*v2*v3*y0 - u3*v0*v3*y1 + u3*v1*v3*y0 + u2*v1*v2*y3 - u2*v2*v3*y1 + u3*v0*v3*y2 - u3*v2*v3*y0 - u3*v1*v3*y2 + u3*v2*v3*y1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	x_6 = (u0*u1*v0*y2 - u0*u2*v0*y1 - u0*u1*v0*y3 - u0*u1*v1*y2 + u0*u3*v0*y1 + u1*u2*v1*y0 + u0*u1*v1*y3 + u0*u2*v0*y3 + u0*u2*v2*y1 - u0*u3*v0*y2 - u1*u2*v2*y0 - u1*u3*v1*y0
		- u0*u2*v2*y3 - u0*u3*v3*y1 - u1*u2*v1*y3 + u1*u3*v1*y2 + u1*u3*v3*y0 + u2*u3*v2*y0 + u0*u3*v3*y2 + u1*u2*v2*y3 - u2*u3*v2*y1 - u2*u3*v3*y0 - u1*u3*v3*y2 + u2*u3*v3*y1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	x_7 = (u0*v1*y2 - u0*v2*y1 - u1*v0*y2 + u1*v2*y0 + u2*v0*y1 - u2*v1*y0 - u0*v1*y3 + u0*v3*y1 + u1*v0*y3 - u1*v3*y0 - u3*v0*y1 + u3*v1*y0
		+ u0*v2*y3 - u0*v3*y2 - u2*v0*y3 + u2*v3*y0 + u3*v0*y2 - u3*v2*y0 - u1*v2*y3 + u1*v3*y2 + u2*v1*y3 - u2*v3*y1 - u3*v1*y2 + u3*v2*y1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	x_8 = (u0*u1*v0*v2*y3 - u0*u1*v0*v3*y2 - u0*u2*v0*v1*y3 + u0*u2*v0*v3*y1 + u0*u3*v0*v1*y2 - u0*u3*v0*v2*y1 - u0*u1*v1*v2*y3 + u0*u1*v1*v3*y2 + u1*u2*v0*v1*y3 - u1*u2*v1*v3*y0 - u1*u3*v0*v1*y2 + u1*u3*v1*v2*y0
		+ u0*u2*v1*v2*y3 - u0*u2*v2*v3*y1 - u1*u2*v0*v2*y3 + u1*u2*v2*v3*y0 + u2*u3*v0*v2*y1 - u2*u3*v1*v2*y0 - u0*u3*v1*v3*y2 + u0*u3*v2*v3*y1 + u1*u3*v0*v3*y2 - u1*u3*v2*v3*y0 - u2*u3*v0*v3*y1 + u2*u3*v1*v3*y0)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	return points(x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8);
}

float getX(float x, float y, points H) {
	return H.x1 * x + H.x2 * y + H.x3 * x * y + H.x4;
}

float getY(float x, float y, points H) {
	return H.x5 * x + H.x6 * y + H.x7 * x * y + H.x8;
}


float MaxX(CImg<float> src, points H) {
	float max_x = getX(0, 0, H);

	if (getX(src.width() - 1, 0, H) > max_x) {
		max_x = getX(src.width() - 1, 0, H);
	}
	if (getX(0, src.height() - 1, H) > max_x) {
		max_x = getX(0, src.height() - 1, H);
	}
	if (getX(src.width() - 1, src.height() - 1, H) > max_x) {
		max_x = getX(src.width() - 1, src.height() - 1, H);
	}

	return max_x;
}

float MinX(CImg<float> src, points H) {
	float min_x = getX(0, 0, H);

	if (getX(src.width() - 1, 0, H) < min_x) {
		min_x = getX(src.width() - 1, 0, H);
	}
	if (getX(0, src.height() - 1, H) < min_x) {
		min_x = getX(0, src.height() - 1, H);
	}
	if (getX(src.width() - 1, src.height() - 1, H) < min_x) {
		min_x = getX(src.width() - 1, src.height() - 1, H);
	}

	return min_x;
}

float MaxY(CImg<float> src, points H) {
	float max_y = getY(0, 0, H);

	if (getY(src.width() - 1, 0, H) > max_y) {
		max_y = getY(src.width() - 1, 0, H);
	}
	if (getY(0, src.height() - 1, H) > max_y) {
		max_y = getY(0, src.height() - 1, H);
	}
	if (getY(src.width() - 1, src.height() - 1, H) > max_y) {
		max_y = getY(src.width() - 1, src.height() - 1, H);
	}

	return max_y;
}

float MinY(CImg<float> src, points H) {
	float min_y = getY(0, 0, H);

	if (getY(src.width() - 1, 0, H) < min_y) {
		min_y = getY(src.width() - 1, 0, H);
	}
	if (getY(0, src.height() - 1, H) < min_y) {
		min_y = getY(0, src.height() - 1, H);
	}
	if (getY(src.width() - 1, src.height() - 1, H) < min_y) {
		min_y = getY(src.width() - 1, src.height() - 1, H);
	}

	return min_y;
}

bool IsBlack(CImg<float> img, int x, int y) {
	if (img(x, y, 0) == 0 && img(x, y, 1) == 0 && img(x, y, 2) == 0)
		return true;
	return false;
}