#include <stdio.h>
#include <math.h>
#include "Hough.h"
#include <string>
#include <vector>
#include "CImg-2.4.0_pre090618\CImg.h"
#include <iostream>
#include <algorithm>
#include <float.h>
#include <sstream> 
#define minRadius 300
using namespace std;
using namespace cimg_library;

Hough::Hough() { 
}

void Hough::LineDetect_Image(string name, int scale, double sigma, int thresh_low, int thresh_high, int lineNum, int linedisturb, int distance) {
  const char *names = name.c_str();
  img.load_bmp(names);
  toGrayScale();
  useFilter(scale, sigma);
  sobel();
  nonMaxSupp();
  threshold(thresh_low, thresh_high);
  Polar_Line();
  CImg<unsigned char> pic = LineDetect(lineNum, linedisturb, distance);
  CImg<unsigned char> pic1 = drawPoint();
  string path = "result1/" + name.substr(9,1) + "_a.bmp";
  string path1 = "result1/" + name.substr(9,1) + "_b.bmp";
  pic.save(path.c_str());
  pic1.save(path1.c_str());
  /*string path = "4.bmp";
  string pathAfter;
  pathAfter = path.substr(0, path.length() - 4) + "_threshold.bmp";
  thres.save(pathAfter.c_str());*/
}

void Hough::CircleDetect_Image(string name, int scale, double sigma, int thresh_low, int thresh_high, int Nums, int min_r, int max_r) {
  const char *names = name.c_str();
  img.load_bmp(names);
  toGrayScale();
  useFilter(scale, sigma);
  sobel();
  nonMaxSupp();
  threshold(thresh_low, thresh_high);
  CircleDetect(Nums, min_r, max_r);
  findpixelCircle();
  /*stringstream s, ss;
  string s1, s2;
  s << thresh_low;
  s >> s1;
  ss << thresh_high;
  ss >> s2;*/
  string path = "result2/" + name.substr(9,1) + "_a.bmp";
  string path1 = "result2/" + name.substr(9,1) + "_b.bmp";
  thres_img.save(path.c_str());
  CircleImg.save(path1.c_str());
}

vector<vector<double> > Hough::createFilter(int scale, double sigma) {
  vector<vector<double> > filter;

  // initialize
  int s = 2 * scale + 1;
  for (int i = 0; i < s; ++i) {
    vector<double> temp;
    for (int j = 0; j < s; j++) {
      temp.push_back(-1);
    }
    filter.push_back(temp);
  }

  float sum = 0;
  for (int x = -s / 2; x <= s / 2; ++x) {
    for (int y = -s / 2; y <= s / 2; ++y) {
      filter[x + s / 2][y + s / 2] = (exp(-(x * x + y * y) / (2.0 * sigma * sigma))) * (1 / (M_PI * (2.0 * sigma * sigma)));
      sum += filter[x + s / 2][y + s / 2];
    }
  }

  // normalize the Filter
  for (int i = 0; i < s; ++i) {
    for (int j = 0; j < s; ++j) {
      filter[i][j] /= sum;
    }
  }

  return filter;
}

vector<vector<double> > Hough::createSobelFilterX() {
  // Sobel x filter
  vector<vector<double> > x_filter(3);
  x_filter[0] = { -1.0, 0, 1.0 };
  x_filter[1] = { -2.0, 0, 2.0 };
  x_filter[2] = { -1.0, 0, 1.0 };
  return x_filter;
}

vector<vector<double> > Hough::createSobelFilterY() {
  // Sobel y filter
  vector<vector<double> > y_filter(3);
  y_filter[0] = { 1.0, 2.0, 1.0 };
  y_filter[1] = { 0, 0, 0 };
  y_filter[2] = { -1.0, -2.0, -1.0 };
  return y_filter;
}

void Hough::toGrayScale()
{
  grayscaled = CImg<unsigned char>(img.width(), img.height(), 1);
  cimg_forXY(img, x, y) {
    unsigned char r = img(x, y, 0);
    unsigned char g = img(x, y, 1);
    unsigned char b = img(x, y, 2);
    double gray = (r * 0.2126 + g * 0.7152 + b * 0.0722);
    grayscaled.atXY(x, y) = gray;
  }
}

void Hough::useFilter(int scale, double sigma) {
  vector<vector<double> > filter = createFilter(scale, sigma);
  int size = (int)filter.size() / 2;
  gFiltered = CImg<unsigned char>(grayscaled.width() - 2 * size, grayscaled.height() - 2 * size, 1);
  for (int i = size; i < grayscaled.height() - size; ++i) {
    for (int j = size; j < grayscaled.width() - size; ++j) {
      double sum = 0;
      for (int x = 0; x < filter.size(); ++x) {
        for (int y = 0; y < filter.size(); ++y) {
          sum += filter[x][y] * (double)(grayscaled.atXY(j + y - size, i + x - size));
        }
      }
      gFiltered.atXY(j - size, i - size) = sum;
    }
  }
}

void Hough::sobel() {
  // Sobel filter
  vector<vector<double> > x_filter = createSobelFilterX();
  vector<vector<double> > y_filter = createSobelFilterY();

  int size = (int)x_filter.size() / 2;  // limit size
  const int width = gFiltered.width() - 2 * size;
  const int height = gFiltered.height() - 2 * size;
  sFiltered = CImg<unsigned char>(width, height, 1);
  angles = CImg<float>(width, height, 1);

  for (int i = size; i < gFiltered.height() - size; ++i) {
    for (int j = size; j < gFiltered.width() - size; ++j) {
      double sumx = 0;
      double sumy = 0;
      for (int x = 0; x < x_filter.size(); ++x) {
        for (int y = 0; y < y_filter.size(); ++y) {
          sumx += x_filter[x][y] * (double)(gFiltered.atXY(j + y - size, i + x - size));
          sumy += y_filter[x][y] * (double)(gFiltered.atXY(j + y - size, i + x - size));
        }
      }
      double sumxsq = sumx * sumx;
      double sumysq = sumy * sumy;
      double sq2 = sqrt(sumxsq + sumysq);
      if (sq2 > 255) sq2 = 255;
      sFiltered.atXY(j - size, i - size) = sq2;
      if (sumx == 0) {
        angles.atXY(j - size, i - size) = 90;
      }    
      else {
        angles.atXY(j - size, i - size) = atan(sumy / sumx);
      }   
    }
  }
}

void Hough::nonMaxSupp() {
  non = CImg<unsigned char>(sFiltered.width() - 2, sFiltered.height() - 2, 1);
  for (int i = 1; i < sFiltered.height() - 1; ++i) {
    for (int j = 1; j < sFiltered.width() - 1; ++j) {
      float tangent = angles.atXY(j, i);
      non.atXY(j - 1, i - 1) = sFiltered.atXY(j, i);
      // horizontal edge
      if (((-22.5 < tangent) && (tangent <= 22.5)) || ((157.5 < tangent) && (tangent <= -157.5))) {
        if ((sFiltered.atXY(j, i) < sFiltered.atXY(j + 1, i)) || (sFiltered.atXY(j, i) < sFiltered.atXY(j - 1, i)))
          non.atXY(j - 1, i - 1) = 0;
      }
      // vertical edge
      if (((-112.5 < tangent) && (tangent <= -67.5)) || ((67.5 < tangent) && (tangent <= 112.5))) {
        if ((sFiltered.atXY(j, i) < sFiltered.atXY(j, i + 1)) || (sFiltered.atXY(j, i) < sFiltered.atXY(j, i - 1)))
          non.atXY(j - 1, i - 1) = 0;
      }
      // -45 degree edge
      if (((-67.5 < tangent) && (tangent <= -22.5)) || ((112.5 < tangent) && (tangent <= 157.5))) {
        if ((sFiltered.atXY(j, i) < sFiltered.atXY(j + 1, i - 1)) || (sFiltered.atXY(j, i) < sFiltered.atXY(j - 1, i + 1)))
          non.atXY(j - 1, i - 1) = 0;
      }
      // 45 degree edge
      if (((-157.5 < tangent) && (tangent <= -112.5)) || ((22.5 < tangent) && (tangent <= 67.5))) {
        if ((sFiltered.atXY(j, i) < sFiltered.atXY(j + 1, i + 1)) || (sFiltered.atXY(j, i) < sFiltered.atXY(j - 1, i - 1)))
          non.atXY(j - 1, i - 1) = 0;
      }
    }
  }
}

void Hough::threshold(int low, int high) {
  thres = CImg<unsigned char>(non);
  for (int i = 0; i < non.height(); ++i) {
    for (int j = 0; j < non.width(); ++j) {
      //cout << (int)non(j,i) << " ";
      if (thres.atXY(j, i) > high) {
        thres.atXY(j, i) = 255;
      } else if (thres.atXY(j, i) < low) {
        thres.atXY(j, i) = 0;
      } else {
        bool anyHigh = false;
        bool anyBetween = false;
        for (int x = i - 1; x < i + 2; ++x) {
          for (int y = j - 1; y < j + 2; ++y) {
            if (x <= 0 || y <= 0 || x > thres.height() || y > thres.width()) {
              continue;
            } else {
              if (thres.atXY(y, x) > high) {
                thres.atXY(j, i) = 255;
                anyHigh = true;
                break;
              } else if (thres.atXY(y, x) <= high && thres.atXY(y, x) >= low) {
                anyBetween = true;
              }               
            }
          }
          if (anyHigh) break;
        }
        if (!anyHigh && anyBetween) {
          for (int x = i - 2; x < i + 3; ++x) {
            for (int y = j - 1; y < j + 3; ++y) {
              if (x < 0 || y < 0 || x > thres.height() || y > thres.width()) {
                continue;
              } else {
                if (thres.atXY(y, x) > high) {
                  thres.atXY(j, i) = 255;
                  anyHigh = true;
                  break;
                }
              }
            }
            if (anyHigh) {
              break;
            }
          }
        }
        if (!anyHigh) {
          thres.atXY(j, i) = 0;
        }
      }
    }
  }
}

void Hough::Polar_Line() {
  CImg<int> pic = thres;
  int rows = pic.width();
  int cols = pic.height();
  int maxLength = (int)sqrt(pow(rows, 2) + pow(cols, 2));
  HoughImg.resize(360, maxLength, 1, 1, 0);
  HoughImg.fill(0);
  cimg_forXY(pic, x, y) {
    if (pic(x,y) != 0) {
      for (int i = 0; i < 360; i++) {
        int r = (int)x * cos(M_PI*i / 180) + y * sin(M_PI*i / 180);
        if (r >= 0 && r < maxLength) {
          HoughImg(i, r)++;
        }
      }
    }
  }
}

CImg<unsigned char> Hough::LineDetect(int lineNum, int linedisturb, int distance) {
  int linenum = lineNum;
  const int Num = lineNum;
  //int linenum = lineNum;
  //HoughImg.display();
  CImg<unsigned char> pic = thres;
   
  pic.resize(pic.width(),pic.height(),1,3);
  vector<parameter> Lines;
  vector<parameter> result;
  vector<int> lineweight;
  cimg_forXY(HoughImg, x, y) {
    parameter a;
    a.x = x;
    a.y = y;
    Lines.push_back(a);
    lineweight.push_back(HoughImg(x, y));
  }
  vector<int> sortlineweight = lineweight;
  sort(sortlineweight.begin(), sortlineweight.end(), greater<int>());
  int fnums = lineNum + linedisturb;
  for (int i = 0; i < fnums; i++) {
    int weight = sortlineweight[i], index;
    vector<int>::iterator iter = find(lineweight.begin(), lineweight.end(), weight);
    index = iter - lineweight.begin();
    int x1 = Lines[index].x, y1 = Lines[index].y;
    //cout << x1 << " " << y1 << endl;
    bool flag = 1;
    for (int k = 0; k < i; k++) {
        int x0 = result[k].x;
        int y0 = result[k].y;
        if (sqrt(pow(x1-x0,2)+pow(y1-y0,2)) < distance) {
          flag = 0;
          break;
        }
    }
    if (flag == 1) {
      result.push_back(Lines[index]);
    }
    else {
      fnums++;
    }
  }
  const unsigned char red[] = { 255, 0, 0 };
  //cout << result.size() << endl;
  
  for (int i = 0; i < result.size(); i++) {
    int a0 = 0;
    nums.push_back(a0);
  }
  for (int i = 0; i < result.size(); i++) {
    //cout << result[i].x << " " << result[i].y << endl;
    //此时sin角度值为0，k为无穷，单独讨论
    if (result[i].x == 0 || result[i].x == 180) {
      int r = 0;
      if (result[i].x == 0) {
        r = result[i].y;
        //cout << "Line " << i+1 << ": x=" << r << endl;
      }
      else {
        r  -= (result[i].y);
        //cout << "Line " << i+1 << ": x=" << r << endl;
      }
      //const int x_min = 0;
      //const int x_max = pic.width() - 1;
      const int y_min = 0;
      const int y_max = pic.height() - 1;
      lineParameter temp;
      temp.k = DBL_MAX;
      temp.b = r;
      Hough_Line.push_back(temp);
      for (int yi = y_min; yi < y_max; yi++) {
        if (thres.atXY(r, yi) != 0) {
          nums[i]++;
        }
      }
      //pic.draw_line(r , y_min, r, y_max, red);
    }

    else {
      double theta = double(result[i].x)*M_PI/180;
      int r = result[i].y;
      //cout << "Line " << i << "  " << theta << "  " << r << endl;
      double k = (double)(-cos(theta) / sin(theta));
      double b = (double) r / sin(theta);
      lineParameter temp;
      temp.k = k;
      temp.b = b;
      Hough_Line.push_back(temp);
      const int x_min = 0;
      const int x_max = pic.width() - 1;
      const int y_min = 0;
      const int y_max = pic.height() - 1;
      const int x0 = (double)(y_min - b) / k;
      const int x1 = (double)(y_max - b) / k;
      const int y0 = x_min * k + b;
      const int y1 = x_max * k + b;
      if (abs(k) > 1) {
        for (int yi = y1; yi < y_max; yi++) {
          int xi = (double)(yi - b) / k;
          if (thres.atXY(xi,yi) != 0 || thres.atXY(xi+1,yi) != 0 || thres.atXY(xi,yi+1) != 0 || thres.atXY(xi+1,yi+1) != 0)
            nums[i]++;
        }
        //pic.draw_line(x0, y_min, x1, y_max, red);
      }
      else {
        for (int xi = x_min; xi < x_max; xi++) {
          int yi = k*xi + b;
          if (thres.atXY(xi,yi) != 0 || thres.atXY(xi+1,yi) != 0 || thres.atXY(xi,yi+1) != 0 || thres.atXY(xi+1,yi+1) != 0)
            nums[i]++;
        }
        //pic.draw_line(x_min, y0, x_max, y1, red);
      }
    }
    //cout << nums[i] << endl;
  }
  vector<int> sortnums = nums;
  sort(sortnums.begin(), sortnums.end(), greater<int>());
  for (int i = 0; i < lineNum; i++) {
    int weight1 = sortnums[i], index1;
    vector<int>::iterator iter1 = find(nums.begin(), nums.end(), weight1);
    index1 = iter1 - nums.begin();
    HoughLine.push_back(Hough_Line[index1]);
  }
  for (int i = 0; i < HoughLine.size(); i++) {
    const int x_min = 0;
    const int x_max = pic.width() - 1;
    const int y_min = 0;
    const int y_max = pic.height() - 1;
    double b = HoughLine[i].b;
    if (HoughLine[i].k < DBL_MAX) {
      double k = HoughLine[i].k;
      const int x0 = (double)(y_min - b) / k;
      const int x1 = (double)(y_max - b) / k;
      const int y0 = x_min * k + b;
      const int y1 = x_max * k + b;
      if (b > 0)
        cout << "Line " << i+1 << ": y=" << k << "*x+" << b << endl;
      else if (b == 0)
        cout << "Line " << i+1 << ": y=" << k << "*x" << endl;
      else
        cout << "Line " << i+1 << ": y=" << k << "*x" << b << endl;
      if (abs(k) > 1) {
        pic.draw_line(x0, y_min, x1, y_max, red);
      }
      else {
        pic.draw_line(x_min, y0, x_max, y1, red);
      }
    }
    else {
      cout << "Line " << i+1 << ": x=" << b << endl;
      pic.draw_line(b , y_min, b, y_max, red);
    }
  }
  //pic.display();
  sort(HoughLine.begin(), HoughLine.end());
  //得出四个角点的坐标
  for (int i = 0; i < lineNum; i++) {
      double k0 = HoughLine[i].k;
      double b0 = HoughLine[i].b;
    for (int j = i+1; j < lineNum; j++) {
      int x = 0, y = 0;
      if (HoughLine[j].k < DBL_MAX) {
        x =  (double)(HoughLine[j].b - b0) / (double)(k0 - HoughLine[j].k);
        y =  (double)(k0*HoughLine[j].b - HoughLine[j].k*b0) / (double)(k0 - HoughLine[j].k);
        
      }
      //K无穷大时的取值
      else {
        x = HoughLine[j].b;
        y = k0*x+b0;
      }
      parameter p;
      p.x = x; p.y = y;
      if (p.x > 0 && p.x < img.width() && p.y > 0 && p.y < img.height()) {
          bool flag1 = 1;
          for (int k = 0; k < HoughPoint.size();k++) {
            if (p.x == HoughPoint[k].x && p.y == HoughPoint[k].y) {
              flag1 = 0;
              break;
            }
          }
          if (flag1) {
            HoughPoint.push_back(p);
          }
        }
    }
  }
  return pic;
}

CImg<unsigned char> Hough::drawPoint() {
  /*for (int i = 0; i < HoughLine.size(); i++) {
    cout << HoughLine[i].k << " " << HoughLine[i].b << endl;
  }*/
  const unsigned char blue[] = { 0, 0, 255 };
  CImg<unsigned char> pic = img;
  for (int i = 0; i < HoughPoint.size(); i++) {
    cout << HoughPoint[i].x << " " << HoughPoint[i].y << endl;
    pic.draw_circle(HoughPoint[i].x, HoughPoint[i].y, 50, blue);
  }
  return pic;
}

void Hough::CircleDetect(int Nums, int min_r, int max_r) {
  thres_img = thres;
  /*thres_img(thres.width(), thres.height(), 1, 3);
  cimg_forXY(thres_img,x,y) {
    thres_img(x,y,0) = thres(x,y);
    thres_img(x,y,1) = thres(x,y);
    thres_img(x,y,2) = thres(x,y);
  }*/
  //thres_img = img;
  thres_img.resize(thres.width(), thres.height(), 1, 3);
  //thres_img.display();
  vector<int> sortCircleWeight;
  vector<pair<int,int>> Circle;
  vector<int> Circleweight;
  vector<pair<int,int>> center;
  int max = 0;
  int width = thres_img.width();
  int height = thres_img.height();
  CImg<int> pic = thres;
  //pic.display();
  vector<pair<int,int>> vote;
  for (int r = min_r; r < max_r; r+=5) {
    HoughImg_Circle.resize(width, height, 1, 1, 0);
    HoughImg_Circle.fill(0);
    max = 0;
    cimg_forXY(pic,x,y) {
      if (pic(x,y) != 0) {
        for (int i = 0; i < 360; i++) {
          int x0 = x - r*cos(i*M_PI/180);
          int y0 = y - r*sin(i*M_PI/180);
          if (x0 > 0 && x0 < width && y0 > 0 && y0 < height)
            HoughImg_Circle(x0,y0)++;
        }
      }
    }
    cimg_forXY(HoughImg_Circle,x,y) {
      if (HoughImg_Circle(x,y) > max) {
        max = HoughImg_Circle(x,y);
      }
    }
    vote.push_back(make_pair(max, r));
  }
  sort(vote.begin(), vote.end(), [](const pair<int, int>& x, const pair<int, int>& y) -> int { return x.first > y.first; });
  //int Nums = 1;
  int knums = 0;
  for (int num = 0; num < Nums; num++) {
    int i = 0;
    HoughImg_Circle.resize(width, height, 1, 1, 0);
    HoughImg_Circle.fill(0);
    int r = vote[num].second;
    cimg_forXY(pic,x,y) {
      if (pic(x,y) != 0) {
        for (int i = 0; i < 360; i++) {
          int x0 = x - r*cos(i*M_PI/180);
          int y0 = y - r*sin(i*M_PI/180);
          if (x0 > 0 && x0 < width && y0 > 0 && y0 < height)
            HoughImg_Circle(x0,y0)++;
        }
      }
    }
    Circle.clear();
    Circleweight.clear();
    cimg_forXY(HoughImg_Circle,x,y) {
      if (HoughImg_Circle(x,y) != 0) {
        Circle.push_back(make_pair(x,y));
        Circleweight.push_back(HoughImg_Circle(x,y));
      }
    }
    unsigned char blue[3] = { 0, 0, 255 };
    sortCircleWeight = Circleweight;
    sort(sortCircleWeight.begin(), sortCircleWeight.end(), greater<int>()); // 将累加矩阵从大到小进行排序
    //同个半径圆的检测
    int count = 0;
    while (1) {
        int weight = sortCircleWeight[count], index;
        vector<int>::iterator iter = find(Circleweight.begin(), Circleweight.end(), weight);
        index = iter - Circleweight.begin();
        int a = Circle[index].first, b = Circle[index].second;
        count++;

        int ii;
        for (ii = 0; ii < center.size(); ii++) {
          // 判断检测出来的圆心坐标是否跟已检测的圆心坐标的距离，如果距离过小，默认是同个圆
            if (sqrt(pow((center[ii].first - a), 2) + pow((center[ii].second - b), 2)) < minRadius) {
                break; 
            }
        }
        if (ii == center.size()) {
            center.push_back(make_pair(a, b));
            thres_img.draw_circle(a, b, r, blue, 1, 1);
            knums++;
            break;
        }
    } 
  }
  cout << "the number of the coins is: " << knums << endl;
  //thres_img.display();
}

void Hough::findpixelCircle() {
  CircleImg = img;
  CircleImg.resize(thres.width(), thres.height());
  cimg_forXY(CircleImg, x, y) {
    if (thres_img(x,y,0) == 0 && thres_img(x,y,1) == 0 && thres_img(x,y,2) == 255 && thres(x,y) == 0) {
      CircleImg(x,y,0) = 255;
      CircleImg(x,y,1) = 0;
      CircleImg(x,y,2) = 0;
    }
  }
  //CircleImg.display();
}