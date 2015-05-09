//============================================================================
// Name        : PedestrianDetection.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

// Setting up: C++11
// http://stackoverflow.com/questions/9131763/eclipse-cdt-c11-c0x-support
// Setting up opencv:
// http://docs.opencv.org/doc/tutorials/introduction/linux_eclipse/linux_eclipse.html
// Matio is required to read and write mat-files.
// http://sourceforge.const double PI  =3.141592653589793238463;net/projects/matio/

#include <iostream>
#include <cv.h>
#include <highgui.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include "PedestrianDetector.hpp"
#include <dirent.h>
#include <string.h>
#include <math.h>
#include <iomanip>

using namespace std;
using namespace cv;
void showBB(vector<bbType> bbs, Mat image, bool wait){
	double alpha = 0.3;
	double threshold = 50;
	for(int iBbs = 0; iBbs < bbs.size(); iBbs++) {
		//printf("n%d, x: %f, %f, %f, %f, %f, Dist: %f \n",iBbs,bbs[iBbs].x1,bbs[iBbs].y2, bbs[iBbs].width3, bbs[iBbs].height4, bbs[iBbs].score5, bbs[iBbs].distance );

		Scalar useColor(0,0,0);
		if(bbs[iBbs].score5 > threshold) {
			alpha = 0.3;
			useColor[1] = 0; // G
			useColor[2] = 255; // R
		}
		else {
			alpha = 0.1;
			useColor[1] = 255; // G
			useColor[2] = 0; // R
		}
		Mat rectangleImage(image.size[0],image.size[1], CV_8UC3, cv::Scalar(0, 0,0));
		rectangle(rectangleImage, Rect(bbs[iBbs].x1,bbs[iBbs].y2,bbs[iBbs].width3, bbs[iBbs].height4), useColor, CV_FILLED, 8, 0 );

		stringstream strsScore; strsScore.precision(3);
		strsScore << bbs[iBbs].score5;
		putText(image, strsScore.str() + "p", Point(bbs[iBbs].x1,bbs[iBbs].y2), FONT_HERSHEY_SIMPLEX, 0.75, Scalar(255,255,255)); //, int thickness=1, int lineType=8, bool bottomLeftOrigin=false

		stringstream strsDistance; strsDistance.precision(3);
		strsDistance << bbs[iBbs].distance;
		putText(image, strsDistance.str() + "m", Point(bbs[iBbs].x1,bbs[iBbs].y2+bbs[iBbs].height4), FONT_HERSHEY_SIMPLEX, 0.75, Scalar(255,255,255)); //, int thickness=1, int lineType=8, bool bottomLeftOrigin=false

		stringstream strsAngle; strsAngle.precision(3);
		strsAngle << bbs[iBbs].angle;
		putText(image, "Angle" + strsAngle.str(), Point(bbs[iBbs].x1,bbs[iBbs].y2+bbs[iBbs].height4 +25), FONT_HERSHEY_SIMPLEX, 0.75, Scalar(255,255,255)); //, int thickness=1, int lineType=8, bool bottomLeftOrigin=false

		//putText(image, to_string((int)(round(bbs[iBbs].score5))), Point(bbs[iBbs].x1,bbs[iBbs].y2), FONT_HERSHEY_SIMPLEX, 1, Scalar(255,255,255)); //, int thickness=1, int lineType=8, bool bottomLeftOrigin=false
		//putText(image, "d:" + to_string((int)(round(bbs[iBbs].distance))), Point(bbs[iBbs].x1,bbs[iBbs].y2+bbs[iBbs].height4), FONT_HERSHEY_SIMPLEX, 1, Scalar(255,255,255)); //, int thickness=1, int lineType=8, bool bottomLeftOrigin=false
		addWeighted(rectangleImage, alpha, image, 1, 0.0, image );
	}

	if(wait)
		namedWindow("BoundingsBox", CV_WINDOW_AUTOSIZE);
	imshow("BoundingsBox",image);

	//printf("Code is done!!");



}
int main(int argc, char** argv) {


	string dirImage = "pedmodels/InriaTest.png";
	//string dirImage = "pedmodels/KimTest.jpg";
	//string dirImage = "pedmodels/MarkTest.jpg";
	vector<bbType> bbs;
	double FOV_verticalDeg = 47; // Vertical field-of-view of camera.
	double FOV_horizontalDeg = 50; //83;
	double angleTiltDegrees = 7; // Downward tilt in degrees.
	double cameraHeight = 1.9; // Height Position of camera.

	double imageResize = 0.5;
	//Mat image = imread(dirImage, 1);

	Mat image, inputImage;
	struct dirent *ent;
	//string dirImages = "/home/pistol/Desktop/DataFolder/2014-11-03-14-37-11/WebCam";
	//string dirImages = "/home/pistol/Desktop/DataFolder/2014-10-16-12-18-30/WebCam";
	string dirImages = "/home/pistol/Desktop/DataFolder/WebCam";

	DIR *dir = opendir(dirImages.data());

	char * pch;


	//image = imread(dirImage, 1);

	//cout << image.cols << "x" << image.rows << endl;

	// Making detector object.
	bool fastDetector = 0;
	string dirDetector;
	if(fastDetector) {
		dirDetector = "pedmodels/AcfInriaDetector.mat";
	}
	else {
		dirDetector = "pedmodels/LdcfInriaDetector.mat";
	}

	PedestrianDetector oPedDetector(dirDetector);
	// Providing camera settings.
	oPedDetector.setCameraSetup(FOV_verticalDeg, FOV_horizontalDeg, angleTiltDegrees, cameraHeight);


	bool useExample = true;
	if(useExample) {
		inputImage = imread(dirImage, 1);
		if (!inputImage.data) {
			printf("No image data \n");
			return -1;
		}
		resize(inputImage, image, Size(), imageResize, imageResize);

		clock_t start, end;
		start = clock();
		bbs = oPedDetector.pedDetector(image);
		end = clock();
		double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
		cout << "\n t1:" << time << " ms\n";

		showBB(bbs, image,0);
		waitKey(0);

	}

	bool useWholeFolder = true;
	if(useWholeFolder) {
		struct dirent **namelist;
		int i;
		string fileName;
		//string pch;
		int n = scandir(dirImages.data(), &namelist, 0, alphasort);
		if (n < 0)
			perror("scandir");
		else {
			for (i = 0; i < n; i++) {
				printf("%s\n", namelist[i]->d_name);


				//while ((ent = readdir (dir)) != NULL) {
				fileName = namelist[i]->d_name;
				string pch=strrchr(namelist[i]->d_name,'.');
				free(namelist[i]);
				if(pch.compare(".jpg") == 0)
				{
					//printf ("%s\n", fileName.data());
					string dirImage = (dirImages+ '/' +fileName);
					inputImage = imread(dirImage.data(), 1);
					resize(inputImage, image, Size(), imageResize,imageResize);

					clock_t start, end;
					start = clock();
					bbs = oPedDetector.pedDetector(image);
					end = clock();
					double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
					cout << "\n t1:" << time << " ms\n";
					showBB(bbs, image,0);waitKey(10);

				}
				//closedir (dir);
			}
		}
		free(namelist);
	}
	return 0;
}

