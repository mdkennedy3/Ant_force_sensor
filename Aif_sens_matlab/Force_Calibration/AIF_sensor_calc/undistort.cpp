#include "mex.h"
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "opencv2/opencv.hpp"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{   
  if (nrhs < 3)
    mexErrMsgTxt("Not enough input arguments");
  if (!(mxIsUint8(prhs[0]) || mxIsInt8(prhs[0])))
    mexErrMsgTxt("Image must be of type uint8 or int8");
  if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != 3 || !mxIsDouble(prhs[1]))
    mexErrMsgTxt("Invalid camera matrix (3x3 double)");      
  if (mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 5 || !mxIsDouble(prhs[2]))
    mexErrMsgTxt("Invalid distortion vector (1x5 double)");      
    
  // Get Image
  char* image_buf = (char*)mxGetData(prhs[0]);    
  int height = mxGetM(prhs[0]);
  int width  = mxGetN(prhs[0]);
  cv::Mat image = cv::Mat::zeros(cv::Size(height, width), CV_8U);  
  memcpy(image.data, image_buf, width*height);
  cv::transpose(image, image);
  
  // Calibration parameters
  double* Aptr = mxGetPr(prhs[1]);
  cv::Mat A = cv::Mat::eye(cv::Size(3,3), CV_32F);
  for (unsigned int j = 0; j < 3; j++)
    for (unsigned int i = 0; i < 3; i++)
      A.at<float>(i,j) = Aptr[3*j+i];
  double* Dptr = mxGetPr(prhs[2]);  
  cv::Mat D = cv::Mat::zeros(cv::Size(1,5), CV_32F);
  for (unsigned int k = 0; k < 5; k++)
    D.at<float>(0,k) = Dptr[k];
  
  // Undistort
  cv::Mat image_undistort;
  cv::undistort(image, image_undistort, A, D); 
  
  // Output
  plhs[0] = mxCreateNumericMatrix(height, width, mxUINT8_CLASS, mxREAL);
  unsigned char* ptr = (unsigned char*)mxGetData(plhs[0]);  
  for (unsigned int c = 0; c < width; c++)
    for (unsigned int r = 0; r < height; r++)
      ptr[c*height+r] = image_undistort.at<unsigned char>(r, c);
  return;  
}
