extern "C" {
#include "mex.h"
#include <stdio.h>
#include <string.h>
}
#include <iostream>
#include <sys/time.h>
#include <opencv2/core/core.hpp>
#include "AprilTags/TagDetector.h"
#include "AprilTags/Tag36h11.h"

// [id, p0, p1, p2, p3, p4] = dapriltags(img);
// this uses the mit apriltag library
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
  if (nrhs < 1)
    mexErrMsgTxt("Not enough input arguments");
  if (!(mxIsUint8(prhs[0]) || mxIsInt8(prhs[0])))
    mexErrMsgTxt("Image must be of type uint8 or int8");

  // Get Image
  char *image_buf = (char *) mxGetData(prhs[0]);
  int height = mxGetM(prhs[0]);
  int width  = mxGetN(prhs[0]);
  cv::Mat image = cv::Mat::zeros(cv::Size(height, width), CV_8U);
  memcpy(image.data, image_buf, width * height);
  cv::transpose(image, image);

  // Initialize a tag detector and detect tags
  static AprilTags::TagDetector tag_detector(AprilTags::tagCodes36h11);
  std::vector<AprilTags::TagDetection> detections = tag_detector.extractTags(image);

  /*
    // Remove tags that is not perfectly matched
    vector<AprilTags::TagDetection> _detections;
    for (unsigned int k = 0; k < detections.size(); k++)
      if (detections[k].hammingDistance == 0)
        _detections.push_back(detections[k]);
    detections = _detections;
  */

  // Output detection, use tag id instead of tag code
  int N = detections.size();
  //plhs[0] = mxCreateNumericMatrix(1, N, mxINT64_CLASS, mxREAL);
  plhs[0] = mxCreateNumericMatrix(1, N, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(2, N, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(2, N, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(2, N, mxREAL);
  plhs[4] = mxCreateDoubleMatrix(2, N, mxREAL);
  plhs[5] = mxCreateDoubleMatrix(2, N, mxREAL);
  // Tag id
  //long long* ptri = (long long*)mxGetData(plhs[0]);
  double *ptri = (double *)mxGetData(plhs[0]);
  for (int i = 0; i < N; i++)
    ptri[i] = detections[i].id;

  for (int k = 1; k < 6; k++) {
    double *ptr = mxGetPr(plhs[k]);
    for (int i = 0; i < N; i++) {
      if (k == 1) {  // Tag position
        ptr[2 * i]   = detections[i].cxy.first;
        ptr[2 * i + 1] = detections[i].cxy.second;
      } else {       // Tag vertices
        ptr[2 * i]   = detections[i].p[k - 2].first;
        ptr[2 * i + 1] = detections[i].p[k - 2].second;
      }
    }
  }
  return;
}



