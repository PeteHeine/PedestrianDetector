/*
 * rgbConvertMex.hpp
 *
 *  Created on: Feb 6, 2015
 *      Author: pistol
 */

#ifndef RGBCONVERTMEX_HPP_
#define RGBCONVERTMEX_HPP_
#include <cv.h>

cv::Mat pdRgbConvert(cv::Mat image);
template<class iT, class oT> oT* rgbConvert( iT *I, int n, int d, int flag, oT nrm );

#endif /* RGBCONVERTMEX_HPP_ */
