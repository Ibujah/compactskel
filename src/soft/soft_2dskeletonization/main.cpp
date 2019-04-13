/*
Copyright (c) 2016 Bastien Durix

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


/**
 *  \brief 2D skeletonization
 *  \author Bastien Durix
 */

#include <boost/program_options.hpp>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <chrono>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>

#include <shape/DiscreteShape.h>
#include <boundary/DiscreteBoundary.h>
#include <skeleton/Skeletons.h>

#include <algorithm/extractboundary/NaiveBoundary.h>
#include <algorithm/skeletonization/SpherePropagation2D.h>
#include <algorithm/skinning/Filling.h>
#include <algorithm/evaluation/ReprojError.h>

#include <displayopencv/DisplayShapeOCV.h>
#include <displayopencv/DisplayBoundaryOCV.h>
#include <displayopencv/DisplaySkeletonOCV.h>

std::tuple<double,double,int,int> EvalSkel(const shape::DiscreteShape<2>::Ptr dissh,
									   const boundary::DiscreteBoundary<2>::Ptr disbnd,
									   const skeleton::GraphSkel2d::Ptr skel)
{
	shape::DiscreteShape<2>::Ptr shp(new shape::DiscreteShape<2>(dissh->getWidth(),dissh->getHeight()));
	algorithm::skinning::Filling(shp,skel);
	
	double res = algorithm::evaluation::SymDiffArea(dissh,shp);
	double res2 = algorithm::evaluation::HausDist(skel,disbnd,dissh->getFrame());
	
	std::list<unsigned int> lnod;
	skel->getAllNodes(lnod);
	unsigned int nbbr = 0;
	for(std::list<unsigned int>::iterator it = lnod.begin(); it != lnod.end(); it++)
	{
		unsigned int deg = skel->getNodeDegree(*it);
		if(deg != 2)
			nbbr += deg;
	}
	nbbr /= 2;
	std::tuple<double,double,int,int> result = std::make_tuple(res*100.0,res2,skel->getNbNodes(),nbbr);
	
	return result;
}

int main(int argc, char** argv)
{
	std::string imgfile, fileimg;
	bool output = false;
	double epsilon;

	boost::program_options::options_description desc("OPTIONS");
	
	desc.add_options()
		("help", "Help message")
		("imgfile", boost::program_options::value<std::string>(&imgfile)->default_value("mask.png"), "Binary image file (*.png)")
		("output", boost::program_options::value<bool>(&output)->implicit_value(true), "Returns output images")
		("epsilon", boost::program_options::value<double>(&epsilon)->default_value(1.0), "Skeleton precision")
		("fileimg", boost::program_options::value<std::string>(&fileimg)->default_value("skeleton.png"), "Skeleton img file")
		;
	
	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
	boost::program_options::notify(vm);
	
	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 0;
	}

	time_t start,end;
	double diff;

	cv::Mat shpimggray = cv::imread(imgfile,cv::ImreadModes::IMREAD_GRAYSCALE);
	cv::Mat shpimg;
	cv::threshold(shpimggray,shpimg,1,255,cv::THRESH_BINARY);
	
	// topological closure of the binary shape
	cv::Mat shpdil;
	cv::Mat element = cv::getStructuringElement(cv::MORPH_RECT,cv::Size(3,3),cv::Point(1,1));
	cv::dilate(shpimg,shpdil,element);
	cv::erode(shpdil,shpimg,element);
	shape::DiscreteShape<2>::Ptr dissh = shape::DiscreteShape<2>::Ptr(new shape::DiscreteShape<2>(shpimg.cols,shpimg.rows));
	cv::Mat cpymat(shpimg.rows,shpimg.cols,CV_8U,&dissh->getContainer()[0]);
	shpimg.copyTo(cpymat);
	cv::Mat image(shpimg.rows,shpimg.cols,CV_8UC3,cv::Scalar(255,255,255));
	
	boundary::DiscreteBoundary<2>::Ptr disbnd = algorithm::extractboundary::NaiveBoundary(dissh);

	if(disbnd->getPrev(0) != disbnd->getNbVertices() -1)
	{
		return -1;
	}
	
	auto start0 = std::chrono::steady_clock::now();
	algorithm::skeletonization::propagation::OptionsSphProp options(2.0*epsilon);
	skeleton::GraphSkel2d::Ptr grskelpropag = algorithm::skeletonization::propagation::SpherePropagation2D(disbnd,options);
	auto duration0 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start0);
	std::tuple<double,double,int,int> respropag = EvalSkel(dissh,disbnd,grskelpropag);
	int t0 = duration0.count();
	double A0 = std::get<0>(respropag); // sym area diff
	double H0 = std::get<1>(respropag); // Hausdorff dist
	int N0 = std::get<2>(respropag); // nb nodes
	int B0 = std::get<3>(respropag); // nb branches

	std::cout << "Skeleton estmated in " << t0 << "ms." << std::endl;
	std::cout << "Symetric difference area over reference area: " << A0 << "\%" << std::endl;
	std::cout << "Hausdorff distance to reference: " << H0 << "px (epsilon=" << epsilon << "px)" << std::endl;
	std::cout << "Number of branches: " << B0 << std::endl;
	std::cout << "Number of nodes: " << N0 << std::endl;
	
	shape::DiscreteShape<2>::Ptr shppropag(new shape::DiscreteShape<2>(dissh->getWidth(),dissh->getHeight()));
	algorithm::skinning::Filling(shppropag,grskelpropag);

	cv::Mat imagepropag;
	image.copyTo(imagepropag);
	displayopencv::DisplayDiscreteShape(dissh,imagepropag,shppropag->getFrame(),cv::Scalar(255,0,0));
	displayopencv::DisplayDiscreteShape(shppropag,imagepropag,shppropag->getFrame(),cv::Scalar(125,125,125));
	displayopencv::DisplayDiscreteBoundary(disbnd,imagepropag,dissh->getFrame(),cv::Scalar(0,0,0));
	displayopencv::DisplayGraphSkeleton(grskelpropag,imagepropag,dissh->getFrame(),cv::Scalar(255,0,0));
	
	if(output)
	{
		cv::imwrite(fileimg, imagepropag);
	}

	return 0;
}
