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
 *  \file DistanceField.h
 *  \brief Defines functions related to distance field
 *  \author Bastien Durix
 */

#ifndef _DISTANCEFIELD_H_
#define _DISTANCEFIELD_H_

#include <boundary/DiscreteBoundary2.h>
#include <shape/DiscreteShape.h>

#include <unordered_map>

/**
 *  \brief Lots of algorithms
 */
namespace algorithm
{
	/**
	 *  \brief skeletonization algorithms
	 */
	namespace skeletonization
	{
		namespace propagation
		{
			/**
			 * @brief Finds the closest points on the boundary to a given center
			 */
			std::list<unsigned int> closestInd(const std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti, const Eigen::Vector2d &C, double &radmin);
			
			/**
			 * @brief Reorders a set of points in direct order around a given center
			 */
			void reorder(const std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti,
						 const std::list<unsigned int> &linds,
						 const Eigen::Vector2d &C,
						 std::vector<unsigned int> &vinds);
			
			/**
			 * @brief Finds a first point inside a given shape, to initialize the algorithm
			 */
			Eigen::Vector2d firstPoint(const boundary::DiscreteBoundary<2>::Ptr disbnd);
			
			/**
			 * @brief Computes the contact sets between each pair of closest points, and gives the next propagation directions
			 */
			void tangencyBoundary(const std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti,
								  const Eigen::Vector2d &C,
								  const double &radmax,
								  const std::vector<unsigned int> &vclo,
								  std::vector<bool> &next,
								  std::list<unsigned int> &lerase);
			
			/**
			 * @brief (re)Computes a center from a set of boundary points
			 */
			void correction(const std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti, Eigen::Vector2d &C, std::list<unsigned int> &lclosest);
			
			/**
			 * @brief (re)Computes a center from a set of 3 boundary points
			 */
			Eigen::Vector2d circleCenter(const Eigen::Vector2d &P1, const Eigen::Vector2d &P2, const Eigen::Vector2d &P3);
		}
	}
}

#endif //_DISTANCEFIELD_H_
