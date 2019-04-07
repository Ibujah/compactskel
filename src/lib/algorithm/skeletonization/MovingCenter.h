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
 *  \file MovingCenter.h
 *  \brief Defines object containing informations about moving centers
 *  \author Bastien Durix
 */

#ifndef _MOVINGCENTER_H_
#define _MOVINGCENTER_H_

#include <Eigen/Dense>
#include <boundary/DiscreteBoundary2.h>
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
			 *  \brief Contains properties about moving centers
			 */
			class MovingCenter
			{
				protected:
					/**
					 *  \brief Position of the center
					 */
					Eigen::Vector2d m_center;

					/**
					 *  \brief Radius of the circle at the center
					 */
					double m_rad;
					
					/**
					 *
					 */
					std::vector<unsigned int> m_closest;

					std::vector<bool> m_next;

					std::list<unsigned int> m_erase;
					
				public:
					/**
					 *  \brief Default constructor
					 */
					MovingCenter();

					/**
					 *  \brief Constructor
					 *
					 *  \param center  position of the center
					 */
					MovingCenter(const Eigen::Vector2d &center);
					
					const Eigen::Vector2d& getCenter() const;

					const double& getRadius() const;

					unsigned int getNbDir() const;
					
					const std::vector<bool>& getNext();

					/**
					 * @brief Computes the contact sets associated to the circle, and the next directions of propagation
					 */
					void computeTangencyData(const std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti, double epsilon);

					/**
					 * @brief Propagates the circle in a given direction
					 */
					bool propagate(const std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti,
								   unsigned int dir,
								   MovingCenter &mov,
								   double epsilon) const;
					
					/**
					 * @brief Gets the indices of the closest points
					 */
					void getIndBnd(std::list<unsigned int>& lsind) const;
					
					/**
					 * @brief Removes points of the boundary that are associated to the given circle: helps to lower the computation time
					 */
					void clean(std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti) const;
					
					/**
					 * @brief Checks if two circles are neighbors
					 */
					static bool neighbors(const MovingCenter &mov1, unsigned int dir1, const MovingCenter &mov2, unsigned int dir2);
			};
		}
	}
}

#endif //_MOVINGCENTER_H_
