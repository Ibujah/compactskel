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
 *  \file SpherePropagation2D.h
 *  \brief Defines functions to compute 2d skeleton with sphere propagation algorithm
 *  \author Bastien Durix
 */

#ifndef _SPHEREPROPAGATION_H_
#define _SPHEREPROPAGATION_H_

#include <skeleton/Skeletons.h>
#include <boundary/DiscreteBoundary2.h>
#include <shape/DiscreteShape.h>

#define _NB_CIRCLES_MAX_ 10000

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
			 *  \brief Sphere propagation options structure
			 */
			struct OptionsSphProp
			{
				/**
				 *  \brief Hausdorff distance
				 */
				double epsilon;
				
				unsigned int iter_max;
				/**
				 *  \brief Default constructor
				 */
				OptionsSphProp(double epsilon_ = 2.1, unsigned int iter_max_=_NB_CIRCLES_MAX_) : epsilon(epsilon_), iter_max(iter_max_)
				{}
			};
			
			/**
			 *  \brief 2D skeletonization, by sphere propagation
			 *
			 *  \param disbnd   discrete boundary of the shape
			 *  \param disshp   discrete shape
			 *  \param options  options of the algorithm
			 *
			 *  \return pointer to the computed 2d graph skeleton
			 */
			skeleton::GraphSkel2d::Ptr SpherePropagation2D(const boundary::DiscreteBoundary<2>::Ptr disbnd, std::map<unsigned int, std::list<unsigned int> > &pt_assoc_skel, const OptionsSphProp &options = OptionsSphProp());
		}
	}
}

#endif //_SPHEREPROPAGATION_H_
