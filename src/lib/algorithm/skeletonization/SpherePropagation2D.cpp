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
 *  @file SpherePropagation2D.cpp
 *  @brief Defines functions to compute 2d skeleton with sphere propagation algorithm
 *  @author Bastien Durix
 */

#include "SpherePropagation2D.h"
#include "BoundaryOperations.h"
#include "MovingCenter.h"
#include <time.h>

#include <iostream>
#include <fstream>

#include <unordered_map>

skeleton::GraphSkel2d::Ptr algorithm::skeletonization::propagation::SpherePropagation2D(const boundary::DiscreteBoundary<2>::Ptr disbnd,
																						std::map<unsigned int, std::list<unsigned int> > &pt_assoc_skel,
																						const OptionsSphProp &options)
{
	OptiBnd optiBnd;
	skeleton::GraphSkel2d::Ptr skel(new skeleton::GraphSkel2d(skeleton::model::Classic<2>()));
	
	// creation of the optimized boundary structure
	createOptiBnd(disbnd,optiBnd);
	
	// estimation of a point inside the shape
	Eigen::Vector2d firstPt = firstPoint(optiBnd);
	
	// estimation of the first Voronoi point
	std::list<unsigned int> lclosest;
	algorithm::skeletonization::propagation::firstVoroPoint(optiBnd,firstPt,lclosest);
	
	// first moving center instance
	algorithm::skeletonization::propagation::MovingCenter mov(firstPt);
	mov.computeContactData(optiBnd,options.epsilon);
	
	// search for a local maximum in the Voronoi diagram to be sure that we have a skeletal point
	bool maximised = true;
	do
	{
		algorithm::skeletonization::propagation::MovingCenter movn;
		maximised = false;
		for(unsigned int i = 0; i < mov.getNext().size(); i++)
		{
			if(mov.getNext()[i])
			{
				algorithm::skeletonization::propagation::MovingCenter movncur;
				mov.propagate(optiBnd,i,movncur,options.epsilon);
				
				if(maximised)
				{
					if(movncur.getRadius() > movn.getRadius())
					{
						movn = movncur;
					}
				}
				else
				{
					if(movncur.getRadius() > mov.getRadius())
					{
						movn = movncur;
						maximised = true;
					}
				}
			}
		}
		if(maximised)
			mov = movn;
	}while(maximised);
	
	mov = algorithm::skeletonization::propagation::MovingCenter(mov.getCenter());
	mov.computeTangencyData(optiBnd,options.epsilon);

	double rad = mov.getRadius();
	unsigned int ind = skel->addNode(Eigen::Vector3d(mov.getCenter().x(),mov.getCenter().y(),rad));
	std::list<unsigned int> bndneigh;
	mov.getIndBnd(bndneigh);
	pt_assoc_skel.insert(std::make_pair(ind,bndneigh));
	mov.clean(optiBnd);
	
	std::list<std::tuple<unsigned int,algorithm::skeletonization::propagation::MovingCenter,unsigned int> > lctr;
	for(unsigned int i = 0; i < mov.getNext().size(); i++)
	{
		if(mov.getNext()[i])
		{
			std::tuple<unsigned int,algorithm::skeletonization::propagation::MovingCenter,unsigned int> cdir;
			cdir = std::make_tuple(ind,mov,i);
			lctr.push_back(cdir);
		}
	}
	
	unsigned int cpt = options.iter_max;
	
	if(lctr.size() != 0)
	{
		do
		{
			std::tuple<unsigned int,algorithm::skeletonization::propagation::MovingCenter,unsigned int> cdir = *(lctr.begin());
			lctr.pop_front();
			
			algorithm::skeletonization::propagation::MovingCenter movn;
			if(std::get<1>(cdir).propagate(optiBnd,std::get<2>(cdir),movn,options.epsilon))
			{
				std::list<unsigned int> bndneighcur;
				movn.getIndBnd(bndneighcur);
				movn.clean(optiBnd);

				unsigned int indcur = skel->addNode(Eigen::Vector3d(movn.getCenter().x(),movn.getCenter().y(),movn.getRadius()));
				pt_assoc_skel.insert(std::make_pair(indcur,bndneighcur));
				skel->addEdge(indcur,std::get<0>(cdir));
				
				for(unsigned int i = 0; i < movn.getNext().size(); i++)
				{
					if(movn.getNext()[i])
					{
						bool skip = false;
						for(std::list<std::tuple<unsigned int,algorithm::skeletonization::propagation::MovingCenter,unsigned int> >::iterator it = lctr.begin(); it != lctr.end() && !skip; it++)
						{
							if(algorithm::skeletonization::propagation::MovingCenter::neighbors(movn,i,std::get<1>(*it),std::get<2>(*it)))
							{
								skel->addEdge(indcur,std::get<0>(*it));
								it = lctr.erase(it);
								skip = true;
							}
						}
						
						if(!skip)
						{
							std::tuple<unsigned int,algorithm::skeletonization::propagation::MovingCenter,unsigned int> cdirn;
							cdirn = std::make_tuple(indcur,movn,i);
							lctr.push_back(cdirn);
						}
					}
				}
			}
			else
			{
				cpt = 1;
			}
			
			if(options.epsilon != 0.0)
				cpt--;
		}while(!lctr.empty() && cpt != 0);
	}
	if(cpt == 0)
		throw std::logic_error("arf");
	
	return skel;
}
