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
 *  \file SpherePropagation2D.cpp
 *  \brief Defines functions to compute 2d skeleton with sphere propagation algorithm
 *  \author Bastien Durix
 */

#include "SpherePropagation2D.h"
#include "DistanceField.h"
#include "MovingCenter.h"
#include <time.h>

#include <thread>
#include <mutex>

#include <iostream>
#include <fstream>

#include <unordered_map>


skeleton::GraphSkel2d::Ptr SpherePropagation2D_helper(const boundary::DiscreteBoundary<2>::Ptr disbnd,
													  std::map<unsigned int, std::list<unsigned int> > &pt_assoc_skel,
													  Eigen::Vector2d &fstpt,
													  const algorithm::skeletonization::propagation::OptionsSphProp &options)
{
	std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > bndopti;
	skeleton::GraphSkel2d::Ptr skel(new skeleton::GraphSkel2d(skeleton::model::Classic<2>()));
	
	for(unsigned int i = 0; i < disbnd->getNbVertices(); i++)
	{
		std::tuple<Eigen::Vector2d,unsigned int,unsigned int> pt = 
			std::make_tuple(disbnd->getCoordinates(i),disbnd->getPrev(i),disbnd->getNext(i)); // be sure it's in trigonometric order
		bndopti.insert(std::make_pair(i,pt));
	}
	
	std::list<unsigned int> lclosest;
	algorithm::skeletonization::propagation::correction(bndopti,fstpt,lclosest);

	algorithm::skeletonization::propagation::MovingCenter mov(fstpt);
	mov.computeTangencyData(bndopti,options.epsilon);
	
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
				mov.propagate(bndopti,i,movncur,options.epsilon);
				
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
	mov.computeTangencyData(bndopti,options.epsilon);

	double rad = mov.getRadius();
	unsigned int ind = skel->addNode(Eigen::Vector3d(mov.getCenter().x(),mov.getCenter().y(),rad));
	std::list<unsigned int> bndneigh;
	mov.getIndBnd(bndneigh);
	pt_assoc_skel.insert(std::make_pair(ind,bndneigh));
	mov.clean(bndopti);
	
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
			if(std::get<1>(cdir).propagate(bndopti,std::get<2>(cdir),movn,options.epsilon))
			{
				std::list<unsigned int> bndneighcur;
				movn.getIndBnd(bndneighcur);
				movn.clean(bndopti);

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

skeleton::GraphSkel2d::Ptr algorithm::skeletonization::propagation::SpherePropagation2D(const boundary::DiscreteBoundary<2>::Ptr disbnd,
																						std::map<unsigned int, std::list<unsigned int> > &pt_assoc_skel,
																						const OptionsSphProp &options)
{
	srand(time(NULL));
	Eigen::Vector2d C = firstPoint(disbnd);
	
	return SpherePropagation2D_helper(disbnd,pt_assoc_skel,C,options);
}
