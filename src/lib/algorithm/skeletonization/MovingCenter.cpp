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
 *  \file MovingCenter.cpp
 *  \brief Defines object containing informations about moving centers
 *  \author Bastien Durix
 */

#include "MovingCenter.h"
#include "DistanceField.h"
#include <set>
#include <iostream>
#include <bitset>

algorithm::skeletonization::propagation::MovingCenter::MovingCenter() {}

algorithm::skeletonization::propagation::MovingCenter::MovingCenter(const Eigen::Vector2d &center) : m_center(center) {}

const Eigen::Vector2d& algorithm::skeletonization::propagation::MovingCenter::getCenter() const
{
	return m_center;
}

const double& algorithm::skeletonization::propagation::MovingCenter::getRadius() const
{
	return m_rad;
}

unsigned int algorithm::skeletonization::propagation::MovingCenter::getNbDir() const
{
	unsigned int nbdir = 0;
	for(unsigned int i = 0; i < m_next.size(); i++)
		if(m_next[i])
			nbdir++;
	return nbdir;
}

const std::vector<bool>& algorithm::skeletonization::propagation::MovingCenter::getNext()
{
	return m_next;
}

void algorithm::skeletonization::propagation::MovingCenter::computeTangencyData(const std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti, double epsilon)
{
	double radmin;
	std::list<unsigned int> lclosest = closestInd(bndopti,m_center,radmin);
	
	reorder(bndopti,lclosest,m_center,m_closest);
	m_next.resize(m_closest.size());
	tangencyBoundary(bndopti,m_center,radmin+epsilon,m_closest,m_next,m_erase);

	if(epsilon != 0.0)
	{
		double dmax = radmin;
		for(std::list<unsigned int>::iterator it = m_erase.begin(); it != m_erase.end(); it++)

		{
			Eigen::Vector2d pt = std::get<0>(bndopti.find(*it)->second);
			double dist = (m_center - pt).norm();
			if(dist > dmax)
				dmax = dist;
		}
		m_rad = (radmin + dmax)/2.0;
	}
	else
	{
		m_rad = radmin;
	}
}

bool algorithm::skeletonization::propagation::MovingCenter::propagate(const std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti,
																		  unsigned int dir,
																		  MovingCenter &mov,
																		  double epsilon) const
{
	unsigned int dir2 = (dir+1)%m_closest.size();
	unsigned int indP1 = m_closest[dir];
	unsigned int indP2 = m_closest[dir2];
	Eigen::Vector2d P1 = std::get<0>(bndopti.find(indP1)->second);
	Eigen::Vector2d P2 = std::get<0>(bndopti.find(indP2)->second);
	Eigen::Vector2d mid = (P1 + P2)*0.5;
	double sqmid = (P1 - mid).squaredNorm();
	Eigen::Vector2d nor(P2.y() - P1.y(),P1.x() - P2.x());
	nor.normalize();

	unsigned int indP3;
	double lambdamin = 0.0;
	bool found1 = false;
	for(auto it = bndopti.begin(); it != bndopti.end(); it++)
	{
		if(std::find(m_closest.begin(),m_closest.end(),it->first) == m_closest.end())
		{
			Eigen::Vector2d pt = std::get<0>(it->second);
			if((mid - pt).dot(nor) != 0.0)
			{
				double lambda = (sqmid - (pt - mid).squaredNorm() - 2.0*(mid - pt).dot(m_center - mid))/(2.0*(mid - pt).dot(nor));
				
				if(lambda > 0.0)
				{
					if(!found1 || lambda < lambdamin)
					{
						lambdamin = lambda;
						indP3 = it->first;
						found1 = true;
					}
				}
			}
		}
	}
	
	if(!found1)
	{
		return false;
	}
	
	Eigen::Vector2d P3 = std::get<0>(bndopti.find(indP3)->second);
	Eigen::Vector2d ctr = circleCenter(P1,P2,P3);
	bool converged = true;
	
	if(converged)
	{
		mov = MovingCenter(ctr);
		mov.computeTangencyData(bndopti,epsilon);

		// direction selection
		unsigned int dirmov = 0;
		bool found = false;
		for(unsigned int i = 0; i < mov.m_closest.size() && !found; i++)
		{
			if(neighbors(mov,i,*this,dir))
			//if(mov.m_closest[i] == m_closest[dir2])
			{
				dirmov = i;
				found = true;
			}
		}
		if(found)
		{
			mov.m_next[dirmov] = false;
		}
	}
	
	return converged;
}

void algorithm::skeletonization::propagation::MovingCenter::getIndBnd(std::list<unsigned int>& lsind) const
{
	lsind.insert(lsind.end(),m_closest.begin(),m_closest.end());
	lsind.insert(lsind.end(),m_erase.begin(),m_erase.end());
}

void algorithm::skeletonization::propagation::MovingCenter::clean(std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti) const
{
	for(std::list<unsigned int>::const_iterator it = m_erase.begin(); it != m_erase.end(); it++)
	{
		bndopti.erase(*it);
	}
}	

bool algorithm::skeletonization::propagation::MovingCenter::neighbors(const MovingCenter &mov1, unsigned int dir1, const MovingCenter &mov2, unsigned int dir2)
{
	unsigned int dir1_2 = (dir1 + 1)%mov1.m_closest.size();
	unsigned int dir2_2 = (dir2 + 1)%mov2.m_closest.size();
	bool eq1 = (mov1.m_closest[dir1] == mov2.m_closest[dir2_2]);
	bool eq2 = (mov2.m_closest[dir2] == mov1.m_closest[dir1_2]);
	
	return (eq1 && eq2);
}
