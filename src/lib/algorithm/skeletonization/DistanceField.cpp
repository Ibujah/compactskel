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
 *  \file DistanceField.cpp
 *  \brief Defines functions related to distance field
 *  \author Bastien Durix
 */

#include "DistanceField.h"
#include <boost/math/distributions/normal.hpp>
#include <set>
#include <nlopt.hpp>
#include <iostream>

std::list<unsigned int> algorithm::skeletonization::propagation::closestInd(const std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti,
																				const Eigen::Vector2d &C,
																				double &radmin)
{
	double distmin = -1.0;
	for(auto it = bndopti.begin(); it != bndopti.end(); it++)
	{
		double dist = (std::get<0>(it->second) - C).norm();
		
		if(dist < distmin || distmin == -1.0)
		{
			distmin = dist;
		}
	}
	
	std::list<unsigned int> inds;
	for(auto it = bndopti.begin(); it != bndopti.end(); it++)
	{
		float dist = (std::get<0>(it->second) - C).norm();

		if(dist == (float)distmin)// || std::abs(dist-distmin) <= 10.0*std::abs(std::min(dist,distmin))*std::numeric_limits<float>::epsilon())
		{
			inds.push_back(it->first);
		}
	}
	
	radmin = distmin;

	return inds;
}

Eigen::Vector2d algorithm::skeletonization::propagation::firstPoint(const boundary::DiscreteBoundary<2>::Ptr disbnd)
{
	Eigen::Vector2d corner = disbnd->getCoordinates(0);
	
	for(unsigned int i = 0; i < disbnd->getNbVertices(); i++)
	{
		Eigen::Vector2d pt = disbnd->getCoordinates(i);
		if(corner.x() > pt.x())
			corner.x() = pt.x();
		if(corner.y() > pt.y())
			corner.y() = pt.y();
	}
	corner.x() -= 1;
	corner.y() -= 1;
	
	double distmin = -1.0;
	unsigned int ind;
	for(unsigned int i = 0; i < disbnd->getNbVertices(); i++)
	{
		Eigen::Vector2d pt = disbnd->getCoordinates(i);
		double dist = (corner - pt).norm();
		if(dist < distmin || distmin == -1.0)
		{
			distmin = dist;
			ind = i;
		}
	}
	
	Eigen::Vector2d ptmin = disbnd->getCoordinates(ind);
	double distintermin = -1.0;
	for(unsigned int i = 0; i < disbnd->getNbVertices(); i++)
	{
		if(i != ind && ptmin != disbnd->getCoordinates(i))
		{
			Eigen::Vector2d pt = disbnd->getCoordinates(i);
			double distinter = (ptmin - pt).norm();
			if(distinter < distintermin || distintermin == -1.0)
			{
				distintermin = distinter;
			}
		}
	}
	
	Eigen::Vector2d C = ptmin + 0.1*distintermin*((ptmin - corner).normalized());
	
	return C;
}

void algorithm::skeletonization::propagation::tangencyBoundary(const std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti,
																   const Eigen::Vector2d &C,
																   const double &radmax,
																   const std::vector<unsigned int> &vclo,
																   std::vector<bool> &next,
																   std::list<unsigned int> &lerase)
{
	for(unsigned int i = 0; i < vclo.size(); i++)
	{
		unsigned int indend = vclo[(i+1)%vclo.size()];
		unsigned int indbeg = vclo[i];
		auto itcur = bndopti.find(indbeg);
		itcur = bndopti.find(std::get<2>(itcur->second));
		bool fini = false;
		bool closed = false;
		std::list<unsigned int> toerase;
		while(!fini)
		{
			if(itcur == bndopti.end())
			{
				fini = true;
			}
			else if(itcur->first == indend)
			{
				closed = true;
				fini = true;
			}
			else if(itcur->first == indbeg)
			{
				closed = false;
				fini = true;
			}
			else
			{
				double dist = (C - std::get<0>(itcur->second)).norm();
				if(dist > radmax)
				{
					fini = true;
				}
				else
				{
					toerase.push_back(itcur->first);
					itcur = bndopti.find(std::get<2>(itcur->second));
				}
			}
		}
		if(closed)
		{
			for(std::list<unsigned int>::iterator it = toerase.begin(); it != toerase.end(); it++)
			{
				lerase.push_back(*it);
			}
			next[i] = false;
		}
		else
		{
			next[i] = true;
		}
	}
}

void algorithm::skeletonization::propagation::correction(const std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti,
															 Eigen::Vector2d &C,
															 std::list<unsigned int> &lclosest)
{
	double radmin;
	lclosest = closestInd(bndopti,C,radmin);

	if(lclosest.size() == 1)
	{
		unsigned int ind1 = *(lclosest.begin());
		Eigen::Vector2d P1 = std::get<0>(bndopti.find(ind1)->second);
		
		Eigen::Vector2d nor = (C - P1).normalized();
		
		double radmin = -1.0;
		double ind2;
		for(auto it = bndopti.begin(); it != bndopti.end(); it++)
		{
			Eigen::Vector2d P2 = std::get<0>(it->second);
			if((P2 - P1).dot(nor) > 0.0 && ind1 != it->first)
			{
				double rad = 0.5*(P2 - P1).norm()/(nor.dot((P2 - P1).normalized()));
				if(rad < radmin || radmin == -1.0)
				{
					radmin = rad;
					ind2 = it->first;
					C = P1 + rad*nor;
				}
			}
		}
		lclosest = closestInd(bndopti,C,radmin);
	}
	
	if(lclosest.size() < 2)
		throw std::logic_error("Not enough points (a)");
	
	if(lclosest.size() == 2)
	{
		std::list<unsigned int>::iterator it = lclosest.begin();
		unsigned int ind1 = *it;
		it++;
		unsigned int ind2 = *it;
		
		Eigen::Vector2d P1 = std::get<0>(bndopti.find(ind1)->second);
		Eigen::Vector2d P2 = std::get<0>(bndopti.find(ind2)->second);

		Eigen::Vector2d nor(P2.y() - P1.y(),P1.x() - P2.x());
		nor.normalize();

		double radmin1 = -1.0;
		double ind3_1;
		double radmin2 = -1.0;
		double ind3_2;
		Eigen::Vector2d C1, C2;

		for(auto it = bndopti.begin(); it != bndopti.end(); it++)
		{
			Eigen::Vector2d P3 = std::get<0>(it->second);
			double scal = (P3 - P1).dot(nor); 
			if(scal > 0.0 && it->first != ind1 && it->first != ind2)
			{
				Eigen::Vector2d ctr = circleCenter(P1,P2,P3);
				double rad = (ctr - P3).norm();
				if(rad < radmin1 || radmin1 == -1.0)
				{
					radmin1 = rad;
					ind3_1 = it->first;
					C1 = ctr;
				}
			}
			if(scal < 0.0 && it->first != ind1 && it->first != ind2)
			{
				Eigen::Vector2d ctr = circleCenter(P1,P2,P3);
				double rad = (ctr - P3).norm();
				if(rad < radmin2 || radmin2 == -1.0)
				{
					radmin2 = rad;
					ind3_2 = it->first;
					C2 = ctr;
				}
			}
		}
		
		if(std::get<1>(bndopti.find(ind1)->second) == ind2 || std::get<2>(bndopti.find(ind1)->second) == ind2) // are neighbors
		{
			if(nor.dot(P1 - C)*nor.dot(P1 - C1) > nor.dot(P1 - C)*nor.dot(P1 - C2))
			{
				C = C1;
			}
			else
			{
				C = C2;
			}
		}
		else
		{
			if(radmin1 != -1.0 && radmin2 != -1.0)
			{
				if(radmin1 > radmin2)
					C = C1;
				else
					C = C2;
			}
			else if(radmin1 != -1.0)
			{
				C = C1;
			}
			else if(radmin2 != -1.0)
			{
				C = C2;
			}
			else
			{
				throw std::logic_error("arf");
			}
		}
		lclosest = closestInd(bndopti,C,radmin);
	}
	
	if(lclosest.size() < 3)
		throw std::logic_error("Not enough points (b)");
}

Eigen::Vector2d algorithm::skeletonization::propagation::circleCenter(const Eigen::Vector2d &P1, const Eigen::Vector2d &P2, const Eigen::Vector2d &P3)
{
	Eigen::Matrix<double,3,2> mat;
	mat.block<1,2>(0,0) = (P2 - P1).transpose();
	mat.block<1,2>(1,0) = (P3 - P2).transpose();
	mat.block<1,2>(2,0) = (P1 - P3).transpose();
	
	Eigen::Vector3d vec;
	vec.x() = 0.5*(P2 - P1).squaredNorm() + (P2 - P1).dot(P1);
	vec.y() = 0.5*(P3 - P2).squaredNorm() + (P3 - P2).dot(P2);
	vec.z() = 0.5*(P1 - P3).squaredNorm() + (P1 - P3).dot(P3);
	
	return mat.colPivHouseholderQr().solve(vec);
}

void algorithm::skeletonization::propagation::reorder(const std::unordered_map<unsigned int,std::tuple<Eigen::Vector2d, unsigned int, unsigned int> > &bndopti,
														  const std::list<unsigned int> &linds,
														  const Eigen::Vector2d &C,
														  std::vector<unsigned int> &vinds)
{
	std::list<std::pair<unsigned int,double> > lind_ang;
	for(std::list<unsigned int>::const_iterator it = linds.begin(); it != linds.end(); it++)
	{
		Eigen::Vector2d vec = std::get<0>(bndopti.find(*it)->second) - C;
		double ang = atan2(vec.y(),vec.x());
		lind_ang.push_back(std::make_pair(*it,ang));
	}

	lind_ang.sort([](const std::pair<unsigned int,double> &lp1, const std::pair<unsigned int,double> &lp2)
			{
				return lp1.second < lp2.second;
			});
	
	vinds.resize(0);
	vinds.reserve(linds.size());
	for(std::list<std::pair<unsigned int,double> >::iterator it = lind_ang.begin(); it != lind_ang.end(); it++)
	{
		vinds.push_back(it->first);
	}
}
