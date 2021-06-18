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
 *  \file SkeletonFile.h
 *  \brief Input/output functions for skeletons
 *  \author Bastien Durix
 */

#include "SkeletonFile.h"
#include <fstream>

namespace fileio
{

void writeSkeleton(const skeleton::GraphSkel2d::Ptr grskel, 
                   const std::map<unsigned int, std::vector<unsigned int> >& deldat, 
                   const std::string& nodname, const std::string& edgname, const std::string& delname)
{
    std::ofstream fnod(nodname);
    std::ofstream fedg(edgname);
    std::ofstream fdel(delname);

    if(!fnod || !fedg || !fdel)
    {
        std::cout << "Problem" << std::endl;
        return;
    }

	std::list<unsigned int> nods;
	grskel->getAllNodes(nods);
	for(std::list<unsigned int>::iterator it = nods.begin(); it != nods.end(); it++)
	{
		mathtools::geometry::euclidian::HyperSphere<2> cir1 = grskel->getNode<mathtools::geometry::euclidian::HyperSphere<2> >(*it);
		Eigen::Vector2d pt = cir1.getCenter().getCoords();
        double rad = cir1.getRadius();
        
        fnod << pt.x() << " " << pt.y() << " " << rad << std::endl;
	}

	std::list<std::pair<unsigned int,unsigned int> > edges;
	grskel->getAllEdges(edges);
	for(std::list<std::pair<unsigned int,unsigned int> >::iterator it = edges.begin(); it != edges.end(); it++)
	{
        fedg << it->first << " " << it->second << std::endl;
	}
    

	for(auto it : deldat)
	{
        fdel << it.first << " : ";
        for(auto itn : it.second)
            fdel << itn << ", ";
        fdel << std::endl;
	}

    fnod.close();
    fedg.close();
    fdel.close();
}

}
