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
 *  \file BoundaryFile2D.cpp
 *  \brief Defines 2D boundary file reader
 *  \author Bastien Durix
 */

#include "BoundaryFile2D.h"

#include <string>
#include <fstream>
#include <sstream>

#include <iostream>

boundary::DiscreteBoundary<2>::Ptr fileio::ReadBoundary2D(const std::string &filename)
{
	std::ifstream file(filename);
	boundary::DiscreteBoundary<2>::Ptr bnd;
	
	if(file)
	{
		bnd = std::make_shared<boundary::DiscreteBoundary<2> >();
		std::list<Eigen::Vector2d> l_pt2d;

		while(!file.eof())
		{
			double x,y;
			file >> x >> y;

			l_pt2d.push_back(Eigen::Vector2d(x,y));
		}
		bnd->addVerticesVector(l_pt2d);
	}
	
	return bnd;
}

void fileio::WriteBoundary2D(const boundary::DiscreteBoundary<2>::Ptr bnd, const std::string &filename)
{
	std::ofstream file(filename);

	if(file)
	{
		file << "\%points" << std::endl;
		for(unsigned int i = 0; i < bnd->getNbVertices(); i++)
		{
			Eigen::Vector2d pt = bnd->getCoordinates(i);
			file << pt.x() << " " << pt.y() << std::endl;
		}
		file << std::endl << std::endl;
		
		file << "\%edges" << std::endl;
		for(unsigned int i = 0; i < bnd->getNbVertices(); i++)
		{
			file << i << " " << bnd->getNext(i) << std::endl;;
		}
		file << std::endl;

		file.close();
	}

}
