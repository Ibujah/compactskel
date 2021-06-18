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

#ifndef _SKELETONFILE_H_
#define _SKELETONFILE_H_

#include <skeleton/Skeletons.h>
#include <mathtools/affine/Frame.h>

/**
 *  \brief Input/output functions
 */
namespace fileio
{
    void writeSkeleton(const skeleton::GraphSkel2d::Ptr grskel, 
                       const std::map<unsigned int, std::vector<unsigned int> >& deldat, 
                       const std::string& nodname, const std::string& edgname, const std::string& delname);
}

#endif //_SKELETONFILE_H_
