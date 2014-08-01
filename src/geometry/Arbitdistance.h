#ifndef __ARBIT_DISTANCE_H__
#define __ARBIT_DISTANCE_H__

// This is a modif of the l2distance to adapt to arbitrary pairewise input distance (Fred Oct 26 2012)

#include <utilities/types.h>
#include <utilities/log.h>

#include <vector>
#include <fstream>
#include <functional>
#include <cmath>
#include <string>
#include <sstream>


typedef     std::vector<double>                                     Point; //Indeed a point store its index in Point[0] and its distance to the other points in the other dimensions of the vector. 
																		   // As a consequence the points indexing starts at 1.
typedef     std::vector<Point>                                      PointContainer;

struct ArbitDistance:
    public std::binary_function<const Point&, const Point&, double>
{
    result_type     operator()(const Point& p1, const Point& p2) const
    {
        AssertMsg(p1.size() == p2.size(), "Points must be in the same dimension (in Arbitrary Distance): dim1=%d, dim2=%d", p1.size(), p2.size());
        result_type sum = 0;
		int index_p2 = (int) p2[0];
		sum = p1[index_p2];
		return sum;
        //for (size_t i = 0; i < p1.size(); ++i)
        //    sum += (p1[i] - p2[i])*(p1[i] - p2[i]);

        //return sqrt(sum);
    }
};

void    read_points2(const std::string& infilename, PointContainer& points)
{
    std::ifstream in(infilename.c_str());
    std::string   line;
	int index = 1;
    while(std::getline(in, line))
    {
        if (line[0] == '#') continue;               // comment line in the file
        std::stringstream linestream(line);
        double x;
        points.push_back(Point());
		points.back().push_back(index);
        while (linestream >> x)
            points.back().push_back(x);
		index++;
    }
}

#endif // __ARBIT_DISTANCE_H__
