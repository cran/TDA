#ifndef __DIONYSUSUTILS_H__
#define __DIONYSUSUTILS_H__

#include <topology/simplex.h>
#include <string>
#include <sstream>
#include <cstdlib>

std::vector<unsigned int> getVertices(const Simplex<unsigned, double> &smp)
{
	std::stringstream sstr;
	std::vector<unsigned int> vertices(smp.dimension()+1);
	sstr << smp;
	std::string vtxStr;
	std::getline(sstr,vtxStr,'<');
	unsigned int vtxIdx;
	for (vtxIdx=0; vtxIdx<smp.dimension();++vtxIdx)
	{
		std::getline(sstr,vtxStr,',');
		vertices[vtxIdx] = (unsigned)std::atoi(vtxStr.c_str());
	}
	std::getline(sstr,vtxStr,'>');
	vertices[vtxIdx] = (unsigned)std::atoi(vtxStr.c_str());
	return vertices;
}

unsigned int getLocation(const Simplex<unsigned, double> &smp, const double * const FUNvaluesInput)
{
	std::vector< unsigned > vertices;
	std::vector< unsigned >::const_iterator vertexItr;
	unsigned int vertex;
	vertices = getVertices( smp );
	vertex = *(vertices.begin());
	for (vertexItr = vertices.begin(); vertexItr != vertices.end(); ++vertexItr)
	{
		if (FUNvaluesInput[ *vertexItr ] > FUNvaluesInput[ vertex ])
		{
			vertex = *vertexItr;
		}
	}
	return vertex+1;
}

# endif // __DIONYSUSUTILS_H__
