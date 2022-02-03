/***********************************************************************
 * Software License Agreement (BSD License)
 *
 * Copyright 2018-2021  Brian R. Phung. All rights reserved.
 * Copyright 2018-2021  Junyan He. All rights reserved.
 * Copyright 2018-2021  Ashley Spear (ashley.spear@utah.edu). 
 *   All rights reserved.
 *
 * THE BSD LICENSE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *************************************************************************/


#include "asciiSTL.h"

namespace std {

asciiSTL::asciiSTL() {
	// TODO Auto-generated constructor stub

}

int asciiSTL::writeIndividualFacets(string folder, GlobalMesh& mesh) {
	// Create an inverse map from the mesh
	map<int, set<int> > inverse_map;
	map<int, set<int> >::iterator it;
	for (it = mesh.facet2grain.begin(); it != mesh.facet2grain.end(); it++) {

		//cout << "FIRST: " <<  it->first << endl;
		set<int> &grainIDs = it->second;
		set<int>::iterator i;
		for (i = grainIDs.begin(); i != grainIDs.end(); ++i) {
			// For our negative modifications we made for the refinement zone, where we subtracted 1 from negative values due to DREAM3D overwriting -1
			// We need to write out the original final stl names
			int modifiedGrainID;
			if (*i < 0) {
				modifiedGrainID = *i - 1;
			} else {
				modifiedGrainID = *i;
			}
			// If Grain ID is in the inverse map, add this facet into the set
			if (inverse_map.count(modifiedGrainID) > 0){
				inverse_map.at(modifiedGrainID).insert(it->first);
			} else {
				set<int> grains;
				grains.insert(it->first);
				inverse_map.insert(pair<int, set<int> >(modifiedGrainID, grains));
			}
		}

	}

	// Cool, now let's go through each grain and write out the STLs
	map<int, set<int> >::iterator it2;
	for (it2 = inverse_map.begin(); it2 != inverse_map.end(); it2++) {
		// Define our new STL name
		string stlout = "indgrain_" + to_string(it2->first) + ".stl";
		ofstream outfile(folder + "/" + stlout);
		outfile << std::setprecision(16);
		outfile << "solid\n";
		set<int> &facetIDs = it2->second;
		set<int>::iterator i;
		for (i = facetIDs.begin(); i != facetIDs.end(); ++i) {

			try{
#if (STLDEBUG > 0)
				cout << "Trying to find element " << *i << endl;
#endif
			
				vector<int> facetNodes = mesh.globalFacets.at(*i);

				// Write a junk facet normal for now
				outfile << " facet normal 0 0 "  << *i << "\n";
				outfile << "  outer loop\n";
#if (STLDEBUG > 1)
				cout << "Trying to find node " << facetNodes[0] << endl;
#endif
				vector<double> xyz1 = mesh.nodes.at(facetNodes[0]);
#if (STLDEBUG > 1)
				cout << "Trying to find node " << facetNodes[1] << endl;
#endif
				vector<double> xyz2 = mesh.nodes.at(facetNodes[1]);
#if (STLDEBUG > 1)
				cout << "Trying to find node " << facetNodes[2] << endl;
#endif
				vector<double> xyz3 = mesh.nodes.at(facetNodes[2]);
#if (STLDEBUG > 1)
				cout << "Done with it all." << endl;
#endif
				// Write vertexes
				outfile << "   vertex " << xyz1[0] << " " << xyz1[1] << " " << xyz1[2] << "\n";
				outfile << "   vertex " << xyz2[0] << " " << xyz2[1] << " " << xyz2[2] << "\n";
				outfile << "   vertex " << xyz3[0] << " " << xyz3[1] << " " << xyz3[2] << "\n";

				// Close out

				outfile << "  endloop\n";
				outfile << " endfacet\n";

			}
			catch(...){
#if (DEBUG > 0)
				cout << "Oops... this facet ID is not valid..." << endl;
#endif
			}

		}
		//
		outfile << "endsolid\n";
		outfile.close();

	}

	//cout << "Written grains." << endl;
	return 0;



}

int asciiSTL::writeGlobalFacets(GlobalMesh& mesh, string stlout) {

		// Define our new STL name
		ofstream outfile(stlout);
		//outfile << std::setprecision(16);
		outfile << "solid\n" << endl;

		typedef map<int, vector<int> >::iterator it_type;
		for (it_type it = mesh.globalFacets.begin(); it != mesh.globalFacets.end(); it++) {
			int facet_id = it->first;
			vector<int> facetNodes = it->second;

			// Write a junk facet normal for now
			outfile << " facet normal 0 0 0\n";
			outfile << "  outer loop\n";

			// Grab xyz from facet node 1
			vector<double> xyz1 = mesh.nodes.at(facetNodes[0]);
			vector<double> xyz2 = mesh.nodes.at(facetNodes[1]);
			vector<double> xyz3 = mesh.nodes.at(facetNodes[2]);

			// Write vertexes
			outfile << "   vertex " << xyz1[0] << " " << xyz1[1] << " " << xyz1[2] << "\n";
			outfile << "   vertex " << xyz2[0] << " " << xyz2[1] << " " << xyz2[2] << "\n";
			outfile << "   vertex " << xyz3[0] << " " << xyz3[1] << " " << xyz3[2] << "\n";

			// Close out

			outfile << "  endloop\n";
			outfile << " endfacet\n";

		}
		outfile << "endsolid\n";
		outfile.close();
		return 0;
}

asciiSTL::~asciiSTL() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
