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


#ifndef MESHIO_GLOBALMESH_H_
#define MESHIO_GLOBALMESH_H_

#include "debug.h"

#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

#include "nanoflann.hpp"
#include "Grain.h"
#include "refine/coarsenSurf.h"
#include "refine/refineSurf.h"
#include "programSettings.h"

using namespace nanoflann;
namespace std {

class GlobalMesh {
public:
	GlobalMesh(programSettings *settings_in);
	virtual ~GlobalMesh();

	// Shrink Factors
	// TODO: Make private
	double Undo_shrink_factorX = 1., Undo_shrink_factorY = 1., Undo_shrink_factorZ = 1.;

	programSettings *settings;

	// Structure of Mesh
	// TODO: Make maps private
	map<int, vector<int> > elements;
	map<int, vector<double> > nodes;
	map<int, vector<int> > globalFacets;
	map<int, set<int> > facet2grain;
	map<int, vector<int> > elementSets;
	int addToFacet2Grain(int facetID, int grainID);
	int addToFacet2Grain(int facetID, set<int> grainID);
	int removeFromFacet2Grain(int facetID, int grainID);
	
	// Process Functions
	int processGrainSurfaceNodes(Grain& grain);
	int processGrainElements(Grain& grain, int &count_nodes, int &count_els);
	int processGrainFacets(Grain &grain, int grainID);

	
	void BuildCrackFrontKD( string FilePath, ofstream &logfile );
	int buildEdgeList(int mode);

	int calculateFacetNormal(vector<int> &nids, vector<double>& normal);
	int calculateFacetAngles(vector<int> &nids, vector<double>& angles);

	// Main coarsen and refine functions
	int coarsen( int rr , ofstream &outfile);
	int refine(ofstream &outfile);
	int refine_iterations = 0;

	// I/O
	void writeMTR( string filename );

	
private:

	int max_gid = 0;
	int max_facet_id;
	int kd_node_id;
	int max_eid;

	// Modification Counters
	set<int> obliteratedNodes;
	set<int> modifiedFacets;
	set<int> facetSplitStatus;
	int resetFaceSplitStatuses() {
		facetSplitStatus.clear();
		return 0;
	}


	// Edges
	struct edge{
		vector<int> n;
		double length;
	};
	vector<edge> edges;
	vector<edge> jEdges;

	// Classification arrays
	set<int> interfaceNodes;
	set<int> BoundaryFacets;

	// Element Sets and Sections
	map<int, int > sections;
	vector<int> materials;
	map< int , double > VoidRadius;

	// Access Maps
	map<int, vector<int> > node2facet;
	map<int, vector<int> > nodeNeighbors;

	// Utility Functions
	double calculateAngleBetweenNormals(vector<double>& u, vector<double>& v);

	// Functions Related to Refinement
	int splitEdge3(vector<int> edge, int nid1, int nid2, double normalTol, double chordTol, double minAngle, double maxAngle);
	int splitAnEdge(int nid1, int nid2, vector<Grain>& grains);

	// Functions Related to Coarsening
	int collapseEdge3( int MODE , vector<int> edge, int nid1, int nid2, double normalTol, double chordTol, double minAngle, double maxAngle );
	int collapseEdge( vector<int> edge, int nid1, int nid2, double normalTol, double chordTol, double minAngle, double maxAngle );

	// Swap Functions
	int multiSwap(vector<vector<int> > proposedNodes, map<int, vector<int> > &proposedFacets, set<int> &relevantNodes, double normalTol, double minAngle);
	int swapPlanar(int nid1, int nid2, vector<int> &facet1, vector<int> &facet2, double normalTol, double minAngle);
	
	// Genreal Utility Functions
	void addFacet2NodeMap(int nodeID, int globalFid);
	void deleteFacetInNodeMap(int nodeID, int globalFid);
	void addNode2NodeNeighbors(int nodeMaster, int nodeSlave);
	void GetSharingEdge( vector<int> &facet_pair , vector<int> &SE );
	int ContainsEdge( vector<int> &SE , vector<int> &facet_pair , vector<int> &target );
	void GetOpposite( vector<int> &facet_pair , vector<int> &fp1 , int ref_node , vector<int> &target_facet );
	int Overlap( vector<int> &FP1 , vector<int> &FP2 );
	void UpdateFacetPair( vector< vector<int> > &facet_pair , vector<int> &fp_new , int idx );
	int CheckNormal( vector<int> &auxFacets , vector<double> &F1 , int nid1 , int nid2 , double normalTol);
	int ExternalLoop( vector<int> &auxFacets , vector<int> &External_Loop , int nid1 , int nid2 , int MODE );
	int FindConcavity( vector<int> &External_Loop , vector<double> &F1 , int num_edges , vector<int> &concave_loc );
	int SwapFacets( vector< vector<int> > &facet_pair , vector< vector<int> > &initial_facet_pair , 
		vector<int> &concave_loc , double minAngle , double maxAngle , double normalTol , vector< vector<int> > &final_facet_pair );
	void UpdateMesh( vector<int> &auxFacets , map<int, vector<int> > &proposedFacets , int nid1 , int nid2 );
	int sizingFunction( double currLen , vector<double> &mp , double& Objective_Length , int Mode );
	void UnifyConvention( int REF , int fid );

	// Mesh Preprocessing Functions
	void CONDENSE3(ofstream &outfile);
	void CONDENSE4( double minAngle , double maxAngle );
	void FixHourglass(ofstream &outfile);
	void CheckDuplicate();
	void RemoveFreeEdge( vector< vector<int> > &tmpFacet , vector<int> &External_Loop , int targetSize );
	bool CheckOverlap( vector< vector<int> > &tmpFacet , vector<int> &External_Loop , vector<int> &AllInternalFacets );
	void FixBoundary(ofstream &outfile);

	// Node Point Cloud
	struct nodeCloudPt {
		vector<vector<double> > nodexyz;
		vector<int> nodeID;

		int kdtree_get_point_count() const {
			return nodexyz.size();
		}

		unsigned long int getSize() {
			return (unsigned long int)(nodexyz.size());
		}

		inline double kdtree_get_pt(const size_t idx, int dim) const {
			return nodexyz[idx][dim];
		}

		template<class BBOX>
		bool kdtree_get_bbox(BBOX& /* bb */) const {
			return false;
		}
	};
	nodeCloudPt cloud;
	typedef KDTreeSingleIndexDynamicAdaptor<L2_Simple_Adaptor<double, nodeCloudPt>  ,
			nodeCloudPt,
			3 /* dim */
	> my_kd_tree_t;
	my_kd_tree_t index;
	// Node Point Cloud

	// Second KD-tree for facets
	nodeCloudPt facecloud;
	my_kd_tree_t facedex;

	// Third KD-tree for Crack Points
	nodeCloudPt CrackPtcloud;
	my_kd_tree_t CrackPtKD;
	

	vector<double> calculateMidpointBetweenTwoDoubleVectors(vector<double> v1, vector<double> v2) {
		vector<double> midpoint;
 		midpoint.push_back((v1[0] + v2[0])/2);
		midpoint.push_back((v1[1] + v2[1])/2);
		midpoint.push_back((v1[2] + v2[2])/2);
		return midpoint;
	}
	double calculateDistanceBetweenTwoDoubleVectors(vector<double> v1, vector<double> v2) {
		return sqrt( (v1[0] - v2[0])*(v1[0] - v2[0]) + (v1[1] - v2[1])*(v1[1] - v2[1]) + (v1[2] - v2[2])*(v1[2] - v2[2]) );
	}

	int isThisADuplicate(vector<double>& xyz, int& duplicateID, my_kd_tree_t &idx);


};

} /* namespace std */

#endif /* MESHIO_GLOBALMESH_H_ */
