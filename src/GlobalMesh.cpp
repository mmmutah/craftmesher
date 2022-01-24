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


#include "GlobalMesh.h"

namespace std {

bool debug(int nid) {
	return false;
}

bool debug(vector<int> edge) {
	return false;
}

GlobalMesh::GlobalMesh(programSettings *settings_in) :
		index(3 /*dim*/, cloud,
				KDTreeSingleIndexAdaptorParams(10 /* max leaf */)), facedex(
				3 /*dim*/, facecloud,
				KDTreeSingleIndexAdaptorParams(10 /* max leaf */)), CrackPtKD(
				3 /*dim*/, CrackPtcloud,
				KDTreeSingleIndexAdaptorParams(10 /* max leaf */)) {
	max_eid = 1, kd_node_id = 1, max_facet_id = 1;
	settings = settings_in;
	Undo_shrink_factorX = settings->shrinkFactorX;
	Undo_shrink_factorY = settings->shrinkFactorY;
	Undo_shrink_factorZ = settings->shrinkFactorZ;
}

void GlobalMesh::addFacet2NodeMap(int nodeID, int globalFid) {
	std::map<int, vector<int> >::iterator it = node2facet.find(nodeID);
	if (it != node2facet.end()) {
		// Make sure the value isn't already in the map
		if (std::find(it->second.begin(), it->second.end(), globalFid)
				!= it->second.end()) {
			/* v contains x */
		} else {
			/* v does not contain x */
			it->second.push_back(globalFid);
		}

	} else {
		vector<int> values;
		values.push_back(globalFid);
		node2facet.insert(pair<int, vector<int> >(nodeID, values));
	}
}

int GlobalMesh::addToFacet2Grain(int facetID, set<int> grainID) {
	std::map<int, set<int> >::iterator it = facet2grain.find(facetID);
	if (it != facet2grain.end()) {
		// It's in this set, just insert it
		facet2grain.at(facetID).insert(grainID.begin(), grainID.end());
	} else {
		set<int> temp;
		temp.insert(grainID.begin(), grainID.end());
		facet2grain.insert(pair<int, set<int> >(facetID, temp));
	}
	return 1;

}

int GlobalMesh::addToFacet2Grain(int facetID, int grainID) {
	std::map<int, set<int> >::iterator it = facet2grain.find(facetID);
	if (it != facet2grain.end()) {
		// It's in this set, just insert it
		facet2grain.at(facetID).insert(grainID);
	} else {
		set<int> temp;
		temp.insert(grainID);
		facet2grain.insert(pair<int, set<int> >(facetID, temp));
	}
	return 1;

}
int GlobalMesh::removeFromFacet2Grain(int facetID, int grainID) {
	facet2grain.at(facetID).erase(grainID);
	return 1;
}

void GlobalMesh::addNode2NodeNeighbors(int nodeMaster, int nodeSlave) {
	std::map<int, vector<int> >::iterator it = nodeNeighbors.find(nodeMaster);
	if (it != nodeNeighbors.end()) {
		// Make sure the value isn't already in the map
		if (std::find(it->second.begin(), it->second.end(), nodeSlave)
				!= it->second.end()) {
			/* v contains x */
		} else {
			/* v does not contain x */
			it->second.push_back(nodeSlave);
		}

	} else {
		vector<int> values;
		values.push_back(nodeSlave);
		nodeNeighbors.insert(pair<int, vector<int> >(nodeMaster, values));
	}
}

void GlobalMesh::deleteFacetInNodeMap(int nodeID, int globalFid) {
	int nidlistindex;
	std::map<int, vector<int> >::iterator it = node2facet.find(nodeID);

	if (it != node2facet.end()) {
		// Check if value is in the map
		vector<int>::iterator it2 = std::find(it->second.begin(),
				it->second.end(), globalFid);
		if (it2 != it->second.end()) {
			/* v contains x */
			nidlistindex = std::distance(it->second.begin(), it2);
		} else {
			/* v does not contain x */
			// We're good
			return;
		}

	} else {
		return;
	}

	it->second.erase(it->second.begin() + nidlistindex);

	return;
}

int GlobalMesh::isThisADuplicate(vector<double> &xyz, int &duplicateID,
		my_kd_tree_t &idx) {

	double query_pt[3] = { xyz[0], xyz[1], xyz[2] };
	const size_t num_results = 1;
	size_t ret_index;
	double out_dist_sqr;
	KNNResultSet<double> resultSet(num_results);
	resultSet.init(&ret_index, &out_dist_sqr);
	idx.findNeighbors(resultSet, query_pt, nanoflann::SearchParams(10));

	duplicateID = ret_index;
	if (out_dist_sqr < std::numeric_limits<double>::epsilon()) {
		return 1;

	} else {
		return 0;
	}

	return 0;
}

int GlobalMesh::buildEdgeList(int mode) {
	// Mode 0: do all the edges
	// Mode -1: do edges only associated with -1

	set<int> nodesEnumerated;

	// Clear our edge list, rebuild from scratch
	edges.clear();
	jEdges.clear();
	map<int, vector<int>>::iterator it;

	// Loop through all the nodes
	for (it = nodeNeighbors.begin(); it != nodeNeighbors.end(); it++) {
		int nid = it->first;

		vector<int> connectedNodes = it->second;

		// Loop through all the node's neighbors
		for (uint i = 0; i < connectedNodes.size(); ++i) {
			// IF the nodes have already been enumerated, skip this iteration
			if (nodesEnumerated.count(connectedNodes[i])) {
				continue;
			} // Else, add an edge to the edge list
			else {
				int &nid1 = nid;
				int &nid2 = connectedNodes[i];

				vector<int> &facetlist1 = node2facet.at(nid1);
				vector<int> &facetlist2 = node2facet.at(nid2);
				// Get length of facet
				double dist = calculateDistanceBetweenTwoDoubleVectors(
						nodes.at(nid), nodes.at(connectedNodes[i]));

				// Matching facets in the list are the facets that share an edge
				vector<int> matchingFacets, auxFacets;

				for (uint i = 0; i < facetlist2.size(); ++i) {
					if (std::find(facetlist1.begin(), facetlist1.end(),
							facetlist2[i]) != facetlist1.end()) {
						matchingFacets.push_back(facetlist2[i]);
					} else {

						auxFacets.push_back(facetlist2[i]);
					}
				}

				bool continueFunction = false;
				// Check our mode:
				switch (mode) {
				case 0: {
					continueFunction = true;
					break;
				}
				case -1: {
					// First, get the GIDs associated with this node
					for (uint i = 0; i < matchingFacets.size(); ++i) {
						set<int> GIDs = facet2grain.at(matchingFacets[i]);
						for (auto gid : GIDs) {
							if (gid < 0) {
								// If there is a GID less than 0, continue with the function
								continueFunction = true;
								break;
							}
						}

					}
					break;
				}
				default: {
					break;
				}
				}
				// Skip the rest if the bool is false
				if (matchingFacets.size() == 0) {
					// if (debug(nid1) && debug(nid2)) {
					if (1) {
						cout
								<< ">>> ERROR: Found zero matching facets in collapsingEdge function for nids "
								<< nid1 << ", " << nid2 << "." << endl;
					}
					//return -1;
					exit(1);
				} else if (matchingFacets.size() == 1) {
					if (1) {
						cout
								<< ">>> ERROR: Found only one matching facet in collapsingEdge function for nids "
								<< nid1 << ", " << nid2
								<< ". This should not happen for manifold geometries."
								<< endl;
					}
					//return -1;
					exit(1);
				} else if ((matchingFacets.size() > 2)
						&& (continueFunction == true)) {
					if (debug(nid1) && debug(nid2)) {
						cout << ">>> DEBUG OUTPUT: Found "
								<< matchingFacets.size()
								<< " matching facets in collapsingEdge function for nids "
								<< nid1 << ", " << nid2 << "." << endl;
					}
					edge temp;
					temp.n = {nid, connectedNodes[i]};
					temp.length = dist;
					jEdges.push_back(temp);
				} else if (continueFunction == true) {
					edge temp;
					temp.n = {nid, connectedNodes[i]};
					temp.length = dist;
					edges.push_back(temp);
				}

			}

		}

		// Add this node into the nodesEnumerated set
		nodesEnumerated.insert(nid);

	}

	return 1;
}

int GlobalMesh::processGrainSurfaceNodes(Grain &grain) {
	for (uint i = 0; i < grain.surfaceNodes.size(); i++) {
		// Extract node from grain
		int closestNID = 0;
		int grainBasisNID = grain.surfaceNodes[i];
		vector<double> xyz = grain.localNodes.at(grainBasisNID);
		int duplicate = isThisADuplicate(xyz, closestNID, index);

		if (duplicate == 0) {
			// If Not Duplicate, place in KD-tree
			cloud.nodeID.push_back(kd_node_id);
			cloud.nodexyz.push_back(xyz);
			index.addPoints(cloud.getSize() - 1, cloud.getSize() - 1);
			// Set in KD-tree map in grain
			grain.kdTreeMap.insert(pair<int, int>(grainBasisNID, kd_node_id));
			// Add to node database:
			nodes.insert(pair<int, vector<double> >(kd_node_id, xyz));
			// Increment kd_node_id
			kd_node_id++;
		} else {
			// Otherwise... Grab the ret index and throw it in the grain's node ID map
			grain.kdTreeMap.insert(
					pair<int, int>(grainBasisNID, closestNID + 1));
		}

	}
	return kd_node_id;
}

int GlobalMesh::processGrainFacets(Grain &grain, int indexID) {
	map<int, vector<int>>::iterator it;
	for (it = grain.localFacets.begin(); it != grain.localFacets.end(); it++) {

		// I have facet IDs and NIDs
		int fid = it->first;
		vector<int> nids = it->second;

		int closestfid = 0;

		// First, convert my nids to global NIDs
		vector<int> globalNIDs;
		globalNIDs.push_back(grain.kdTreeMap.at(nids[0]));
		globalNIDs.push_back(grain.kdTreeMap.at(nids[1]));
		globalNIDs.push_back(grain.kdTreeMap.at(nids[2]));

		vector<double> xyz0 = grain.localNodes.at(nids[0]);
		vector<double> xyz1 = grain.localNodes.at(nids[1]);
		vector<double> xyz2 = grain.localNodes.at(nids[2]);

		// Build a facet centroid, used in the KD tree search
		vector<double> xyz =
				{ (xyz0[0] + xyz1[0] + xyz2[0]) / 3.0, (xyz0[1] + xyz1[1]
						+ xyz2[1]) / 3.0, (xyz0[2] + xyz1[2] + xyz2[2]) / 3.0 };

		// Check if this facet matches any facets I currently have in my facet list
		int duplicate = isThisADuplicate(xyz, closestfid, facedex);
		if (duplicate == 0) {
			// If Not Duplicate, place in KD-tree
			facecloud.nodeID.push_back(max_facet_id);
			facecloud.nodexyz.push_back(xyz);
			facedex.addPoints(facecloud.getSize() - 1, facecloud.getSize() - 1);
			// Set in KD-tree map in grain
			globalFacets.insert(
					pair<int, vector<int> >(max_facet_id, globalNIDs));
			for (uint n = 0; n < globalNIDs.size(); ++n) {
				addFacet2NodeMap(globalNIDs[n], max_facet_id);
			}

			// Populate Node Neighbors list
			addNode2NodeNeighbors(globalNIDs[0], globalNIDs[1]);
			addNode2NodeNeighbors(globalNIDs[0], globalNIDs[2]);
			addNode2NodeNeighbors(globalNIDs[1], globalNIDs[0]);
			addNode2NodeNeighbors(globalNIDs[1], globalNIDs[2]);
			addNode2NodeNeighbors(globalNIDs[2], globalNIDs[0]);
			addNode2NodeNeighbors(globalNIDs[2], globalNIDs[1]);

			grain.globalToGrainFacets.insert(pair<int, int>(max_facet_id, fid));
			grain.grainToGlobalFacets.insert(pair<int, int>(fid, max_facet_id));

			addToFacet2Grain(max_facet_id, grain.giveRawGrainNumber());
			max_facet_id += 1;

		} else {
			// Otherwise... Grab the ret index and throw it in the grain's node ID map
			grain.globalToGrainFacets.insert(
					pair<int, int>(closestfid + 1, fid));
			grain.grainToGlobalFacets.insert(
					pair<int, int>(fid, closestfid + 1));
			addToFacet2Grain(closestfid + 1, grain.giveRawGrainNumber());
		}

	}
	return max_facet_id;
}


int GlobalMesh::processGrainElements(Grain &grain, int &count_nodes,
		int &count_els) {
	// Run through all the elements in the grain
	map<int, vector<int>>::iterator it;

	for (it = grain.localElements.begin(); it != grain.localElements.end();
			it++) {
		int eid = it->first;
		vector<int> nids = it->second;

		// Go through the nids. If they exist in the grain's local KD map, replace it:
		for (uint i = 0; i < nids.size(); i++) {
			int originalNID = nids[i];
			int surfaceNodeBool = grain.getKDMap(originalNID, nids[i]);
			// if the surfacenodebool is zero, we have an internal node:
			if (surfaceNodeBool == 0) {
				// GLOBAL NID INCREMENTS IN FUNCTION:
				int internalNodeBool = grain.inInternalNodeMap(originalNID,
						nids[i], kd_node_id);
				if (internalNodeBool == 0) {
					vector<double> xyz = grain.localNodes.at(originalNID);
					nodes.insert(pair<int, vector<double> >(kd_node_id, xyz));
					kd_node_id++;
				}

			}
		}

		// Great, now add element to map:
		elements.insert(pair<int, vector<int> >(max_eid, nids));

		// Add this element into the corresponding set:
		if (elementSets.find(grain.givePositiveGrainNumber())
				== elementSets.end()) {
			// not found
			vector<int> elementSet { max_eid };
			elementSets.insert(
					pair<int, vector<int> >(grain.givePositiveGrainNumber(),
							elementSet));
		} else {
			// found
			elementSets.at(grain.givePositiveGrainNumber()).push_back(max_eid);
		}

		// Add this element into the corresponding set:
		if (sections.find(grain.givePositiveGrainNumber()) == sections.end()) {
			// not found
			int material = grain.giveBaseGrainNumber();
			sections.insert(
					pair<int, int>(grain.givePositiveGrainNumber(), material));

			if (find(materials.begin(), materials.end(), material)
					== materials.end()) {
				materials.push_back(material);
				if (material > max_gid) {
					max_gid = material;
				}
			}
		}

		max_eid++;
	}
	count_nodes = kd_node_id;
	count_els = max_eid;
	return 0;
}


int GlobalMesh::calculateFacetNormal(vector<int> &nids,
		vector<double> &normal) {

	// Grab the xyzs of the facet nodes
	vector<double> &xyz1 = nodes.at(nids[0]);
	vector<double> &xyz2 = nodes.at(nids[1]);
	vector<double> &xyz3 = nodes.at(nids[2]);

	// Separate into different variables
	double &x1 = xyz1[0];
	double &y1 = xyz1[1];
	double &z1 = xyz1[2];
	double &x2 = xyz2[0];
	double &y2 = xyz2[1];
	double &z2 = xyz2[2];
	double &x3 = xyz3[0];
	double &y3 = xyz3[1];
	double &z3 = xyz3[2];

	// Calculate normal vector components
	double n1 = (y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1);
	double n2 = (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1);
	double n3 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
	double nsum = sqrt(n1 * n1 + n2 * n2 + n3 * n3);

	normal[0] = n1 / nsum;
	normal[1] = n2 / nsum;
	normal[2] = n3 / nsum;

	return 1;
}

int GlobalMesh::calculateFacetAngles(vector<int> &nids,
		vector<double> &angles) {
	// Grab the xyzs of the facet nodes
	vector<double> &xyz1 = nodes.at(nids[0]);
	vector<double> &xyz2 = nodes.at(nids[1]);
	vector<double> &xyz3 = nodes.at(nids[2]);

	// Separate into different variables
	double &x1 = xyz1[0];
	double &y1 = xyz1[1];
	double &z1 = xyz1[2];
	double &x2 = xyz2[0];
	double &y2 = xyz2[1];
	double &z2 = xyz2[2];
	double &x3 = xyz3[0];
	double &y3 = xyz3[1];
	double &z3 = xyz3[2];

	// Vector from n1 to n2
	vector<double> n1_2 = { x2 - x1, y2 - y1, z2 - z1 };
	// Vector from n1 to n3
	vector<double> n1_3 = { x3 - x1, y3 - y1, z3 - z1 };
	double n1_2_mag = sqrt(
			n1_2[0] * n1_2[0] + n1_2[1] * n1_2[1] + n1_2[2] * n1_2[2]);
	double n1_3_mag = sqrt(
			n1_3[0] * n1_3[0] + n1_3[1] * n1_3[1] + n1_3[2] * n1_3[2]);
	n1_2 = {(x2 - x1) / n1_2_mag, (y2 - y1) / n1_2_mag, (z2 - z1) / n1_2_mag};
	n1_3 = {(x3 - x1) / n1_3_mag, (y3 - y1) / n1_3_mag, (z3 - z1) / n1_3_mag};

	// Vector from n2 to n1
	vector<double> n2_1 = { -n1_2[0], -n1_2[1], -n1_2[2] };

	// Vector from n2 to n3
	vector<double> n2_3 = { x3 - x2, y3 - y2, z3 - z2 };
	double n2_3_mag = sqrt(
			n2_3[0] * n2_3[0] + n2_3[1] * n2_3[1] + n2_3[2] * n2_3[2]);
	n2_3 = {(x3 - x2) / n2_3_mag, (y3 - y2) / n2_3_mag, (z3 - z2) / n2_3_mag};

	// Angles
	angles[0] = calculateAngleBetweenNormals(n1_2, n1_3);
	angles[1] = calculateAngleBetweenNormals(n2_1, n2_3);

	angles[2] = 180 - angles[1] - angles[0];

	std::sort(angles.begin(), angles.end());

	return 1;

}

double GlobalMesh::calculateAngleBetweenNormals(vector<double> &u,
		vector<double> &v) {
	double intermediate = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
	double u_norm = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
	double v_norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	intermediate = intermediate / (u_norm * v_norm);

	if (intermediate < 0) {
		intermediate += 1e-6;
	} else if (intermediate > 1) {
		intermediate -= 1e-6;
	}
	return acos(intermediate) * 180.0 / 3.14159265;
}


int GlobalMesh::coarsen(int rr, ofstream &outfile) {
	coarsenSurf coarsen;
	// Classify all edges in the surface mesh and build edge list for ALL edges.
	buildEdgeList(0);

	// Sort them based on increasing length
	std::sort(edges.begin(), edges.end(), [](edge const &a, edge const &b) {
		return a.length < b.length;
	});
	vector<edge>::iterator it;

	// Switches to precondition the mesh
	bool RemoveHourglass = settings->RemoveHourglass;
	
	bool RemoveWavyBoundary = settings->SmoothBoundary;
	
	// Remove any nodes that are only connected to 3 neighbors
	// This operation detects a Y shape inside a triangle and replaces the 3 smaller ones with the outer triangle
	bool Condense3 = settings->Condense3;
	
	// Junyan's note on 07/07/20
	// The CONDENSE4 function does not seem to work properly, don't use this switch
	// Remove any nodes that are only connected to 4 neighbors and adjust local facets
	// This operation detects an internal node within a parallelogram and attempts to replace it with 2 triangles (instead of 4)
	bool Condense4 = false;

	// Loop through all interface edges, store all nodes on the GB interfaces
	for (it = jEdges.begin(); it != jEdges.end(); it++) {
		vector<int> edge = (*it).n;
		// Store all interface nodes
		interfaceNodes.insert(edge[0]);
		interfaceNodes.insert(edge[1]);
	}

	// Step 1: Precondition the surface mesh based on user inputs 
	if (RemoveHourglass && rr == 0)
		FixHourglass(outfile);
	if (Condense3)
		CONDENSE3(outfile);
	if (RemoveWavyBoundary && rr == 0)
		FixBoundary(outfile);
	if (Condense4)
		CONDENSE4(coarsen.minAngle, coarsen.maxAngle);


	// Step 2: Loop through all boundary edges, coarsen as necessary
	int i = 0;
	int Good_Swap = 0;
	// Loop through all the boundary edges
	for (it = jEdges.begin(); it != jEdges.end(); it++) {
		vector<int> edge = (*it).n;
		#if (COARSENDEBUG > 0)
			cout << "----------- BEGIN ITERATION " << i << " -----------" << endl;
		#endif
		i++;

		// Get nodes on this edge
		int nid1 = edge[0];
		int nid2 = edge[1];

		// Check if this edge is connected to something that we have already modified
		if (obliteratedNodes.count(nid1) || obliteratedNodes.count(nid2)) {
			#if (COARSENDEBUG > 1)
				cout << "We hit a previously eliminated node, let's skip this" << endl;
			#endif
			continue;
		}

		// Attempt to coarsen 
		int success = 0;
		success = collapseEdge3(1, edge, nid1, nid2, coarsen.normalTol,
				coarsen.chordTol, coarsen.minAngle, coarsen.maxAngle);
		if (success == 0) {
			Good_Swap += 1;
		}
		if (success == -2) {
			// For debugging purpose, if it fails at weird places, write the surface mesh as is
			return -1;
		}
	}


	// Step 3: Loop through all regular edges, coarsen as necessary
	for (it = edges.begin(); it != edges.end(); it++) {
		vector<int> edge = (*it).n;
		if (debug(edge))
			cout << "----------- BEGIN ITERATION " << i << " -----------" << endl;
		i++;

		// Get nodes on this edge
		int nid1 = edge[0];
		int nid2 = edge[1];

		// Check if this edge is connected to something that we have already modified
		if (obliteratedNodes.count(nid1) || obliteratedNodes.count(nid2)) {
			#if (COARSENDEBUG > 0)
				cout << "We hit a previously eliminated node, let's skip this" << endl;
			#endif
			continue;
		}

		// Check if this edge is directly connected to a boundary node
		if (interfaceNodes.count(nid1) || interfaceNodes.count(nid2)) {
			#if (COARSENDEBUG > 0)
				cout << "This edge is directly connected to the interface, let's skip this" << endl;
			#endif
			continue;
		}

		// Attempt to coarsen 
		int success = 0;
		success = collapseEdge(edge, nid1, nid2, coarsen.normalTol,
				coarsen.chordTol, coarsen.minAngle, coarsen.maxAngle);
		if (success == 0) {
			Good_Swap += 1;
		}
		if (success == -2) {
			// For debugging purpose, if it fails at weird places, write the surface mesh as is
			return -1;
		}
	}

	// If we reach here, a coarsening pass has finished
	outfile << ">>> Coarsening pass made " << Good_Swap << " good swaps" << endl;


	// Step 4: Second sanity check, check again (this time globally) for any duplicate facets 
	outfile << ">>> Checking for duplicate facets..." << endl;
	CheckDuplicate();


	// Step 5: Clear everything
	obliteratedNodes.clear();
	modifiedFacets.clear();
	interfaceNodes.clear();

	return 0;
}


int GlobalMesh::refine(ofstream &outfile) {
	refineSurf refine;

	// DEBUG INT
	int i = 0;
	facetSplitStatus.clear();
	vector<edge>::iterator it;

	// Classify all edges in the surface mesh and build edge list for ALL edges.
	buildEdgeList(0);

	// Edge refinement can be done with a single function, so we can combine the boundary and regular edges
	vector<edge> combinedEdges = jEdges;
	combinedEdges.insert(combinedEdges.end(), edges.begin(), edges.end());

	// Sort them based on increasing lengths
	std::sort(combinedEdges.begin(), combinedEdges.end(),
			[](edge const &a, edge const &b) {
				return a.length > b.length;
			});

	// Loop through all edges
	int good_refine = 0;
	for (it = combinedEdges.begin(); it != combinedEdges.end(); it++) {
		vector<int> edge = (*it).n;
		if (debug(edge))
			cout << "----------- BEGIN ITERATION " << i << " -----------"
					<< endl;
		i++;

		int nid1 = edge[0];
		int nid2 = edge[1];

		// Calculate the edge length
		vector<double> &xyz1 = nodes.at(nid1);
		vector<double> &xyz2 = nodes.at(nid2);
		double dist = refine.calculateDistanceBetweenTwoDoubleVectors(xyz1, xyz2);
		vector<double> mid_pt = { 0.5 * (xyz1[0] + xyz2[0]), 0.5
				* (xyz1[1] + xyz2[1]), 0.5 * (xyz1[2] + xyz2[2]) };

		// Check if we should refine this edge
		int success = -1;
		double ph;
		if ( sizingFunction(dist, mid_pt, ph, 0) == -1 ) {
			// If the sizing function decides that this edge needs to be refined, attempt to refine it
			success = splitEdge3(edge, nid1, nid2, refine.normalTol, refine.chordTol, refine.minAngle, refine.maxAngle);

			if (success == -1) {
				if (debug(edge))
					cout << ">>> Edge refinement failed!" << endl;
			} else
				good_refine += 1;
		}
	}

	// If we reach here, we have finished one refinement iteration
	outfile << ">>> Refinement pass made " << good_refine << " good refinements."
			<< endl;
	return 0;
}


int GlobalMesh::splitEdge3(vector<int> edge, int nid1, int nid2,
		double normalTol, double chordTol, double minAngle, double maxAngle) {
	// Determine where the vertex will be placed
	if (debug(edge))
		cout << "****************** BEGIN of GlobalMesh splitEdge function." << endl;

	// Status
	bool split = true;

	// First, grab all the facets connected to the verticies
	vector<int> &facetlist1 = node2facet.at(nid1);
	vector<int> &facetlist2 = node2facet.at(nid2);
	if (debug(edge)) {
		cout << "Facet List for NID" << nid1 << ": ";
		for (auto q : facetlist1) {
			cout << "FID" << q << " ";
		}
		cout << endl;
		cout << "Facet List for NID" << nid2 << ": ";
		for (auto q : facetlist2) {
			cout << "FID" << q << " ";
		}
		cout << endl;

		cout << "Edge: " << nid1 << " " << nid2 << endl;
	}

	// Matching facets in the list are the facets that share an edge
	vector<int> matchingFacets, auxFacets;
	for (uint i = 0; i < facetlist2.size(); ++i) {
		// Check if this facet has been modified before
		if (modifiedFacets.count(facetlist2[i]) || modifiedFacets.count(facetlist1[i])) {
			if (debug(edge))
				cout << "Split failed due to attempt to change an already modifed facet." << endl;
			return -1;
		}

		// Classify the facets
		if (std::find(facetlist1.begin(), facetlist1.end(), facetlist2[i]) != facetlist1.end()) {
			matchingFacets.push_back(facetlist2[i]);
		} else {
			auxFacets.push_back(facetlist2[i]);
		}
	}

	// Check if what we get makes sense
	if (matchingFacets.size() < 2) {
		if (debug(edge))
			cout << "*** DEBUG OUTPUT: Found " << matchingFacets.size()
					<< " matching facets in collapsingEdge3 function for nids "
					<< nid1 << ", " << nid2 << "." << endl;
		return -1;
	}

	if (debug(edge)) {
		cout << "All Facets Being Worked With: " << endl;
		cout << "Matched: ";
		for (uint i = 0; i < matchingFacets.size(); ++i) {
			cout << matchingFacets[i] << " ";
		}
		cout << endl;
		cout << "Aux: ";
		for (uint i = 0; i < auxFacets.size(); ++i) {
			cout << auxFacets[i] << " ";
		}
		cout << endl;
	}


	// As a first attempt, place a new vertex at the midpoint of the edge
	if (debug(edge))
		cout << "Temporarily setting proposed nid." << endl;
	vector<double> splitvertexxyz = calculateMidpointBetweenTwoDoubleVectors(nodes.at(nid1), nodes.at(nid2));

	// This is the proposed id, don't increment kd_node_id until this is accepted
	int addedNID = kd_node_id;
	// Temporarily add this nid into the database
	nodes[kd_node_id] = splitvertexxyz;

	map<int, vector<int> > proposedFacets;
	map<int, set<int> > proposedFacet_GIDs;

	vector<int> facetsCondition1, facetsCondition2;

	// Loop through the matching facets and split them
	for (uint i = 0; i < matchingFacets.size(); ++i) {
		if (debug(edge))
			cout << "Looping through matching facet: " << matchingFacets[i] << endl;

		if (facetSplitStatus.count(matchingFacets[i])) {
			split = false;
			break;
		}

		// Get a copy of this facet.
		vector<int> currentFacet1 = globalFacets.at(matchingFacets[i]);
		vector<int> currentFacet2 = globalFacets.at(matchingFacets[i]);

		if (debug(edge)){
			cout << "Obtained matching facet: " << matchingFacets[i] << endl;
			cout << "Original Facet FID" << matchingFacets[i] << ": NID"
					<< currentFacet1[0] << ", NID" << currentFacet1[1]
					<< ", NID" << currentFacet1[2] << endl;
		}

		// Build the proposed facets by replacing nodes
		replace(currentFacet1.begin(), currentFacet1.end(), nid1, addedNID);
		replace(currentFacet2.begin(), currentFacet2.end(), nid2, addedNID);

		// Place the new facets in the proposed facets
		proposedFacets[matchingFacets[i]] = currentFacet1;
		proposedFacet_GIDs[matchingFacets[i]] = facet2grain.at(matchingFacets[i]);
		facetsCondition1.push_back(matchingFacets[i]);

		proposedFacets[max_facet_id + i] = currentFacet2;
		proposedFacet_GIDs[max_facet_id + i] = facet2grain.at(matchingFacets[i]);
		facetsCondition2.push_back(max_facet_id + i);
		// Dont increment max_facet_id until the propsoed facets are finalzied

		// Calculate the normal for the original facet
		if (debug(edge))
			cout << "Grabbing normal for original facet: " << matchingFacets[i]
					<< endl;
		vector<double> currentNormal = { 0, 0, 0 };
		calculateFacetNormal(globalFacets.at(matchingFacets[i]), currentNormal);// Dont use currentFacet becasue it's been replaced
		
		// Calculate the normal for the proposed facets
		if (debug(edge))
			cout << "Grabbing Normal for Facet 1" << endl;
		vector<double> proposedNormal1 = { 0, 0, 0 };
		calculateFacetNormal(proposedFacets.at(matchingFacets[i]),proposedNormal1);
		if (debug(edge))
			cout << "Grabbing Normal for Facet 2" << endl;
		vector<double> proposedNormal2 = { 0, 0, 0 };
		calculateFacetNormal(proposedFacets.at(max_facet_id), proposedNormal2);

		if (debug(edge))
			cout << "Grabbing Angle Between Normals" << endl;
		double angleBtweenN1 = calculateAngleBetweenNormals(currentNormal,
				proposedNormal1);
		double angleBtweenN2 = calculateAngleBetweenNormals(currentNormal,
				proposedNormal2);

		if (isnan(angleBtweenN1)) {
			cout << "*** u: " << currentNormal[0] << " " << currentNormal[1]
					<< " " << currentNormal[2] << endl;
			cout << "*** v1: " << proposedNormal1[0] << " "
					<< proposedNormal1[1] << " " << proposedNormal1[2] << endl;
			vector<double> u = currentNormal;
			vector<double> v = proposedNormal1;
			cout << "u.v = " << u[0] * v[0] + u[1] * v[1] + u[2] * v[2] << endl;
			split = false;
			break;
		}

		if (isnan(angleBtweenN2)) {
			cout << "*** u: " << currentNormal[0] << " " << currentNormal[1]
					<< " " << currentNormal[2] << endl;
			cout << "*** v2: " << proposedNormal2[0] << " "
					<< proposedNormal2[1] << " " << proposedNormal2[2] << endl;
			vector<double> u = currentNormal;
			vector<double> v = proposedNormal2;
			cout << "u.v = " << u[0] * v[0] + u[1] * v[1] + u[2] * v[2] << endl;
			split = false;
			break;
		}

		if ((angleBtweenN1 > normalTol) && (angleBtweenN2 > normalTol)) {
			if (debug(edge))
				cout << "Failed due to angle between normals: " << angleBtweenN1
						<< " " << angleBtweenN2 << endl;
			split = false;
			break;
		}
	}

	// If the splitting function failed, undo everything 
	if (!split) {
		// Erase the node I temporarily created
		nodes.erase(kd_node_id);
		if (debug(edge))
			cout << "Split failed due to mismatched normal or chord." << endl;
		return -1;
	} else {
		if (debug(edge))
			cout << "Collapsed one." << endl;
	}

	// We will skip facet swapping in refinement for now... 

	// Update mesh information
	// Copy proposed facets over
	set<int> modifiedNodes;
	set<int> auxNodes;
	map<int, vector<int>>::iterator it;
	for (it = proposedFacets.begin(); it != proposedFacets.end(); it++) {
		if (debug(edge))
			cout << "Writing facet FID" << it->first << ": NID" << it->second[0]
					<< " NID" << it->second[1] << " NID" << it->second[2]
					<< endl;
		globalFacets[it->first] = it->second;
		//		modifiedFacets.insert(it->first);
		modifiedNodes.insert(it->second.begin(), it->second.end());
		addToFacet2Grain(it->first, proposedFacet_GIDs.at(it->first));
	}

	for (auto n : modifiedNodes) {
		if (n != nid1 || n != nid2 || n != kd_node_id) {
			auxNodes.insert(n);
		}
	}

	// Add new facets to node2facet
	vector<int> facets4n2f;

	for (uint i = 0; i < matchingFacets.size(); ++i) {
		facets4n2f.push_back(max_facet_id + i);
		facets4n2f.push_back(matchingFacets[i]);
	}
	node2facet[kd_node_id] = facets4n2f;

	for (uint d = 0; d < facetsCondition1.size(); d++) {
		vector<int> nidn1f = node2facet.at(nid1);
		nidn1f.erase(
				std::remove(nidn1f.begin(), nidn1f.end(), facetsCondition1[d]),
				nidn1f.end());
		node2facet[nid1] = nidn1f;
	}

	for (uint d = 0; d < facetsCondition2.size(); d++) {
		vector<int> nidn1f = node2facet.at(nid1);
		nidn1f.push_back(facetsCondition2[d]);
		node2facet[nid1] = nidn1f;
	}

	for (auto n : auxNodes) {
		vector<int> nidnnf = node2facet.at(n);
		for (auto f : facets4n2f) {
			vector<int> fsnid = globalFacets.at(f);
			if (count(fsnid.begin(), fsnid.end(), n) > 0) {
				if (count(nidnnf.begin(), nidnnf.end(), f) == 0)
					nidnnf.push_back(f);
			} else {
				nidnnf.erase(std::remove(nidnnf.begin(), nidnnf.end(), f),
						nidnnf.end());
			}
		}
		node2facet[n] = nidnnf;
	}

	// Update Node Neighbors
	// If successful, nid1 will no longer be neighbors with nid2
	vector<int> nn1f = nodeNeighbors.at(nid1);
	nn1f.erase(std::remove(nn1f.begin(), nn1f.end(), nid2), nn1f.end());
	nn1f.push_back(kd_node_id);
	nodeNeighbors[nid1] = nn1f;
	vector<int> nn2f = nodeNeighbors.at(nid2);
	nn2f.erase(std::remove(nn2f.begin(), nn2f.end(), nid1), nn2f.end());
	nn2f.push_back(kd_node_id);
	nodeNeighbors[nid2] = nn2f;

	// Add the new node into nodeNeighbors
	vector<int> nn_newNode(modifiedNodes.begin(), modifiedNodes.end());
	nodeNeighbors[kd_node_id] = nn_newNode;

	for (auto n : auxNodes) {
		vector<int> nnnf = nodeNeighbors.at(n);
		nnnf.push_back(kd_node_id);
		nodeNeighbors[n] = nnnf;
	}

	// If we are here, we have successfully completed a refinement operation
	if (debug(edge))
		cout << "****************** END of GlobalMesh coarsenEdge function."
				<< endl;

	// Iterate relavant parameters:
	kd_node_id++;
	for (uint i = 0; i < matchingFacets.size(); ++i) {
		facetSplitStatus.insert(max_facet_id);
		max_facet_id++;
		facetSplitStatus.insert(matchingFacets[i]);
	}

	return 0;
}

GlobalMesh::~GlobalMesh() {
	// TODO Auto-generated destructor stub
}


void GlobalMesh::GetSharingEdge(vector<int> &facet_pair, vector<int> &SE) {
	vector<int> f1 = { facet_pair[0], facet_pair[1], facet_pair[2] };
	vector<int> f2 = { facet_pair[3], facet_pair[4], facet_pair[5] };
	for (auto n : f1) {
		if (std::find(f2.begin(), f2.end(), n) != f2.end()) {
			SE.push_back(n);
		}
	}

	// 2 facets should only share one edge, which has 2 nodes
	if (SE.size() != 2) {
		cout << "ERROR: Did not find the right number of sharing edge for facets:" << endl;
		cout << "F1: " << f1[0] << "," << f1[1] << "," << f1[2] << endl;
		cout << "F2: " << f2[0] << "," << f2[1] << "," << f2[2] << endl;
		exit(1); // This is fatal, just end the program
	}
}

int GlobalMesh::ContainsEdge(vector<int> &SE, vector<int> &facet_pair, vector<int> &target) {
	// Unpack the facet pair into two facets
	vector<int> f1 = { facet_pair[0], facet_pair[1], facet_pair[2] };
	vector<int> f2 = { facet_pair[3], facet_pair[4], facet_pair[5] };

	if (std::find(f2.begin(), f2.end(), SE[0]) != f2.end() && std::find(f2.begin(), f2.end(), SE[0]) != f2.end()) {
		// If facet 2 contains the edge we provide, this is the target facet we are looking for
		target = vector<int>(f2.begin(), f2.end());
	} else {
		target = vector<int>(f1.begin(), f1.end());
	}

	// Also, return the node in the target facet that's not on the edge we provide
	int tn = -1;
	for (auto n : target) {
		if (std::find(SE.begin(), SE.end(), n) == SE.end()) {
			tn = n;
			break;
		}
	}
	return tn;
}

void GlobalMesh::GetOpposite(vector<int> &facet_pair, vector<int> &fp1,
		int ref_node, vector<int> &target_facet) {
	// Unpack the facet pair into two facets
	vector<int> f1 = { facet_pair[0], facet_pair[1], facet_pair[2] };
	vector<int> f2 = { facet_pair[3], facet_pair[4], facet_pair[5] };

	// Get the edge connected to ref_node and in the positive convention, following the orientation convention of facet 1
	vector<int> SE; // Sharing edge, signed
	if (fp1[0] == ref_node)
		SE = {fp1[1] , fp1[2]};
	else if ( fp1[1] == ref_node )
		SE = {fp1[2] , fp1[0]};
	else 
		SE = {fp1[0] , fp1[1]};

	// Check if SE (signed) exists in facet 1. If so, then facet 1 is on the SAME side of new_pair1, return facet 2 instead
	int loc = -1;
	for (int idx = 0; idx < 3; idx++) {
		if (std::find(SE.begin(), SE.end(), f1[idx]) == SE.end()) {
			loc = idx;
			break;
		}
	}
	bool return_f2 = false;
	if (loc == 0) {
		if (f1[1] == SE[0] && f1[2] == SE[1]) {
			return_f2 = true;
		}
	} else if (loc == 1) {
		if (f1[2] == SE[0] && f1[0] == SE[1]) {
			return_f2 = true;
		}
	} else {
		if (f1[0] == SE[0] && f1[1] == SE[1]) {
			return_f2 = true;
		}
	}

	if (return_f2) {
		target_facet = vector<int>(f2.begin(), f2.end());
	} else {
		target_facet = vector<int>(f1.begin(), f1.end());
	}
}

int GlobalMesh::Overlap(vector<int> &FP1, vector<int> &FP2) {
	// Check if the two facets share the same exact nodes
	for (int i = 0; i < 2; i++) {
		vector<int> cf = { FP1[2 * i], FP1[2 * i + 1], FP1[2 * i + 2] };
		for (int j = 0; j < 2; j++) {
			vector<int> of = { FP2[2 * j], FP2[2 * j + 1], FP2[2 * j + 2] };
			if (cf[0] == of[0] && cf[1] == of[1] && cf[2] == of[2]) {
				return 1;
			}
		}
	}
	return 0;
}

void GlobalMesh::UpdateFacetPair(vector<vector<int> > &facet_pair, vector<int> &fp_new,int idx) {
	vector<int> SE, new_pair1, new_pair2;
	GetSharingEdge(facet_pair[idx], SE);
	int ref_node = ContainsEdge(SE, fp_new, new_pair1);
	GetOpposite(facet_pair[idx], new_pair1, ref_node, new_pair2);
	vector<int> tmp = { new_pair1[0], new_pair1[1], new_pair1[2], new_pair2[0],
			new_pair2[1], new_pair2[2] };
	for (int i = 0; i < 6; i++) {
		facet_pair[idx][i] = tmp[i];
	}
}

int GlobalMesh::CheckNormal(vector<int> &auxFacets, vector<double> &F1,
		int nid1, int nid2, double normalTol) {
	for (uint i = 0; i < auxFacets.size(); ++i) {
		// Original normal
		vector<double> F = { 0, 0, 0 };
		vector<int> cf = globalFacets.at(auxFacets[i]);
		calculateFacetNormal(cf, F);

		// Compute averaged normal of all surrounding facets
		for (int ii = 0; ii < 3; ii++) {
			F1[ii] += F[ii] / (float) auxFacets.size();
		}

		// Modified facet by using the temporary node (node 0)
		replace(cf.begin(), cf.end(), nid1, 0);
		replace(cf.begin(), cf.end(), nid2, 0);

		// Modified normal
		vector<double> F2 = { 0, 0, 0 };
		calculateFacetNormal(cf, F2);

		double angleBtweenN = calculateAngleBetweenNormals(F, F2);

		if (isnan(angleBtweenN)) {
			#if (DEBUG > 0 )
				cout << "Failed due to one of the normals is NAN" << endl;
			#endif
			return -1;
		}

		if (angleBtweenN > normalTol) {
			#if (DEBUG > 0 )
				cout << "Failed due to angle between matching facets is too large: " << angleBtweenN << endl;
			#endif
			return -2;
		}
	}
	return 0;
}

int GlobalMesh::ExternalLoop(vector<int> &auxFacets, vector<int> &External_Loop,
		int nid1, int nid2, int MODE) {
	// To get all the external edges that form the boundary of the closed polygon, we go through all the aux facets
	set<int> Unique_nodes;
	vector<vector<int> > temp_set;
	for (auto cfid : auxFacets) {
		vector<int> fn = globalFacets.at(cfid); // Get all nodes on this facet
		vector<int> curr_vector;

		// Order the current external edge using default positive convention
		if (fn[0] == nid1 || fn[0] == nid2) {
			curr_vector = {fn[1] , fn[2]};
		}
		else if ( fn[1] == nid1 || fn[1] == nid2 ) {
			curr_vector = {fn[2] , fn[0]};
		}
		else {
			curr_vector = {fn[0] , fn[1]};
		}

		// Store this external edge
		temp_set.push_back(curr_vector);
	}
	int num_edges = temp_set.size();

	// Fill the external edge loop, connect the external edges in a tip-to-tail manner
	External_Loop = temp_set[0]; // Put the first edge in to get started
	Unique_nodes.insert(temp_set[0][0]);
	Unique_nodes.insert(temp_set[0][1]);
	int max_try = 100; // Iteratively try max_try times before aborting
	int count = 0;

	// If MODE = 0, we are coarsening a regular edge. In this case, the external edges form a closed polygon
	// The number of nodes in the loop = the number of edges in the loop

	// If MODE = 1, we are coarsening a boundary edge. In this case, the non-boundary external edges do not form a closed polygon
	// The number of nodes in the loop = the number of edges in the loop + 1

	while (External_Loop.size() != num_edges + MODE && count < max_try) {
		count += 1;
		// Loop through all remaining external edges
		for (int i = 1; i < num_edges; i++) {
			vector<int> cv = temp_set[i]; // Get the current external edge

			// If the end of the current edge = begining of the queue, put this edge at the begining of the queue 
			if (cv[1] == External_Loop[0]) {
				External_Loop.insert(External_Loop.begin(), cv[0]);
				Unique_nodes.insert(cv[0]);
				break;
			}

			// If the begining of the current edge = end of the queue, put this edge at the end of the queue 
			if (cv[0] == External_Loop.back()) {
				External_Loop.push_back(cv[1]);
				Unique_nodes.insert(cv[1]);
				break;
			}
		}
	}

	// Sanity check
	// Check 1: See if we indeed found all the nodes on the external edge loop
	if (External_Loop.size() != num_edges + MODE) {
		#if (COARSENDEBUG > 0)
		cout << ">>> WARNING: Didn't find all external edges, skipping" << endl;
		#endif
		num_edges = -2;
	}
	// Check 2: See if all nodes in the external loop are unique
	// If they are not, this indicates that the external loop has a self-intersection.
	// This case is too complex and is not handled
	else if (External_Loop.size() != Unique_nodes.size()) {
		#if (COARSENDEBUG > 0)
		cout << ">>> WARNING: At least one self-intersection on the external loop, skipping" << endl;
		#endif
		num_edges = -2;
	}
	return num_edges;
}

int GlobalMesh::FindConcavity(vector<int> &External_Loop, vector<double> &F1,
		int num_edges, vector<int> &concave_loc) {
	// This implementation is only an approximate. The 'right' way to do this check is to first project the nodes on the 3D polygon onto a
	// nominal plane, then check the angles on the projected, planar polaygon. This implementation does not project the 3D (possibly non-planar)
	// polygon into a planar one.

	int local_concave = 0;
	for (int i = 0; i < num_edges; i++) {
		// Make a temporary local facet using the 3 nodes
		int idx2 = i + 1;
		int idx3 = i + 2;
		if (idx3 >= num_edges)
			idx3 -= num_edges;
		if (idx2 >= num_edges)
			idx2 -= num_edges;
		// cout << "Checking local facet with nodes: " << External_Loop[i] << " , " << External_Loop[idx2] << " , " << External_Loop[idx3] << endl;

		// Get normal of this temp facet
		vector<int> test_facet = { External_Loop[i], External_Loop[idx2], External_Loop[idx3] };
		vector<double> TN = { 0, 0, 0 };
		calculateFacetNormal(test_facet, TN);

		// Compare this normal to the reference normal. If this is a local concavity, the normal of the test facet points into the opposite direction
		// of the reference normal (which is the averaged normal of all aux facets)
		if (F1[0] * TN[0] + F1[1] * TN[1] + F1[2] * TN[2] < 0.) {
			local_concave += 1;
			concave_loc.push_back(idx2); // Also store the location of the local concavity
		}
	}	
	return local_concave;
}

int GlobalMesh::SwapFacets(vector<vector<int> > &facet_pair,
		vector<vector<int> > &initial_facet_pair, vector<int> &concave_loc,
		double minAngle, double maxAngle, double normalTol,
		vector<vector<int> > &final_facet_pair) {
	bool Initial_Good = true;
	int total_swap = 0;

	#if (SWAPDEBUG > 0)
		cout << "Total of " << facet_pair.size() << " facet pairs" << endl;
	#endif
	
	// Only swap the facet pairs when there are more than 1 pair, otherwise, just use the initial guess
	if ( facet_pair.size() > 1 ){
		for (int i = 0; i < facet_pair.size(); i++) {
			// cout << "Checking facet pair " << i << endl;

			// Get current facet pair
			vector<int> fp = facet_pair[i];

			// Check if the current sharing edge is connected to a concavity
			// Do not swap this edge if it is directly connected to a local concavity
			// Swapping this will in general, introduce overlapping facets
			vector<int> E;
			GetSharingEdge(facet_pair[i], E);
			if (std::find(concave_loc.begin(), concave_loc.end(), E[0])!= concave_loc.end()
					|| std::find(concave_loc.begin(), concave_loc.end(), E[1])!= concave_loc.end()) {
				#if (SWAPDEBUG > 0)
					cout << ">>> WARNING: The sharing edge of current facet pair is connected to a local concavity, do not swap this edge" << endl;
				#endif
				continue;
			}

			// Get initial angle metrics
			// Initial facet 1
			vector<double> if1_a = { 0, 0, 0 };
			vector<int> if1 = { fp[0], fp[1], fp[2] };
			calculateFacetAngles(if1, if1_a);
			// Initial facet 2
			vector<double> if2_a = { 0, 0, 0 };
			vector<int> if2 = { fp[3], fp[4], fp[5] };
			calculateFacetAngles(if2, if2_a);

			double initial_min = std::min(if1_a[0], if2_a[0]);
			double initial_max = std::max(if1_a[2], if2_a[2]);

			// Check how good is the initial guess triangulation
			// The initial guess is disqualified, when the min face angle is smaller than tolerance, or when the max face angle is larger than tolerance
			// Otherwise, the initial guess is qualified to be used. This does not mean that the initial guess will ALWAYS be used.
			if (initial_min < minAngle || initial_max > maxAngle) {
				Initial_Good = false;
			}

			// Swap the facet pair by swapping the diagonal edge
			vector<int> fp_new = { fp[0], fp[1], fp[5], fp[5], fp[1], fp[2] };

			// Get swapped angle metrics
			// Swapped facet 1
			vector<double> sf1_a = { 0, 0, 0 };
			vector<int> sf1 = { fp_new[0], fp_new[1], fp_new[2] };
			calculateFacetAngles(sf1, sf1_a);
			// Swapped facet 2
			vector<double> sf2_a = { 0, 0, 0 };
			vector<int> sf2 = { fp_new[3], fp_new[4], fp_new[5] };
			calculateFacetAngles(sf2, sf2_a);

			double swapped_min = std::min(sf1_a[0], sf2_a[0]);
			double swapped_max = std::max(sf1_a[2], sf2_a[2]);
			bool can_swap = false;

			// Check if the swapped facets have better metrics than the original ones
			if (swapped_min > initial_min && swapped_min > minAngle
					&& swapped_max < maxAngle) {
				// We passed all angle metric checks, now let's also check the normals
				// Initial normals
				vector<double> IN1 = { 0, 0, 0 };
				calculateFacetNormal(if1, IN1);
				vector<double> IN2 = { 0, 0, 0 };
				calculateFacetNormal(if2, IN2);

				// New normals
				vector<double> SN1 = { 0, 0, 0 };
				calculateFacetNormal(sf1, SN1);
				vector<double> SN2 = { 0, 0, 0 };
				calculateFacetNormal(sf2, SN2);

				if (calculateAngleBetweenNormals(IN1, SN1) < normalTol
						&& calculateAngleBetweenNormals(IN2, SN2) < normalTol) {
					// Great, we passed the normal check!
					can_swap = true;
					total_swap += 1;

					// If we get to here, we have succcessfully passed all checks. We will then use the swapped facets
					#if (SWAPDEBUG > 0)
						cout << "Passed all checks, swapping!" << endl;
					#endif

					// Update mesh information
					// Update neighboring facet pairs (i+1 and i-1)
					// Backward update
					if (i != 0) {
						UpdateFacetPair(facet_pair, fp_new, i - 1);
					}
					// Forward update
					if (i != facet_pair.size() - 1) {
						UpdateFacetPair(facet_pair, fp_new, i + 1);
					}

					// Check if the current old facet pair is connected to n-2th or n+2th facet pair				
					if (i >= 2) {
						int test1 = Overlap(facet_pair[i], facet_pair[i - 2]);
						if (test1 == 1) {
							UpdateFacetPair(facet_pair, fp_new, i - 2);
						}
					}

					if (i < facet_pair.size() - 2) {
						int test1 = Overlap(facet_pair[i], facet_pair[i + 2]);
						if (test1 == 1) {
							UpdateFacetPair(facet_pair, fp_new, i + 2);
						}
					}

					// Copy current swapped facet
					facet_pair[i] = vector<int>(fp_new.begin(), fp_new.end());
					#if (SWAPDEBUG > 0)
						cout << "Swapped one" << endl;
					#endif
				}
			}

			// If the original facets are qualified and are better than the swapped ones, use them instead
			if (can_swap == 0 && Initial_Good == 1) {
				#if (SWAPDEBUG > 0)
					cout << "Local facet pair is already Delaunay" << endl;
				#endif
			}
		}
	}else{
		#if (SWAPDEBUG > 0)
			cout << "Only one facet pair in here! Not swapping anything!" << endl;
		#endif
	}

	#if (SWAPDEBUG > 0)
		cout << "All swapping completed!" << endl;
	#endif

	// In the case that we didn't make any successful swaps, check if the old ones can be used
	if (total_swap == 0) {
		// There are no successful swaps, check if we can use the initial guess
		if (Initial_Good == 1) {
			#if (SWAPDEBUG > 0)
				cout << "No successful swaps but our initial guess passed all checks, use this instead" << endl;
			#endif
			final_facet_pair = vector<vector<int> >(initial_facet_pair.begin(),
					initial_facet_pair.end());
		} else {
			#if (SWAPDEBUG > 0)
				cout << "We cannot swap and cannot use the initial guess, too bad, just leave" << endl;
			#endif
			return -1;
		}
	} else {
		// We made some successful swaps, use this one
		#if (SWAPDEBUG > 0)
			cout << "Using swapped facets" << endl;
		#endif
		final_facet_pair = vector<vector<int> >(facet_pair.begin(),
				facet_pair.end());
	}
	return 0;
}

void GlobalMesh::UpdateMesh(vector<int> &auxFacets,
		map<int, vector<int> > &proposedFacets, int nid1, int nid2) {
	// Step 7.1: node2facet and nodeNeighbors
	// Remove old facets in node2facet, remove old nodes in nodeNeighbors
	vector<int> nn2f;
	for (auto fid : auxFacets) {
		vector<int> cf = globalFacets.at(fid);
		for (int i = 0; i < 3; i++) {
			nn2f = node2facet.at(cf[i]);
			nn2f.erase(std::remove(nn2f.begin(), nn2f.end(), fid), nn2f.end());
			node2facet[cf[i]] = nn2f;

			if (cf[i] != nid1 && cf[i] != nid2) {
				nn2f = nodeNeighbors.at(cf[i]);
				// Remove all internal connections to nid1 and nid2
				nn2f.erase(std::remove(nn2f.begin(), nn2f.end(), nid1),
						nn2f.end());
				nn2f.erase(std::remove(nn2f.begin(), nn2f.end(), nid2),
						nn2f.end());
				nodeNeighbors[cf[i]] = nn2f;
			}
		}
	}
	// Remove nid2
	nodes.erase(nid2);
	node2facet.erase(nid2);
	nodeNeighbors.erase(nid2);
	// Clear nid1
	node2facet[ nid1 ] = {};
	nodeNeighbors[ nid1 ] = {};

	// Step 7.2: globalFacets
	// Remove old facets in globalFacets
	for (auto fid : auxFacets) {
		globalFacets.erase(fid);
	}

	// Step 7.3: proposedFacets
	map<int, vector<int>>::iterator it;
	for (it = proposedFacets.begin(); it != proposedFacets.end(); it++) {
		int cfid = it->first;
		vector<int> cf = it->second;

		// Add new facets
		globalFacets[cfid] = cf;

		// Update node2facet
		for (int i = 0; i < 3; i++) {
			nn2f = node2facet.at(cf[i]);
			if (std::find(nn2f.begin(), nn2f.end(), cfid) == nn2f.end()) { // Add the current facet if it does not exist
				nn2f.push_back(cfid);
			}
			node2facet[cf[i]] = nn2f;

			// Update nodeNeighbors
			nn2f = nodeNeighbors.at(cf[i]);
			for (int j = 0; j < 3; j++) {
				// Loop through all other nodes on the same facet, add to list if they don't already exist
				if ((i != j)
						&& std::find(nn2f.begin(), nn2f.end(), cf[j])
								== nn2f.end()) {
					nn2f.push_back(cf[j]);
				}
			}
			nodeNeighbors[cf[i]] = nn2f;
		}
	}

	// Step 7.4: Record all nodes that are removed 
	obliteratedNodes.insert(nid2);
	#if (COARSENDEBUG > 0)
		cout << "DONE updating mesh" << endl;
	#endif
}

// Collapse edges, for regular edges
int GlobalMesh::collapseEdge(vector<int> edge, int nid1, int nid2,
		double normalTol, double chordTol, double minAngle, double maxAngle) {
	#if (COARSENDEBUG > 0)
		cout << "****************** BEGIN of GlobalMesh coarsenEdge function."<< endl;
	#endif
	
	vector<double> &xyznid1 = nodes.at(nid1);
	vector<double> &xyznid2 = nodes.at(nid2);
	vector<double> mid_pt = { 0.5 * (xyznid1[0] + xyznid2[0]), 0.5
			* (xyznid1[1] + xyznid2[1]), 0.5 * (xyznid1[2] + xyznid2[2]) };
	// Make a temp node at the mid-point of the target edge with node ID of 0
	nodes[0] = mid_pt;
	
	// Store old nodal coords of nid1
	vector<double> old_nc = { xyznid1[0], xyznid1[1], xyznid1[2] };

	#if (COARSENDEBUG > 0)
		cout << "NID1 " << nid1 << " Coords: " << xyznid1[0] << " "
				<< xyznid1[1] << " " << xyznid1[2] << endl;
		cout << "NID2 " << nid2 << " Coords: " << xyznid2[0] << " "
				<< xyznid2[1] << " " << xyznid2[2] << endl;
	#endif

	// First, get the facets connected to the vertices. Luckily, we have a node2facet database to use
	vector<int> &facetlist1 = node2facet.at(nid1);
	vector<int> &facetlist2 = node2facet.at(nid2);

	// Step 1: Check if dragging nodes will significantly influence the geometry
	vector<int> matchingFacets, auxFacets;
	for (uint i = 0; i < facetlist2.size(); ++i) {
		if (std::find(facetlist1.begin(), facetlist1.end(), facetlist2[i])
				!= facetlist1.end()) {
			matchingFacets.push_back(facetlist2[i]);
		} else {
			auxFacets.push_back(facetlist2[i]);
		}
	}
	for (uint i = 0; i < facetlist1.size(); ++i) {
		if (std::find(facetlist2.begin(), facetlist2.end(), facetlist1[i])
				== facetlist2.end()) {
			auxFacets.push_back(facetlist1[i]);
		}
	}
	// Perform some initial sanity check on the matching facets
	if (matchingFacets.size() == 0) {
		if (debug(edge))
			cout
					<< "*** ERROR: Found zero matching facets in collapsingEdge function for nids "
					<< nid1 << ", " << nid2 << "." << endl;
		return -1;
	} else if (matchingFacets.size() == 1) {
		if (debug(edge))
			cout
					<< "*** ERROR: Found only one matching facet in collapsingEdge function for nids "
					<< nid1 << ", " << nid2
					<< ". This should not happen for manifold geometries."
					<< endl;
		return -1;
	} else if (matchingFacets.size() > 2) {
		if (debug(edge))
			cout << "*** DEBUG OUTPUT: Found " << matchingFacets.size()
					<< " matching facets in collapsingEdge function for nids "
					<< nid1 << ", " << nid2 << "." << endl;
		return -1;
	}

	if (debug(edge)) {
		cout << "All Facets Being Worked With: " << endl;
		cout << "Matched: " << matchingFacets[0] << " " << matchingFacets[1]
				<< endl;
		cout << "Aux: ";
		for (uint i = 0; i < auxFacets.size(); ++i) {
			cout << auxFacets[i] << " ";
		}
		cout << endl;
	}

	// Now that we have all aux facets, check the normal
	vector<double> F1 = { 0, 0, 0 }; // Reference normal , used later
	int pass = CheckNormal(auxFacets, F1, nid1, nid2, normalTol);
	if (pass != 0) {
		return -1;
	}
	// If we get to here, we are good to collapse this edge
	#if (COARSENDEBUG > 0)
		cout << "Normal check passed! OK to collapse edge!" << endl;
	#endif


	// From past experence, the facet ordering conventions in the input mesh may not be consistent,
	// especially near sharp edges. For the code to work properly, we locally unify the facet ordering 
	// conventions of all aux facets to the default positive convention.  
	// Step 1.5: Uinify node ordering conventions
	vector<int> AF = auxFacets;
	AF.push_back(matchingFacets[1]);
	int REF = matchingFacets[0];
	int ccc = 0;
	while (AF.size() != 0 && ccc < 100) {
		ccc += 1;
		vector<int> rf = globalFacets.at(REF);
		for (auto f : AF) {
			vector<int> cf = globalFacets.at(f);
			int t1 = std::find(rf.begin(), rf.end(), cf[0]) != rf.end();
			int t2 = std::find(rf.begin(), rf.end(), cf[1]) != rf.end();
			int t3 = std::find(rf.begin(), rf.end(), cf[2]) != rf.end();
			if (t1 + t2 + t3 == 2) { // These two share an edge
				UnifyConvention(REF, f);
				REF = f;
				AF.erase(std::remove(AF.begin(), AF.end(), f), AF.end());
				break;
			}
		}
	}

	// Check if we indeed unified the convention for ALL aux facets
	if (AF.size() != 0) {
		#if (COARSENDEBUG > 0)
		cout << ">>> WARNING: Facet convention unification failed, skipping this edge" << endl;
		#endif
		return -1;
	}


	// Step 2: Get nodes on the external edge loop, organize them in CCW order
	vector<int> External_Loop;
	int num_edges = ExternalLoop(auxFacets, External_Loop, nid1, nid2, 0);

	// Abort if weird cases are encountered
	if (num_edges == -1) {
		exit(1);
	} else if (num_edges == -2) {
		return -1;
	}

	// Check average edge length of the initial guess
	double avg_edge_length = 0.;
	vector<double> EL_center = { 0., 0., 0. };
	double LEN = (double) External_Loop.size();
	for (auto nnn : External_Loop) {
		vector<double> p = nodes.at(nnn);
		EL_center[0] += p[0] / LEN;
		EL_center[1] += p[1] / LEN;
		EL_center[2] += p[2] / LEN;
		avg_edge_length += calculateDistanceBetweenTwoDoubleVectors(p, mid_pt);
	}
	double ph;
	// Use the sizing function to decide if this needs to be coarsened or not
	if (sizingFunction(avg_edge_length / (float) num_edges, mid_pt, ph, 0) != 1) {
		return -1;
	}

	if (debug(edge)) {
		cout << "External edge loop: " << endl;
		for (auto q : External_Loop) {
			vector<double> c = nodes.at(q);
			cout << "Node " << q << ": " << c[0] << "," << c[1] << "," << c[2] << endl;
		}
	}

	// Step 3: Check for local concave polys
	vector<int> concave_loc;
	int local_concave = FindConcavity(External_Loop, F1, num_edges, concave_loc);


	// Step 4: Reorder External_Loop to start with the first local concavity
	if (local_concave > 0 && concave_loc[0] != 0) {
		#if (COARSENDEBUG > 0)
			cout << "Reordering External_Loop starting at position " << concave_loc[0] << endl;
		#endif
		vector<int> t = vector<int>(External_Loop.begin() + concave_loc[0],
				External_Loop.end());
		for (int i = 0; i < concave_loc[0]; i++) {
			t.push_back(External_Loop[i]);
		}
		External_Loop = vector<int>(t.begin(), t.end());
		if (debug(edge)) {
			cout << "Reordered external edge loop: ";
			for (auto q : External_Loop) {
				cout << q << " , ";
			}
			cout << endl;
		}
	}

	// Step 5: Form all initial surface pairs
	// Step 5.1: Move node 1 to the mid-point of the edge
	nodes[nid1] = mid_pt;

	// Step 5.2: Make initial triangulations
	vector<vector<int> > facet_pair, initial_facet_pair;
	for (int i = 1; i < num_edges; i++) {
		vector<int> curr_pair;
		if (i == num_edges - 1) {
			curr_pair = {nid1 , External_Loop[i - 1] , External_Loop[i] , nid1 , External_Loop[i] , External_Loop[0]};
		}
		else {
			curr_pair = {nid1 , External_Loop[i - 1] , External_Loop[i] , nid1 , External_Loop[i] , External_Loop[i + 1]};
		}
		facet_pair.push_back(curr_pair);
		initial_facet_pair.push_back(curr_pair); // Keep an unmodified copy
	}

	// Step 6: Attempt to swap those surface pairs
	vector<vector<int> > final_facet_pair;
	int swapped = SwapFacets(facet_pair, initial_facet_pair, concave_loc,
			minAngle, maxAngle, normalTol, final_facet_pair);
	if (swapped == -1) {
		nodes[nid1] = old_nc;
		return -1;
	}

	// Step 7: Update connectivity
	// Collect all old facets
	auxFacets.push_back(matchingFacets[0]);
	auxFacets.push_back(matchingFacets[1]);

	// Just use the old facet IDs on the new facets, we should have enough
	if (auxFacets.size() > final_facet_pair.size() + 1) {
		// Build proposedFacets
		map<int, vector<int> > proposedFacets;
		vector<vector<int> > tmpFacet;
		for (int i = 0; i < final_facet_pair.size(); i++) {
			vector<int> tf1 = { final_facet_pair[i][0], final_facet_pair[i][1],
					final_facet_pair[i][2] };
			vector<int> tf2 = { final_facet_pair[i][3], final_facet_pair[i][4],
					final_facet_pair[i][5] };
			// Check if tf1 is in tmpFacet
			if (std::find(tmpFacet.begin(), tmpFacet.end(), tf1)
					== tmpFacet.end()) {
				tmpFacet.push_back(tf1);
			}
			// Check if tf2 is in tmpFacet
			if (std::find(tmpFacet.begin(), tmpFacet.end(), tf2)
					== tmpFacet.end()) {
				tmpFacet.push_back(tf2);
			}
		}
		#if (COARSENDEBUG > 0)
			cout << "Found " << tmpFacet.size() << " unique facets" << endl;
		#endif

		// Sanity check, see if we created duplicate facets
		if (CheckOverlap(tmpFacet, External_Loop, auxFacets)) {
			#if (COARSENDEBUG > 0)
			cout << ">>> WARNING: Swapping algorithm created overlaping facets! Skipping!" << endl;
			#endif
			// If we fail, we need to restore location of node 1 to its original location
			nodes[nid1] = old_nc;
			return -1;
		}

		// Sanity check, see if the swapping algorithm produces unexpected results
		// The number of UNIQUE facets after swapping, should equal to the number of facets in the initial guess (although the facets are different)
		if (tmpFacet.size() != initial_facet_pair.size() + 1) {
			#if (COARSENDEBUG > 0)
			cout << ">>> WARNING: Initially there are " << initial_facet_pair.size() + 1
					<< " unique facets, now we have " << tmpFacet.size() << endl;
			#endif
			// If the two numbers do not agree, this is most likely because we have a facet with a free edge (that partially overlaps with other facets)
			// This is potentially very bad, but we can attempt to fix it by removing the problematic facet
			RemoveFreeEdge(tmpFacet, External_Loop, initial_facet_pair.size());
		}

		for (int i = 0; i < tmpFacet.size(); i++) {
			proposedFacets[auxFacets[i]] = tmpFacet[i];
			#if (COARSENDEBUG > 0)
				cout << "New facet " << auxFacets[i] << " has nodes: " << tmpFacet[i][0] << "," << tmpFacet[i][1] << "," << tmpFacet[i][2] << endl;
			#endif
		}

		UpdateMesh(auxFacets, proposedFacets, nid1, nid2);
		for (auto q : External_Loop) {
			obliteratedNodes.insert(q);
		}

	} else {
		cout << "*** ERROR: Somehow we ran out of facet IDs...." << endl;
		cout << "Old: " << auxFacets.size() << " vs. New: "
				<< final_facet_pair.size() + 1 << endl;
		exit(1);
	}

	return 0;
}


// Collapse edges, for interface edges
int GlobalMesh::collapseEdge3(int MODE, vector<int> edge, int nid1, int nid2,
		double normalTol, double chordTol, double minAngle, double maxAngle) {
	#if (COARSENDEBUG > 0)
		cout << "****************** BEGIN of GlobalMesh coarsenEdge3 function."<< endl;
	#endif

	vector<double> xyznid1 = nodes.at(nid1);
	vector<double> xyznid2 = nodes.at(nid2);
	vector<double> mid_pt = { 0.5 * (xyznid1[0] + xyznid2[0]), 0.5
			* (xyznid1[1] + xyznid2[1]), 0.5 * (xyznid1[2] + xyznid2[2]) };
	// Make a temp node at the center of the target edge and assign node ID of 0
	nodes[0] = mid_pt;

	// Store old nodal coords of nid1
	vector<double> old_nc = { xyznid1[0], xyznid1[1], xyznid1[2] };

	// First, get the facets connected to the vertices. Luckily, we have a node2facet database to use
	vector<int> &facetlist1 = node2facet.at(nid1);
	vector<int> &facetlist2 = node2facet.at(nid2);

	// Step 1: Check if dragging nodes will significantly influence the geometry
	vector<int> matchingFacets, auxFacets;
	for (uint i = 0; i < facetlist2.size(); ++i) {
		if (std::find(facetlist1.begin(), facetlist1.end(), facetlist2[i])
				!= facetlist1.end()) {
			matchingFacets.push_back(facetlist2[i]);
		} else {
			auxFacets.push_back(facetlist2[i]);
		}
	}
	for (uint i = 0; i < facetlist1.size(); ++i) {
		if (std::find(facetlist2.begin(), facetlist2.end(), facetlist1[i])
				== facetlist2.end()) {
			auxFacets.push_back(facetlist1[i]);
		}
	}

	// Perform some initial sanity check on the matching facets
	if (matchingFacets.size() == 0) {
		if (debug(edge))
			cout << "*** ERROR: Found zero matching facets in collapsingEdge function for nids " << nid1 << ", " << nid2 << "." << endl;
		return -1;
	} else if (matchingFacets.size() == 1) {
		if (debug(edge))
			cout << "*** ERROR: Found only one matching facet in collapsingEdge function for nids " << nid1 << ", " << nid2
			<< ". This should not happen for manifold geometries." << endl;
		return -1;
	}
	if (matchingFacets.size() <= 2) {
		if (debug(edge))
			cout << "*** DEBUG OUTPUT: Found " << matchingFacets.size() << " matching facets in collapsingEdge function for nids "
			<< nid1 << ", " << nid2 << "." << endl;
		return -1;
	}
	#if (COARSENDEBUG > 0)
		cout << "All Facets Being Worked With: " << endl;
		cout << "Matched: ";
		for (uint i = 0; i < matchingFacets.size(); ++i) {
			cout << matchingFacets[i] << " ";
		}
		cout << endl;
		cout << "Aux: ";
		for (uint i = 0; i < auxFacets.size(); ++i) {
			cout << auxFacets[i] << " ";
		}
		cout << endl;
	#endif

	// Now that we have all aux facets, check the normal
	vector<double> F1 = { 0, 0, 0 }; // Reference normal , used later
	int pass = CheckNormal(auxFacets, F1, nid1, nid2, normalTol);
	if (pass != 0) {
		return -1;
	}

	// If we get to here, we are good to collapse this edge
	#if (COARSENDEBUG > 0)
		cout << "Normal check passed! OK to collapse edge!" << endl;
	#endif

	
	// Step 2: Get all nodes on the external edge loop, organize them in CCW order
	// To do this, we go through all the aux facets
	vector<vector<int> > Split_auxFacets;
	// Step 2.1: Find me all other boundary edges (should be 2 more, other than the edge we are currently working on) that contains nid1 or nid2
	vector<int> Boundary_Nodes;
	for (int ii = 0; ii < jEdges.size(); ii++) {
		vector<int> edge = jEdges[ii].n;
		if ( (edge[0] == nid1 && edge[1] != nid2) || (edge[0] == nid2 && edge[1] != nid1) ) {
			Boundary_Nodes.push_back(edge[1]);
		} else if ( (edge[1] == nid1 && edge[0] != nid2) || (edge[1] == nid2 && edge[0] != nid1) ) {
			Boundary_Nodes.push_back(edge[0]);
		}
	}

	// Check the number of boundary nodes, as of now, we can only handle the case where there's 2 more
	if (Boundary_Nodes.size() != 2) {
		// This simply means that you are on a T-shaped region, where you can find more than 2 extra nodes that are also on the interface
		// This case is a bit too complex, just skip this 
		return -1;
	}

	// Do not coarsen this edge if any of the other boundary nodes have been modified before
	for (auto q : Boundary_Nodes) {
		if (obliteratedNodes.count(q))
			return -1;
	}

	#if (COARSENDEBUG > 0)
		cout << "The two extra nodes on the sharing edge are:" << endl;
		for ( auto q : Boundary_Nodes ){
			vector<double> c = nodes.at( q );
			cout << c[0] << " " << c[1] << " " << c[2] << endl;
		}
	#endif


	// Step 2.2: Split all aux facets based on the branch they are in
	int total_found = 0;
	for (auto cfid : matchingFacets) { // Each matching facet is a branch
		vector<int> curr_branch;
		// Step 2.2.1: Check if this facet has 2 edges on the sharing boundary
		vector<int> cfn = globalFacets.at(cfid);
		if (std::find(cfn.begin(), cfn.end(), Boundary_Nodes[0]) != cfn.end()
				|| std::find(cfn.begin(), cfn.end(), Boundary_Nodes[1]) != cfn.end() ) {
			#if (COARSENDEBUG > 0)
				cout << "Facet " << cfid << " has 2 edges on the sharing boundary...." << endl;
			#endif
			// A special case needs to be implemented to take care of this, let's skip this for now
			return -1;
		}

		// Step 2.2.2: Find me the 2 facets that are directly connected to the facet we are working on now
		vector<int> Node_Memory;
		// Give me the node on the current facet that's not nid1 and nid2
		int a_node_on_sharing_facet;
		for (auto n : cfn) {
			if (n != nid1 && n != nid2) {
				a_node_on_sharing_facet = n;
				Node_Memory.push_back(n);
				break;
			}
		}
		// Search among the connected facets of a_node_on_sharing_facet
		for (auto fid : node2facet.at(a_node_on_sharing_facet)) {
			if (fid != cfid) {
				vector<int> n = globalFacets.at(fid);
				if (std::find(n.begin(), n.end(), nid1) != n.end()
						|| std::find(n.begin(), n.end(), nid2) != n.end()) {
					curr_branch.push_back(fid);
					UnifyConvention(cfid, fid); // Attempt to check convention
				}
			}
		}
		// If everthing goes well, we should find 2 facets connected to the current facet that we are working on, that are not in matchingFacets
		if (curr_branch.size() != 2) {
			cout << "*** ERROR: At this point we found " << curr_branch.size()
					<< " facets on this branch!" << endl;
			exit(1);
		}

		// Step 2.2.3: Continue our search in 2 directions (left & right)
		for (int ii = 0; ii < 2; ii++) {
			// For each direction, initiate our search on a_node_on_sharing_facet
			int a_node = a_node_on_sharing_facet;

			// Check if we hit the interface edge already with this current facet
			int on_this_facet = curr_branch[ii];
			vector<int> sfn = globalFacets.at(on_this_facet); // Sharing facet nodes
			if (std::find(sfn.begin(), sfn.end(), Boundary_Nodes[0]) != sfn.end()
					|| std::find(sfn.begin(), sfn.end(), Boundary_Nodes[1]) != sfn.end()) {
				// We hit the interface edge in this direction, end
				continue;
			}

			// Ok, let's search along this direction for more facets
			int another_node;
			bool HIT_END = false;
			int ccc = 0;
			while (1 && ccc < 100) {
				ccc += 1;
				for (auto nn : globalFacets.at(on_this_facet)) {
					// Find the node that's not nid1/2 and not a_node
					if (nn != a_node && nn != nid1 && nn != nid2
							&& std::find(Node_Memory.begin(), Node_Memory.end(),nn) == Node_Memory.end()) {
						another_node = nn;
						Node_Memory.push_back(nn);
						break;
					}
				}

				// Search in the connected facets of this newly found node
				for (auto fid : node2facet.at(another_node)) {
					if (fid == cfid) {
						continue;
					}
					vector<int> cfn = globalFacets.at(fid); // Get nodes on this current facet
					if ((std::find(cfn.begin(), cfn.end(), nid1) != cfn.end()
							|| std::find(cfn.begin(), cfn.end(), nid2) != cfn.end()) && // If it contains nid1/2
							std::find(curr_branch.begin(), curr_branch.end(), fid) == curr_branch.end()) { // If it is not currently in the branch
						curr_branch.push_back(fid); // Add it
						UnifyConvention(on_this_facet, fid); // Attempt to check convention
						on_this_facet = fid; // Move on to this facet

						// Check if this is the end
						if (std::find(cfn.begin(), cfn.end(), Boundary_Nodes[0]) != cfn.end()
							|| std::find(cfn.begin(), cfn.end(),Boundary_Nodes[1]) != cfn.end()) {
							// We hit the interface edge in this direction, end
							HIT_END = true;
							break;
						}
						break; // Should only find one at a time
					}
				}
				if (HIT_END) {
					break;
				}

				// Move on to the next node
				a_node = another_node;
			}
			if (HIT_END == false) {
				cout << "Failed to find all branches!" << endl;
				return -1;
			}
		}

		Split_auxFacets.push_back(curr_branch);
		total_found += curr_branch.size();
		
		#if (COARSENDEBUG > 0)
			cout << "I found " << curr_branch.size() << " facets on the current branch" << endl;
			for ( auto q : curr_branch ){
				cout << q << " " ;
			}
			cout << endl;
		#endif
	}

	// Check if we found all aux facets
	if (total_found != auxFacets.size()) {
		#if (COARSENDEBUG > 0)
		cout
				<< ">>> WARNING: We didn't find all aux facets! Skipping this edge..."
				<< endl;
		#endif
		return -1;
	}

	// Continue Step 2
	map<int, vector<int> > proposedFacets;
	set<int> Temp_Store;
	// Loop through all branches that we sorted
	for (int iter = 0; iter < Split_auxFacets.size(); iter++) {
		vector<int> LoopMe = Split_auxFacets[iter];
		#if (COARSENDEBUG > 0)
			cout << "********************************************* Working on branch " << iter << endl;
		#endif
		// Get all nodes on the external edge loop, organize them in CCW order
		vector<int> External_Loop;
		int num_edges = ExternalLoop(LoopMe, External_Loop, nid1, nid2, 1);
		if (num_edges == -1) {
			exit(1);
		} else if (num_edges == -2) {
			return -1;
		}

		// Check average edge length of the initial guess
		double avg_edge_length = 0.;
		for (auto nnn : External_Loop) {
			vector<double> p = nodes.at(nnn);
			avg_edge_length += calculateDistanceBetweenTwoDoubleVectors(p,
					mid_pt);
		}
		double ph;
		// Use the sizing function to decide if we should coarsen or not
		if (sizingFunction(avg_edge_length / (float) num_edges, mid_pt, ph, 0) != 1) {
			nodes[nid1] = old_nc;
			return -1;
		}

		if (debug(edge)) {
			cout << "External edge loop: " << endl;
			for (auto q : External_Loop) {
				vector<double> c = nodes.at(q);
				cout << "Node " << q << ": " << c[0] << "," << c[1] << ","
						<< c[2] << endl;
			}
		}

		// Step 3: Check for local concave polys
		vector<int> concave_loc;
		int local_concave = FindConcavity(External_Loop, F1, num_edges,concave_loc);

		// Step 5: Form all initial surface pairs
		// Step 5.1: Move node 1 to mid-side location
		nodes[nid1] = mid_pt;
		// Step 5.2: Make internal triangulations
		vector<vector<int> > facet_pair, initial_facet_pair;
		for (int i = 1; i < External_Loop.size() - 1; i++) {
			vector<int> curr_pair = { nid1, External_Loop[i - 1],
					External_Loop[i], nid1, External_Loop[i], External_Loop[i
							+ 1] };
			facet_pair.push_back(curr_pair);
			initial_facet_pair.push_back(curr_pair); // Keep an unmodified copy
		}

		// Step 6: Attempt to swap those surface pairs
		vector<vector<int> > final_facet_pair;
		int swapped = SwapFacets(facet_pair, initial_facet_pair, concave_loc,
				minAngle, maxAngle, normalTol, final_facet_pair);
		if (swapped == -1) {
			nodes[nid1] = old_nc;
			return -1;
		}

		// Step 7: Update connectivity
		// Collect all old facets
		LoopMe.push_back(matchingFacets[iter]);

		// Just use the old facet IDs on the new facets, we should have enough
		if (LoopMe.size() > final_facet_pair.size() + 1) {
			// Build proposedFacets
			vector<vector<int> > tmpFacet;
			for (int i = 0; i < final_facet_pair.size(); i++) {
				vector<int> tf1 = { final_facet_pair[i][0],
						final_facet_pair[i][1], final_facet_pair[i][2] };
				vector<int> tf2 = { final_facet_pair[i][3],
						final_facet_pair[i][4], final_facet_pair[i][5] };
				// Check if tf1 is in tmpFacet
				if (std::find(tmpFacet.begin(), tmpFacet.end(), tf1)
						== tmpFacet.end()) {
					tmpFacet.push_back(tf1);
				}
				// Check if tf2 is in tmpFacet
				if (std::find(tmpFacet.begin(), tmpFacet.end(), tf2)
						== tmpFacet.end()) {
					tmpFacet.push_back(tf2);
				}
			}
			#if (COARSENDEBUG > 0)
				cout << "Found " << tmpFacet.size() << " unique facets" << endl;
			#endif

			// Sanity check, see if the swapping algorithm created something unexpected
			if (tmpFacet.size() != initial_facet_pair.size() + 1) {
				#if (COARSENDEBUG > 0)
				cout << ">>> WARNING: Initially there are "
						<< initial_facet_pair.size() + 1
						<< " unique facets, now we have " << tmpFacet.size()
						<< endl;
				#endif
				RemoveFreeEdge(tmpFacet, External_Loop,
						initial_facet_pair.size());
			}

			for (int i = 0; i < tmpFacet.size(); i++) {
				proposedFacets[LoopMe[i]] = tmpFacet[i];
				#if (COARSENDEBUG > 0)
					cout << "New facet " << LoopMe[i] << " has nodes: " << tmpFacet[i][0] << "," << tmpFacet[i][1] << "," << tmpFacet[i][2] << endl;
				#endif
			}

			for (auto q : External_Loop) {
				Temp_Store.insert(q);
			}
		} else {
			cout << "*** ERROR: Somehow we ran out of facet IDs...." << endl;
			cout << "Old: " << auxFacets.size() << " vs. New: "
					<< final_facet_pair.size() + 1 << endl;
			exit(1);
		}
	}
	#if (COARSENDEBUG > 0)
		cout << "Total of " << proposedFacets.size() << " proposed facets" << endl;
	#endif

	// Step 7: Continue to update connectivity
	if (proposedFacets.size() > 0) {
		// Collect all old facets
		for (auto f : matchingFacets) {
			auxFacets.push_back(f);
		}

		UpdateMesh(auxFacets, proposedFacets, nid1, nid2);
		for (auto q : Temp_Store) {
			obliteratedNodes.insert(q);
		}
	}
	return 0;
}

void GlobalMesh::CONDENSE3(ofstream &outfile) {
	int num_condensed = 0;
	std::map<int, vector<int> >::iterator itr;
	std::map<int, vector<int> > node2facetIter(node2facet);
	for (itr = node2facetIter.begin(); itr != node2facetIter.end(); itr++) {
		int cnid = itr->first;
		// Check if the condense operation has already been run.
		if (obliteratedNodes.find(cnid) != obliteratedNodes.end()) {
			continue;
		}
		vector<int> FL = itr->second;
		vector<int> nn = nodeNeighbors.at(cnid);
		if (FL.size() == 3 && nn.size() == 3) { // This is the case if a node has only 3 neighbors
		// If it is on the interface, let's keep the center node, but adjust to a better position
			if (interfaceNodes.count(cnid)) {
				continue;
			}
			int new_fid = FL[0];

			// Unify ordering convention
			UnifyConvention(new_fid, FL[1]);
			UnifyConvention(new_fid, FL[2]);

			// Compute the external edge loop in default positive convention
			// There shoud only be 3 nodes, this will form the new facet that we will add to the mesh
			// To do this, we go through all the connected facets
			vector<vector<int> > temp_set;
			for (auto cfid : FL) {
				vector<int> fn = globalFacets.at(cfid); // Get all nodes on this facet
				vector<int> curr_vector;
				if (fn[0] == cnid) {
					curr_vector = {fn[1] , fn[2]};
				}
				else if ( fn[1] == cnid ) {
					curr_vector = {fn[2] , fn[0]};
				}
				else {
					curr_vector = {fn[0] , fn[1]};
				}
				temp_set.push_back(curr_vector);
			}
			// Fill the external edge loop
			vector<int> External_Loop = temp_set[0];
			int max_try = 10;
			int count = 0;
			while (External_Loop.size() != 3 && count <= max_try) {
				count += 1;
				for (int i = 1; i < 3; i++) {
					vector<int> cv = temp_set[i];
					if (cv[1] == External_Loop[0]) {
						External_Loop.insert(External_Loop.begin(), cv[0]);
						break;
					}
					if (cv[0] == External_Loop.back()) {
						External_Loop.push_back(cv[1]);
						break;
					}
				}
			}
			if (External_Loop.size() != 3) {
				cout << ">>> ERROR: We didn't find all 3 external nodes, skipping" << endl;
				continue;
			}
			num_condensed += 1;

			// Check if the facet formed by the 3 external nodes already exists
			bool facet_exists = false;
			vector<int> suspected_facets = node2facet.at(External_Loop[0]);
			for (auto sf : suspected_facets) {
				vector<int> sfn = globalFacets.at(sf);
				if (std::find(sfn.begin(), sfn.end(), External_Loop[1])
						!= sfn.end()
						&& std::find(sfn.begin(), sfn.end(), External_Loop[2])
								!= sfn.end()) {
					facet_exists = true;
					break;
				}
			}

			// Update mesh information
			for (int i = 0; i < 3; i++) {
				int myNeighbor = nn[i];

				// Update node2facet
				vector<int> cn2f = node2facet.at(myNeighbor);
				cn2f.erase(std::remove(cn2f.begin(), cn2f.end(), FL[1]),
						cn2f.end());
				cn2f.erase(std::remove(cn2f.begin(), cn2f.end(), FL[2]),
						cn2f.end());
				if (facet_exists == false
						&& std::find(cn2f.begin(), cn2f.end(), new_fid)
								== cn2f.end()) {
					cn2f.push_back(new_fid);
				}
				node2facet[myNeighbor] = cn2f;

				// Update nodeNeighbors
				vector<int> cnn = nodeNeighbors.at(myNeighbor);
				cnn.erase(std::remove(cnn.begin(), cnn.end(), cnid), cnn.end());
				nodeNeighbors[myNeighbor] = cnn;

				globalFacets.erase(FL[i]); // Remove the 3 central facets
			}

			// Replace 3 small central facets with one facet
			if (facet_exists == false) {
				globalFacets[new_fid] = External_Loop;
			}

			// Remove cnid
			nodes.erase(cnid);
			node2facet.erase(cnid);
			nodeNeighbors.erase(cnid);
			obliteratedNodes.insert(cnid);
		}
	}
	outfile << ">>> Removed " << num_condensed << " unnecessary caps" << endl;
}

// Junyan's note on 07/08/20
// This function intends to condense a quad with an internal node into 2 triangles
// This is not tested thoroughly and it's known to fail from time to time
// Needs to be revisited before using this function
void GlobalMesh::CONDENSE4(double minAngle, double maxAngle) {
	int num_condensed = 0;
	std::map<int, vector<int> >::iterator itr;
	for (itr = node2facet.begin(); itr != node2facet.end(); itr++) {
		int cnid = itr->first;
		vector<int> FL = itr->second;
		vector<int> nn = nodeNeighbors.at(cnid);
		if (FL.size() == 4 && nn.size() == 4) { // This is the case if a node has only 4 neighbors
			if (interfaceNodes.count(cnid)) {
				continue;
			}
			num_condensed += 1;
			int new_fid1 = FL[0];
			int new_fid2 = FL[1];

			// First, we attempt to remove the internal node and replace with two triangles
			// Compute a reference average normal
			vector<double> F1 = { 0, 0, 0 };
			for (auto q : FL) {
				vector<double> cf = { 0, 0, 0 };
				calculateFacetNormal(globalFacets.at(q), cf);
				F1[0] += cf[0] / 4.;
				F1[1] += cf[1] / 4.;
				F1[2] += cf[2] / 4.;
			}

			// Compute External_Loop with the right connectivity
			// To do this, we go through all the connected facets
			vector<vector<int> > temp_set;
			for (auto cfid : FL) {
				vector<int> fn = globalFacets.at(cfid); // Get all nodes on this facet
				vector<int> curr_vector;
				if (fn[0] == cnid) {
					curr_vector = {fn[1] , fn[2]};
				}
				else if ( fn[1] == cnid ) {
					curr_vector = {fn[2] , fn[0]};
				}
				else {
					curr_vector = {fn[0] , fn[1]};
				}
				temp_set.push_back(curr_vector);
			}

			// Fill the external edge loop
			vector<int> External_Loop = temp_set[0];
			while (External_Loop.size() != 4) {
				for (int i = 1; i < 4; i++) {
					vector<int> cv = temp_set[i];
					if (cv[1] == External_Loop[0]) {
						External_Loop.insert(External_Loop.begin(), cv[0]);
						break;
					}
					if (cv[0] == External_Loop.back()) {
						External_Loop.push_back(cv[1]);
						break;
					}
				}
			}

			// Step 3: Check for local concave polys
			int num_edges = 4;
			int local_concave = 0;
			vector<int> concave_loc;
			for (int i = 0; i < num_edges; i++) {
				// Make a local facet
				int idx2 = i + 1;
				int idx3 = i + 2;
				if (idx3 >= num_edges) {
					idx3 -= num_edges;
				}
				if (idx2 >= num_edges) {
					idx2 -= num_edges;
				}
				vector<int> test_facet = { External_Loop[i],
						External_Loop[idx2], External_Loop[idx3] };
				vector<double> TN = { 0, 0, 0 };
				calculateFacetNormal(test_facet, TN);

				// Compare this normal to the reference normal
				if (F1[0] * TN[0] + F1[1] * TN[1] + F1[2] * TN[2] < 0.) {
					local_concave += 1;
					concave_loc.push_back(idx2);
				}
			}
			// Step 4: Reorder External_Loop to start with the first local concavity
			if (local_concave > 0 && concave_loc[0] != 0) {
				vector<int> t = vector<int>(
						External_Loop.begin() + concave_loc[0],
						External_Loop.end());
				for (int i = 0; i < concave_loc[0]; i++) {
					t.push_back(External_Loop[i]);
				}
				External_Loop = vector<int>(t.begin(), t.end());
			}

			// Attempt to swap edges, there are only 2 options
			vector<int> o1f1, o1f2, o2f1, o2f2;
			vector<double> tf1_a = { 0, 0, 0 };
			vector<double> tf2_a = { 0, 0, 0 };

			// Check which one is better
			bool Opt1 = true;
			bool Opt2 = false;

			// Option 1
			o1f1 = {External_Loop[0] , External_Loop[1] , External_Loop[2]};
			o1f2 = {External_Loop[0] , External_Loop[2] , External_Loop[3]};
			calculateFacetAngles(o1f1, tf1_a);
			calculateFacetAngles(o1f2, tf2_a);
			double o1_min = std::min(tf1_a[0], tf2_a[0]);
			double o1_max = std::max(tf1_a[2], tf2_a[2]);

			if (local_concave == 0) { // Only check how good the second option is, when there is no concavity
				// Check how good is this option
				if (o1_min < minAngle || o1_max > maxAngle) {
					Opt1 = false;
				}

				// Option 2
				o2f1 = {External_Loop[0] , External_Loop[1] , External_Loop[3]};
				o2f2 = {External_Loop[1] , External_Loop[2] , External_Loop[3]};
				calculateFacetAngles(o2f1, tf1_a);
				calculateFacetAngles(o2f2, tf2_a);
				double o2_min = std::min(tf1_a[0], tf2_a[0]);
				double o2_max = std::max(tf1_a[2], tf2_a[2]);

				if (o2_min > o1_min && o2_min > minAngle && o2_max < maxAngle) {
					// Choose option 2
					Opt2 = true;
					Opt1 = false;
				}
			}

			// See if we can use any of the two
			if (Opt1 || Opt2) {
				vector<int> nf1, nf2;
				if (Opt1) {
					nf1 = o1f1;
					nf2 = o1f2;
				} else {
					nf1 = o2f1;
					nf2 = o2f2;
				}

				// Update mesh information
				// Step 7.1: node2facet and nodeNeighbors
				// Remove old facets in node2facet, remove old nodes in nodeNeighbors
				vector<int> nn2f;
				for (auto fid : FL) {
					vector<int> cf = globalFacets.at(fid);
					for (int i = 0; i < 3; i++) {
						nn2f = node2facet.at(cf[i]);
						nn2f.erase(std::remove(nn2f.begin(), nn2f.end(), fid),
								nn2f.end());
						node2facet[cf[i]] = nn2f;

						if (cf[i] != cnid) {
							nn2f = nodeNeighbors.at(cf[i]);
							// Remove all internal connections to cnid
							nn2f.erase(
									std::remove(nn2f.begin(), nn2f.end(), cnid),
									nn2f.end());
							nodeNeighbors[cf[i]] = nn2f;
						}
					}
				}
				// Step 7.2: globalFacets
				// Remove old facets in globalFacets
				for (auto fid : FL) {
					globalFacets.erase(fid);
				}
				// Step 7.3: proposedFacets
				for (int ii = 0; ii < 2; ii++) {
					int cfid;
					vector<int> cf;
					if (ii == 0) {
						cfid = new_fid1;
						cf = nf1;
					} else {
						cfid = new_fid2;
						cf = nf2;
					}

					// Add new facets
					globalFacets[cfid] = cf;

					// Update node2facet
					for (int i = 0; i < 3; i++) {
						nn2f = node2facet.at(cf[i]);
						if (std::find(nn2f.begin(), nn2f.end(), cfid)
								== nn2f.end()) { // Add the current facet if it does not exist
							nn2f.push_back(cfid);
						}
						node2facet[cf[i]] = nn2f;

						// Update nodeNeighbors
						nn2f = nodeNeighbors.at(cf[i]);
						for (int j = 0; j < 3; j++) {
							// Loop through all other nodes on the same facet, add to list if they don't already exist
							if ((i != j)
									&& std::find(nn2f.begin(), nn2f.end(),
											cf[j]) == nn2f.end()) {
								nn2f.push_back(cf[j]);
							}
						}
						nodeNeighbors[cf[i]] = nn2f;
					}
				}
				// Remove cnid
				nodes.erase(cnid);
				node2facet.erase(cnid);
				nodeNeighbors.erase(cnid);
				obliteratedNodes.insert(cnid);
			} else {
				// Too bad, none of those two options satisfy our criteria
				// Simply adjusting the nodal location of the interior node
				vector<double> avg_loc = { 0., 0., 0. };
				for (auto q : nn) {
					vector<double> cp = nodes.at(q);
					avg_loc[0] += cp[0] / 4.;
					avg_loc[1] += cp[1] / 4.;
					avg_loc[2] += cp[2] / 4.;
				}
				nodes[cnid] = avg_loc;
			}
		}
	}
	cout << "Adjusted " << num_condensed << " quadrilateral facets" << endl;
}

void GlobalMesh::FixHourglass(ofstream &outfile) {
	int num_hourglass = 0;
	int curr_extra_nid = -1;
	std::map<int, vector<double> >::iterator itr;

	// Loop through all the nodes
	for (itr = nodes.begin(); itr != nodes.end(); itr++) {
		int cnid = itr->first;
		vector<double> nc = itr->second;

		// Exclude boundary nodes
		if (interfaceNodes.count(cnid)) {
			continue;
		}

		// Get all facets connected to this node
		vector<int> facet_list = node2facet.at(cnid);

		// Attempt to divide all facets into two groups
		vector<int> AF = vector<int>(facet_list.begin() + 1, facet_list.end());
		int REF = facet_list[0];
		vector<int> FL1, FL2;

		// Fill the first group
		FL1.push_back(REF); // Give it a starting point
		while (1) {
			vector<int> rf = globalFacets.at(REF);
			int Old_size = FL1.size();
			// Check all facets in the list
			for (auto f : AF) {
				vector<int> cf = globalFacets.at(f);
				int t1 = std::find(rf.begin(), rf.end(), cf[0]) != rf.end();
				int t2 = std::find(rf.begin(), rf.end(), cf[1]) != rf.end();
				int t3 = std::find(rf.begin(), rf.end(), cf[2]) != rf.end();
				
				// If the current facet that we are checking shares an edge with the reference facet, they belong to the same group
				if (t1 + t2 + t3 == 2) {
					FL1.push_back(f);
					REF = f;
					AF.erase(std::remove(AF.begin(), AF.end(), f), AF.end());
					break;
				}
			}
			// We have found all the facets in this group when the size of FL1 does not change
			if (FL1.size() == Old_size) {
				break;
			}
		}

		// Fill the second group. Anything that's not in the first group will be in the second group
		for (auto f : facet_list) {
			if (std::find(FL1.begin(), FL1.end(), f) == FL1.end()) {
				FL2.push_back(f);
			}
		}

		// If we encounter a hourglass, we will find some facets in FL2. Because not all facets in a hourglass share an edge, some only share a node
		if (FL2.size() > 0) {
			// We found a possible hourglass
			num_hourglass += 1;

			// Get all nodes on the external loop of each group of facets
			vector<int> External_Loop, External_Loop2;
			for (auto f : FL1) {
				for (auto n : globalFacets.at(f)) {
					if (n != cnid
							&& std::find(External_Loop.begin(),
									External_Loop.end(), n)
									== External_Loop.end()) {
						External_Loop.push_back(n);
					}
				}
			}
			for (auto f : FL2) {
				for (auto n : globalFacets.at(f)) {
					if (n != cnid
							&& std::find(External_Loop2.begin(),
									External_Loop2.end(), n)
									== External_Loop2.end()) {
						External_Loop2.push_back(n);
					}
				}
			}

			// Calculate new positions for the two copies of the center node
			// We place them at the average location of the respective external polygons
			vector<double> avg1 = { 0., 0., 0. };
			vector<double> avg2 = { 0., 0., 0. };
			double n1 = (float) External_Loop.size();
			for (auto q : External_Loop) {
				vector<double> c = nodes.at(q);
				avg1[0] += c[0] / n1;
				avg1[1] += c[1] / n1;
				avg1[2] += c[2] / n1;
			}
			double n2 = (float) External_Loop2.size();
			for (auto q : External_Loop2) {
				vector<double> c = nodes.at(q);
				avg2[0] += c[0] / n2;
				avg2[1] += c[1] / n2;
				avg2[2] += c[2] / n2;
			}
			// Duplicate the node
			nodes[cnid] = avg1;
			nodes[curr_extra_nid] = avg2;

			// Update mesh information
			nodeNeighbors[cnid] = External_Loop;
			node2facet[cnid] = FL1;
			nodeNeighbors[curr_extra_nid] = External_Loop2;
			node2facet[curr_extra_nid] = FL2;
			for (auto f : FL2) {
				vector<int> fn = globalFacets.at(f);
				replace(fn.begin(), fn.end(), cnid, curr_extra_nid);
				globalFacets[f] = fn;
			}
			for (auto n : External_Loop2) {
				vector<int> nn = nodeNeighbors.at(n);
				replace(nn.begin(), nn.end(), cnid, curr_extra_nid);
				nodeNeighbors[n] = nn;
			}

			// Update node count
			curr_extra_nid -= 1;
		}
	}
	outfile << ">>> Done pre-conditioning, removed " << num_hourglass
			<< " hourglasses." << endl;
}

void GlobalMesh::CheckDuplicate() {
	nodeCloudPt facet_cloud;
	my_kd_tree_t dupliFacetKD(3 /*dim*/, facet_cloud,
			KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	int num_duplicate = 0;

	map<int, vector<int>>::iterator it;
	for (it = globalFacets.begin(); it != globalFacets.end(); it++) {
		int fid = it->first;
		vector<int> nids = it->second;
		vector<double> xyz0 = nodes.at(nids[0]);
		vector<double> xyz1 = nodes.at(nids[1]);
		vector<double> xyz2 = nodes.at(nids[2]);

		// Build a facet centroid, used in the KD tree search
		vector<double> xyz =
				{ (xyz0[0] + xyz1[0] + xyz2[0]) / 3.0, (xyz0[1] + xyz1[1]
						+ xyz2[1]) / 3.0, (xyz0[2] + xyz1[2] + xyz2[2]) / 3.0 };

		int closest_index = 0;
		int duplicate = isThisADuplicate(xyz, closest_index, dupliFacetKD);

		if (duplicate == 0) {
			// If Not Duplicate, place in KD-tree
			facet_cloud.nodeID.push_back(fid);
			facet_cloud.nodexyz.push_back(xyz);
			dupliFacetKD.addPoints(facet_cloud.getSize() - 1,
					facet_cloud.getSize() - 1);
		} else {
			// We found a duplicate facet
			num_duplicate += 1;
			// cout << xyz0[0] << "," << xyz0[1] << "," << xyz0[2] << endl;
			// cout << xyz1[0] << "," << xyz1[1] << "," << xyz1[2] << endl;
			// cout << xyz2[0] << "," << xyz2[1] << "," << xyz2[2] << endl;
			// exit(1);
		}
	}
	if (num_duplicate > 0) {
		cout << ">>> ERROR: Found " << num_duplicate << " duplicate facets!"
				<< endl;
		exit(1);
	}
}

// Split a string with a delimiter
// Function author: Vincenzo Pii
// Extracted from: https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
vector<string> split(string s, string delimiter) {
	size_t pos_start = 0, pos_end, delim_len = delimiter.length();
	string token;
	vector<string> res;

	while ((pos_end = s.find(delimiter, pos_start)) != string::npos) {
		token = s.substr(pos_start, pos_end - pos_start);
		pos_start = pos_end + delim_len;
		res.push_back(token);
	}

	res.push_back(s.substr(pos_start));
	return res;
}

void GlobalMesh::BuildCrackFrontKD(string FilePath, ofstream &outfile) {
	ifstream AdvCrk;
	string Line;
	int PtIDX = 0;

	AdvCrk.open(FilePath);
	vector<double> cp;
	vector<double> adv_centroid;
	while (getline(AdvCrk, Line)) {
		vector<string> tmp = split(Line, ",");

		// If we are reading the advancing crack file
		if (settings->defectType == 0) {
			cp = {stod(tmp[7])*Undo_shrink_factorX , stod(tmp[8])*Undo_shrink_factorY , stod(tmp[9])*Undo_shrink_factorZ};
			adv_centroid = { stod(tmp[1])*Undo_shrink_factorX, stod(tmp[2])*Undo_shrink_factorX, stod(tmp[3])*Undo_shrink_factorX  };
		} else {
			cp = {stod(tmp[1])*Undo_shrink_factorX , stod(tmp[2])*Undo_shrink_factorY , stod(tmp[3])*Undo_shrink_factorZ};
			adv_centroid = cp;
			double cr = stod(tmp[4])*((Undo_shrink_factorX + Undo_shrink_factorY + Undo_shrink_factorZ)/3.);
			VoidRadius.insert( std::pair< int , double >( PtIDX , cr ) );
		}

		CrackPtcloud.nodeID.push_back(PtIDX);
		CrackPtcloud.nodexyz.push_back(cp);
		CrackPtcloud.crack_init_pt.push_back(adv_centroid);
		CrackPtKD.addPoints(CrackPtcloud.getSize() - 1,
				CrackPtcloud.getSize() - 1);
		PtIDX += 1;
	}
	AdvCrk.close();
	if (PtIDX <= 1) {
		cout << ">>> ERROR: Did not read any points! Please double check your file path!"
				<< endl;
		exit(1);
	}
	if (settings->defectType == 0) {
		outfile << ">>> Read " << PtIDX << " crack front points" << endl;
	} else {
		outfile << ">>> Read " << PtIDX << " void centers" << endl;
	}

}

void GlobalMesh::UnifyConvention(int REF, int fid) {
	// Check which edge that the two shares
	const vector<int> f1 = globalFacets.at(REF);
	vector<int> f2 = globalFacets.at(fid);

	vector<int> SE; // Signed 
	if (std::find(f2.begin(), f2.end(), f1[0]) != f2.end()
			&& std::find(f2.begin(), f2.end(), f1[1]) != f2.end()) {
		SE = {f1[0] , f1[1]};
	}
	else if ( std::find( f2.begin() , f2.end() , f1[1] ) != f2.end() && std::find( f2.begin() , f2.end() , f1[2] ) != f2.end() ) {
		SE = {f1[1] , f1[2]};
	}
	else {
		SE = {f1[2] , f1[0]};
	}

	vector<int> SE1; // Signed 
	if (std::find(f1.begin(), f1.end(), f2[0]) != f1.end()
			&& std::find(f1.begin(), f1.end(), f2[1]) != f1.end()) {
		SE1 = {f2[0] , f2[1]};
	}
	else if ( std::find( f1.begin() , f1.end() , f2[1] ) != f1.end() && std::find( f1.begin() , f1.end() , f2[2] ) != f1.end() ) {
		SE1 = {f2[1] , f2[2]};
	}
	else {
		SE1 = {f2[2] , f2[0]};
	}

	// Check the ordering of this edge
	if (SE[0] == SE1[0] && SE[1] == SE1[1]) { // If the two facets share the same convention, the order of the sharing edge will be different in each facet
		globalFacets[ fid ] = {f2[2] , f2[1] , f2[0]}; // Invert the order
	}
}

// Define global sizing
// Mode: 0 : Surface mesh edge length;  1 : Volume mesh edge length
int GlobalMesh::sizingFunction(double currLen, vector<double> &mp,
		double& Objective_Length, int Mode) {
	// Flags:
	// -1 : Refine
	// 0 : No operation
	// 1 : Coarsen

	// Grad Scheme
	int GS = settings->gradationScheme;

	// No-coarsening zone size
	double R0 = settings->R0; // Radius of the refinement zone

	// Transition zone size
	double Rt = settings->Rt;

	// Desired max length
	double MaxEdge = settings->coar_targetSize;
	// Original edge length
	double OriginalEdge = settings->OriginalEdge;
	// Desired min length
	double MinEdge = OriginalEdge;
	if (refine_iterations > 0) {
		MinEdge = OriginalEdge / (2 * refine_iterations);
	}

	// For writing mtr file
	if (Mode == 1) { // If used in volume mesh mode
		MaxEdge = settings->coar_volTargetSize; // Give it some extra volume mesh gradation
		MinEdge = settings->ref_volTargetSize; // Give it some extra volume mesh refinement
	}

	// Search current point in KD tree
	const size_t num_results = 1;
	size_t ret_index;
	double out_dist_sqr;
	KNNResultSet<double> resultSet(num_results);
	resultSet.init(&ret_index, &out_dist_sqr);
	double query_pt[3] = { mp[0], mp[1], mp[2] };
	CrackPtKD.findNeighbors(resultSet, query_pt, nanoflann::SearchParams(10));
	double r = sqrt(out_dist_sqr);

	// Check ahead or behind crack
	// Check distance of current point and distance of adv point
	vector<double> init_pnt = CrackPtcloud.crack_init_pt.at(ret_index);
	vector<double> closest_cfn_pnt = CrackPtcloud.nodexyz.at(ret_index);
	double r_cfn = sqrt( pow(init_pnt[0] - closest_cfn_pnt[0],2) + pow(init_pnt[1] - closest_cfn_pnt[1],2) );
	double r_qyr = sqrt( pow(query_pt[0] - init_pnt[0],2) + pow(query_pt[1] - init_pnt[1],2) );

	// cout << "init_pnt: " << init_pnt[0] << " " << init_pnt[1] << " " << init_pnt[2] << endl;
	// cout << "closest_pnt: " << closest_pnt[0] << " " << closest_pnt[1] << " " << closest_pnt[2] << endl;
	// cout << "query: " << query_pt[0] << " " << query_pt[1] << " " << query_pt[2] << endl;
	// cout << r_cfn << " " << r_qyr << endl;

	// exit(1); 
	double ratio;

	// Gradation method for voids
	if (settings->defectType == 1) {
		// In this case, let's rescale R0 and Rt by the radius of void
		double void_radius = VoidRadius.at(ret_index);
		R0 *= void_radius;
		Rt *= void_radius;
	}

	if (GS == 1) {
		// Define refinement zone
		if (r <= R0) {
			// Refinement objective function
				ratio = pow(r / R0, settings->powerExp);
				Objective_Length = ratio * (OriginalEdge - MinEdge) + MinEdge;
			if (currLen > Objective_Length) {
				return -1;
			} // Refine
			else
				return 1; // Allow coarsening in the refinement region
		} else {
			ratio = min(1., abs(r - R0) / (Rt - R0));
			Objective_Length = ratio * (MaxEdge - OriginalEdge) + OriginalEdge;

			if (currLen < Objective_Length) {
				return 1;
			} // Coarsen
		}
	}
	else if (GS == 2) {
		if ((r_qyr < r_cfn) && ( (r_cfn - r_qyr) < R0  ) ){
			ratio = pow((r_cfn - r_qyr) / R0, settings->powerExp);
			Objective_Length = ratio * (MaxEdge - MinEdge) + MinEdge;
			if (currLen < Objective_Length) {
				return 1;
			} // Coarsen
			else
				return -1; // Allow coarsening in the refinement region
		} else {
			if (r <= R0) {
				Objective_Length = MinEdge;

			if (currLen > Objective_Length) {
				return -1;
			} // Refine
			else
				return 1; // Allow coarsening in the refinement region

			} else {
				ratio = min(1., abs(r - R0) / (Rt - R0));
				Objective_Length = ratio * (MaxEdge - OriginalEdge) + OriginalEdge;

				if (currLen < Objective_Length) {
					return 1;
				} // Coarsen
			}



		}
	} else {
		cerr << "Gradation scheme unrecognized: " << GS << endl;
		exit(1);
	}
	return 0; // No operation
}

void GlobalMesh::writeMTR(string filename, string folder) {
	size_t lastindex = filename.find_last_of(".");
	string FilePath = filename.substr(0, lastindex);

	vector<double> NodalEdgeLength;

	// Read from the .node file
	ifstream Node;
	string Line;
	int LineIDX = 0;
	int num_nodes = 0;

	Node.open(folder + '/' + FilePath + ".1.node");
	while (getline(Node, Line)) {
		LineIDX += 1;
		vector<string> tmp;
		tmp = split(Line, "  ");

		if (LineIDX == 1) {
			tmp = split(Line, "  ");
			num_nodes = stoi(tmp[0]);
		} else if (LineIDX != num_nodes + 2) {
			tmp = split(Line, "  ");
			int s = tmp.size();
			if (s < 3) {
				break;
			}
			// Read coordinates of nodes in the background volume mesh
			vector<double> cp = { stod(tmp[s - 3]), stod(tmp[s - 2]), stod(
					tmp[s - 1]) };
			double DesiredLen = 0.;

			// Use the sizing function to calculate the desired edge length at this point
			int ph = sizingFunction(0., cp, DesiredLen, 1);
			NodalEdgeLength.push_back(DesiredLen);
		}
	}
	Node.close();

	// Write to mtr file
	ofstream mtr;
	mtr.open(folder + '/' + FilePath + ".1.mtr");

	mtr << num_nodes << "  1" << endl;
	for (auto l : NodalEdgeLength) {
		mtr << l << endl;
	}
	mtr.close();
}

void GlobalMesh::RemoveFreeEdge(vector<vector<int> > &tmpFacet,
		vector<int> &External_Loop, int targetSize) {
	// Check if there is one with a free edge. If so, remove it
	map<string, vector<int> > edgeMap;
	for (int ii = 0; ii < tmpFacet.size(); ii++) {
		vector<int> f = tmpFacet[ii];
		vector<vector<int> > fe = { { f[0], f[1] }, { f[1], f[2] },
				{ f[2], f[0] } };
		for (auto e : fe) {
			string edge_name;
			if (e[0] > e[1]) {
				edge_name = to_string(e[0]) + "_" + to_string(e[1]);
			} else {
				edge_name = to_string(e[1]) + "_" + to_string(e[0]);
			}

			map<string, vector<int> >::iterator it = edgeMap.find(edge_name);
			if (it == edgeMap.end()) {
				vector<int> fl = { ii };
				edgeMap.insert(std::pair<string, vector<int> >(edge_name, fl));
			} else {
				edgeMap[edge_name].push_back(ii);
			}
		}
	}

	vector<string> External_Edge;
	External_Loop.push_back(External_Loop[0]);
	for (int ii = 0; ii < External_Loop.size() - 1; ii++) {
		vector<int> e = { External_Loop[ii], External_Loop[ii + 1] };
		string edge_name;
		if (e[0] > e[1]) {
			edge_name = to_string(e[0]) + "_" + to_string(e[1]);
		} else {
			edge_name = to_string(e[1]) + "_" + to_string(e[0]);
		}
		External_Edge.push_back(edge_name);
	}

	map<string, vector<int> >::iterator it;
	for (it = edgeMap.begin(); it != edgeMap.end(); it++) {
		string edge_name = it->first;
		vector<int> List = it->second;
		if (List.size() == 1
				&& std::find(External_Edge.begin(), External_Edge.end(),
						edge_name) == External_Edge.end()) {
			tmpFacet.erase(tmpFacet.begin() + List[0]);
		#if (COARSENDEBUG > 0)
			cout
					<< ">>> WARNING: Removed one facet with free edge, make sure you know what you are doing!"
					<< endl;
		#endif
		}
	}

	if (tmpFacet.size() != targetSize + 1) {
		cout << ">>> ERROR: Size mismatch after removing free edges: "
				<< targetSize + 1 << " vs. " << tmpFacet.size() << endl;
		exit(1);
	}
}

bool GlobalMesh::CheckOverlap(vector<vector<int> > &tmpFacet,
		vector<int> &External_Loop, vector<int> &AllInternalFacets) {
	set<int> Checked;
	for (auto n : External_Loop) {
		// Get all facets connected to this node
		vector<int> cf;
		cf = node2facet.at(n);
		
		for (auto cfid : cf) {
			if (Checked.count(cfid)) {
				continue;
			}
			if (std::find(AllInternalFacets.begin(), AllInternalFacets.end(),
					cfid) != AllInternalFacets.end()) {
				continue;
			}
			vector<int> fn;
			fn = globalFacets.at(cfid);
			
			int t1 = std::find(External_Loop.begin(), External_Loop.end(),
					fn[0]) != External_Loop.end();
			int t2 = std::find(External_Loop.begin(), External_Loop.end(),
					fn[1]) != External_Loop.end();
			int t3 = std::find(External_Loop.begin(), External_Loop.end(),
					fn[2]) != External_Loop.end();
			if (t1 + t2 + t3 == 3) {
				for (auto tf : tmpFacet) {
					int t1_ = std::find(tf.begin(), tf.end(), fn[0])
							!= tf.end();
					int t2_ = std::find(tf.begin(), tf.end(), fn[1])
							!= tf.end();
					int t3_ = std::find(tf.begin(), tf.end(), fn[2])
							!= tf.end();
					if (t1_ + t2_ + t3_ == 3) {
						return true;
					}
				}
			}
			Checked.insert(cfid);
		}
	}
	return false;
}

void GlobalMesh::FixBoundary(ofstream &outfile) {
	outfile << ">>> Begin sawtooth edge smoothing" << endl;

	// User inputs:
	int N_smooth_iters = settings->EdgeSmoothing;

	// Build noundary node map
	map<int, set<int> > nodeMap;
	map<int, set<int> >::iterator iter;
	for (int i = 0; i < jEdges.size(); i++) {
		vector<int> en = jEdges[i].n;
		for (int i = 0; i < 2; i++) {
			iter = nodeMap.find(en[i]);
			if (iter == nodeMap.end()) {
				nodeMap.insert(std::pair<int, set<int> >(en[i], { en[1 - i] }));
			} else {
				nodeMap[en[i]].insert(en[1 - i]);
			}
		}
	}

	set<int> Processed_Nodes;

	// Preprocess
	for (iter = nodeMap.begin(); iter != nodeMap.end(); iter++) {
		int nid = iter->first;
		set<int> myNeighbors = iter->second;
		if (myNeighbors.size() != 2) {
			continue;
		} // A normal edge node has two neighbors

		vector<int> tmp;
		for (auto n : myNeighbors) {
			tmp.push_back(n);
		}
		int nn1 = tmp[0];
		int nn2 = tmp[1];
		set<int> neighborSet1 = nodeMap.at(nn1);
		set<int> neighborSet2 = nodeMap.at(nn2);

		vector<int> fl = node2facet.at(nid);
		vector<int> fl1 = node2facet.at(nn1);
		vector<int> fl2 = node2facet.at(nn2);

		int shared_facet = 0;
		for (auto f : fl) {
			if (std::find(fl1.begin(), fl1.end(), f) != fl1.end()
					&& std::find(fl2.begin(), fl2.end(), f) != fl2.end()) {
				shared_facet = f;
				break;
			}
		}

		if (shared_facet != 0) {
			int other_facet, other_node = 0;
			for (auto f : fl1) {
				vector<int> fn = globalFacets.at(f);
				if (std::find(fn.begin(), fn.end(), nn2) != fn.end()
						&& f != shared_facet) {
					other_facet = f;
					for (auto n : fn) {
						if (n != nn1 && n != nn2) {
							other_node = n;
							break;
						}
					}
					break;
				}
			}
			if (other_facet == 0 || other_node == 0) {
				continue;
			}

			// Update mesh
			// Remove from node2facet
			vector<int> n2f;
			vector<int> involvedNodes = { nid, nn1, nn2, other_node };
			for (auto in : involvedNodes) {
				n2f = node2facet.at(in);
				n2f.erase(std::remove(n2f.begin(), n2f.end(), shared_facet),
						n2f.end());
				n2f.erase(std::remove(n2f.begin(), n2f.end(), other_facet),
						n2f.end());
				node2facet[in] = n2f;
			}

			// Remove from nodeNeighbors
			vector<int> nn;
			nn = nodeNeighbors.at(nn1);
			nn.erase(std::remove(nn.begin(), nn.end(), nn2), nn.end());
			nodeNeighbors[nn1] = nn;
			nn = nodeNeighbors.at(nn2);
			nn.erase(std::remove(nn.begin(), nn.end(), nn1), nn.end());
			nodeNeighbors[nn2] = nn;

			// Swap the facets
			globalFacets[ shared_facet ] = {nn1 , nid , other_node};
			globalFacets[ other_facet ] = {other_node , nid , nn2};

			// Add to node2facet
			vector<int> involvedFacets = { shared_facet, other_facet };
			for (auto f : involvedFacets) {
				for (auto n : globalFacets.at(f)) {
					n2f = node2facet.at(n);
					if (std::find(n2f.begin(), n2f.end(), f) == n2f.end()) {
						n2f.push_back(f);
					}
					node2facet[n] = n2f;
				}
			}

			// Add to nodeNeighbors
			nn = nodeNeighbors.at(nid);
			if (std::find(nn.begin(), nn.end(), other_node) == nn.end()) {
				nn.push_back(other_node);
			}
			nodeNeighbors[nid] = nn;
			nn = nodeNeighbors.at(other_node);
			if (std::find(nn.begin(), nn.end(), nid) == nn.end()) {
				nn.push_back(nid);
			}
			nodeNeighbors[other_node] = nn;
		}
	}

	// Ok, actually move the nodes here 
	// Make a map to store all the temp locations
	map<int, vector<double> > proposedNodeCoords;
	map<int, vector<double> >::iterator propNodeIter;

	for ( int smooth_iter = 0 ; smooth_iter < N_smooth_iters ; smooth_iter ++ ){
		outfile << ">>> Smoothing iteration " << smooth_iter + 1 << endl;
		for (iter = nodeMap.begin(); iter != nodeMap.end(); iter++) {
			int nid = iter->first;
			set<int> myNeighbors = iter->second;
			if (myNeighbors.size() != 2 || Processed_Nodes.count(nid)) {
				continue;
			} // A normal edge node has two neighbors

			vector<int> p2c; // All points to check for colinearity
			for (auto n : myNeighbors) {
				p2c.push_back(n);
			}
			int nn1 = p2c[0];
			int nn2 = p2c[1];

			vector<double> p0 = nodes.at(nid);
			vector<double> p1 = nodes.at(nn1);
			vector<double> p2 = nodes.at(nn2);

			// Marked modified nodes
			Processed_Nodes.insert(nid);

			// If we pass the test, compute new coords here
			// This is actually Laplacian smoothing of the edges...
			vector<double> tmp_pt = { 0.25 * (p1[0] + p2[0]) + 0.5 * p0[0],
					0.25 * (p1[1] + p2[1]) + 0.5 * p0[1],
					0.25 * (p1[2] + p2[2]) + 0.5 * p0[2] };

			// Check for potential facet inversion
			bool facet_inverted = false;
			// Iterate all connected facets
			for ( auto f : node2facet.at( nid ) ){
				// Original normal
				vector<double> F = { 0, 0, 0 };
				vector<int> cf = globalFacets.at( f );
				calculateFacetNormal(cf, F);

				// Modified facet
				nodes[0] = tmp_pt;
				replace(cf.begin(), cf.end(), nid, 0);
				// Modified normal
				vector<double> F2 = { 0, 0, 0 };
				calculateFacetNormal(cf, F2);

				if ( F[0]*F2[0] + F[1]*F2[1] + F[2]*F2[2] < 0. ){
					facet_inverted = true;
					break;
				}
			}
			if ( facet_inverted ){ continue; }

			proposedNodeCoords.insert( std::pair<int, vector<double> >( nid , tmp_pt ) );
		}
		
		for ( propNodeIter = proposedNodeCoords.begin() ; propNodeIter != proposedNodeCoords.end() ; propNodeIter ++ ){
			nodes[ propNodeIter->first ] = propNodeIter->second;
		}
		outfile << "	Adjusted " << proposedNodeCoords.size() << " node coordinates" << endl;

		Processed_Nodes.clear();
		proposedNodeCoords.clear();
	}

	outfile << ">>> Sawtooth edge smoothing done!" << endl;
}

} /* namespace std */
