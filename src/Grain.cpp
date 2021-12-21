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

#include "Grain.h"

namespace std {

Grain::Grain(std::string fname) {
	std::vector<std::string> splitted;

	boost::split(splitted, fname, boost::is_any_of("._"));

	filename = fname;
	grainNumber = std::stoi(splitted[1]);
	prefix = splitted[0];

}

Grain::Grain(string stlFilename, int grainNumber, string stlPrefix) {

	filename = stlFilename;
	this->grainNumber = grainNumber;
	prefix = stlPrefix;

}

int Grain::convertStl(string folder) {
	string original_cwd = boost::filesystem::current_path().generic_string();
	int status = chdir(folder.c_str());
	string com = "python ../stl_binary2ascii.py -j " + filename
			+ " > /dev/null";
	const char *command = com.c_str();
	status = system(command);
	string out = "..";
	status = chdir(out.c_str());

	filename = "ascii" + filename;
	prefix = "ascii" + prefix;

	status = chdir(original_cwd.c_str());

	return 0;
}

void Grain::writeStl() {
	ofstream stlfile;
	stlfile.open(filename);
}

int Grain::inInternalNodeMap(int localNID, int &outputNID,
		int &global_nid_count) {
	// We'll set outputNID if the value exists in the internalNodeMap (i.e. already assigned a global value).
	// We output 1.
	// Otherwise, we return 0 (We should set this internal node ID).

	map<int, int>::iterator it = internalNodeMap.find(localNID);
	if (it != internalNodeMap.end()) {
		//element found;
		outputNID = it->second;
		return 1;
	}

	// Otherwise, this is an internal node not already accounted for.
	outputNID = global_nid_count;
	// Add it to the map
	internalNodeMap.insert(pair<int, int>(localNID, global_nid_count));

	return 0;
}

int Grain::getKDMap(int localNID, int &outputNID) {
	// We'll set outputNID if the value exists in the KDTreeMap (i.e. is a surface node).
	// If it is a surface node, we output 1.
	// Otherwise, we return 0 (this is an internal node and should get a new NID).

	map<int, int>::iterator it = kdTreeMap.find(localNID);
	if (it != kdTreeMap.end()) {
		//element found;
		outputNID = it->second;
		return 1;
	}

	// Otherwise, this is an internal node.
	//outputNID = 0;
	return 0;
}

int Grain::meshGrain(string folder, string cont_str, string ref_str, int Mode,
		programSettings &settings) {
	map<int, vector<vector<double> > > hole_list;
	int junk = 0;
	return this->meshGrain(folder, cont_str, ref_str, Mode, settings, junk,
			hole_list);
}

int Grain::meshGrain(string folder, string cont_str, string ref_str, int Mode,
		programSettings &settings, int &detectIsland,
		map<int, vector<vector<double> > > &hole_list) {
	//string ref = "1.1";
	bool exitfunction = false;
	string ref = cont_str;
	if (grainNumber < 0) {
		ref = ref_str + "a" + to_string(settings.ref_volTargetSize);
	}

	string com;
	size_t lastindex = filename.find_last_of(".");
	string base_file = filename.substr(0, lastindex);

	if ((detectIsland == 1) && (Mode == 2)) {
		// See if this grain ID is in the set:
		int rawGrainNum = this->giveRawGrainNumber();
		std::map<int, vector<vector<double> > >::iterator it = hole_list.find(
				rawGrainNum);

		if (it != hole_list.end()) {
			string sourcefile = folder + "/" + base_file + ".1.smesh";
			string sinkfile = folder + "/" + base_file + ".temp.smesh";

			// Open files
			std::ifstream source;
			source.open(sourcefile);
			std::ofstream sink;
			sink.open(sinkfile);

			string line;
			string holeline = "# part 3: hole list.";
			while (getline(source, line)) {
				if (line == holeline) {
					sink << line << "\n";
					vector<vector<double> > cent = hole_list.at(rawGrainNum);
					int sizeOfList = hole_list.at(rawGrainNum).size();
					sink << sizeOfList << "\n";
					for (int h = 0; h < sizeOfList; h++) {
						sink << 10000000 + h << " " << cent[h][0] << " "
								<< cent[h][1] << " " << cent[h][2] << "\n";
					}
					sink
							<< "\n# part 4: region list.\n0\n# Generated by tetgen -Y ./stls/asciiSV1Feature_2.stl and crackMesher";

					break;

				} else {
					sink << line << "\n";
				}

			}
			sink.close();
			source.close();

			boost::filesystem::path sourcepath = sourcefile;
			boost::filesystem::path sinkpath = sinkfile;
			boost::filesystem::remove(sourcefile);
			boost::filesystem::rename(sinkpath, sourcepath);
		}

	}

	while (exitfunction == false) {

		if (Mode == 2) {
			com = "tetgen -pmY " + folder + "/" + base_file + ".1.smesh -q"
					+ ref + " -V > " + base_file + ".txt";
		} else {
			com = "tetgen -pY " + folder + "/" + filename + " -q" + ref
					+ " -V > " + base_file + ".txt";
		}

		// cout << com << endl;
		const char *command = com.c_str();
		int status = system(command);

		bool failed = false;
		// First, check if the expected files even exist.
		if (!boost::filesystem::exists(
				folder + "/" + base_file + "." + to_string(Mode) + ".edge")) {
			std::cout << endl << "*** VOLUME MESHING ERROR: STL FILE "
					<< filename << " IS MISSING ITS EDGE FILE.";
			failed = true;
		}

		if (!boost::filesystem::exists(
				folder + "/" + base_file + "." + to_string(Mode) + ".ele")) {
			std::cout << endl << "*** VOLUME MESHING ERROR: STL FILE "
					<< filename << " IS MISSING ITS ELEMENT FILE.";
			failed = true;
		}

		if (!boost::filesystem::exists(
				folder + "/" + base_file + "." + to_string(Mode) + ".face")) {
			std::cout << endl << "*** VOLUME MESHING ERROR: STL FILE "
					<< filename << " IS MISSING ITS FACE FILE.";
			failed = true;
		}

		if (!boost::filesystem::exists(
				folder + "/" + base_file + "." + to_string(Mode) + ".node")) {
			std::cout << endl << "*** VOLUME MESHING ERROR: STL FILE "
					<< filename << " IS MISSING ITS NODE FILE.";
			failed = true;
		}

		if (failed == true) {
			cout << endl
					<< "*** Check your STL file to ensure it is volume meshable."
					<< endl;
			cout << "Failed command: " << command << endl;
			string response;
			cout
					<< "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Select from the following !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
					<< endl;
			cout << "S (skip meshing grain): If S, Grain " << grainNumber
					<< " will NOT be part of the mesh." << endl;
			cout
					<< "T (terminate): Cancel mesh generation completely and exit program."
					<< endl;
			cout
					<< "New quality measure (e.g. 2/0a1): enter another potential -q constraint for THIS GRAIN ONLY."
					<< endl;
			cout
					<< "         Logic dictates you enter laxer quality constraints such that the volume meshing succeeds."
					<< endl;
			cout << "Response: ";
			cin >> response;
			if (boost::iequals(response, "S")) {
				return -1;
			} else if (boost::iequals(response, "T")) {
				cout << "*** VOLUME MESHING FAILED. Exiting." << endl;
				exit(-1);
			} else {
				ref = response;
				cout << "New Command: "
						<< "tetgen -pY " + folder + "/" + filename + " -q" + ref
								+ " > /dev/null" << endl;
			}

		} else {
			if (Mode == 2) {
				boost::filesystem::copy_file(
						folder + "/" + base_file + ".1.smesh",
						folder + "/" + base_file + ".2.smesh");
			}
			exitfunction = true;
		}
	}

	return 0;
}

int Grain::grabSurfaceNodes(string folder, int Mode) {
	// cout << "Reading: " << prefix << "_" << to_string(grainNumber) << ".2" << endl;
	TetGen tetgenOut(
			prefix + "_" + to_string(grainNumber) + "." + to_string(Mode),
			folder);
	tetgenOut.sMeshReader(surfaceNodes, localFacets);
	tetgenOut.nodeReader(localNodes);
	tetgenOut.elementReader(localElements);
	return 0;
}

int Grain::writeSmshAndRead(string folder) {
	string com = "tetgen -Y " + folder + "/" + filename + " > /dev/null";
	//cout << com << endl;
	const char *command = com.c_str();
	int status = system(command);

	// Generated temporary smsh
	TetGen tetgenOut(prefix + "_" + to_string(grainNumber) + ".1", folder);
	tetgenOut.sMeshReader(surfaceNodes, localFacets);
	tetgenOut.nodeReader(localNodes);

	return 0;
}

int Grain::importBinary(ifstream &in) {
	// Grab the sizes of the arrays we are working with
	int sizeSurfaceNodes;
	int sizeLocalFacets;
	int sizeLocalNodes;
	in.read(reinterpret_cast<char*>(&sizeSurfaceNodes), sizeof(sizeSurfaceNodes));
	in.read(reinterpret_cast<char*>(&sizeLocalFacets), sizeof(sizeLocalFacets));
	in.read(reinterpret_cast<char*>(&sizeLocalNodes), sizeof(sizeLocalNodes));

	// Initialize our vector of surface nodes
	for (int i = 0; i < sizeSurfaceNodes; i++) {
		int temp;
		in.read(reinterpret_cast<char*>(&temp), sizeof(temp));
		surfaceNodes.push_back(temp);
	}

	// Read our local facets
	for (int i = 0; i < sizeLocalFacets; i++) {
		int tempFID;
		int nid1, nid2, nid3;
		in.read(reinterpret_cast<char*>(&tempFID), sizeof(tempFID));
		in.read(reinterpret_cast<char*>(&nid1), sizeof(nid1));
		in.read(reinterpret_cast<char*>(&nid2), sizeof(nid2));
		in.read(reinterpret_cast<char*>(&nid3), sizeof(nid3));
		vector<int> temp{ nid1, nid2, nid3 };
		localFacets.insert( pair<int, vector<int> >(tempFID, temp)  );
	}

	// Read our local nodes
	for (int i = 0; i < sizeLocalNodes; i++) {
		int tempNID;
		double x, y, z;
		in.read(reinterpret_cast<char*>(&tempNID), sizeof(tempNID));
		in.read(reinterpret_cast<char*>(&x), sizeof(x));
		in.read(reinterpret_cast<char*>(&y), sizeof(y));
		in.read(reinterpret_cast<char*>(&z), sizeof(z));
		vector<double> temp{x,y,z};
		localNodes.insert( pair<int, vector<double> >(tempNID, temp)  );
	}

	return 0;
}

int Grain::exportBinary(string filename) {
	std::ofstream out(filename, ios_base::out | ios_base::app );
	// Write: Grain ID, Length of Filename, Filename, Length of Prefix, Prefix
	out.write(reinterpret_cast<const char*>(&grainNumber), sizeof(grainNumber));
	int filenameSize = this->filename.size();
	out.write(reinterpret_cast<const char*>(&filenameSize), sizeof(filenameSize));
	out.write(&this->filename[0], filenameSize);
	int prefixSize = prefix.size();
	out.write(reinterpret_cast<const char*>(&prefixSize), sizeof(prefixSize));
	out.write(&prefix[0], prefixSize);
	// Write: Number of Surface Nodes, Number of Local Facets, Number of Local Nodes
	int sizeSurfaceNodes = surfaceNodes.size();
	int sizeLocalFacets = localFacets.size();
	int sizeLocalNodes = localNodes.size();

	out.write(reinterpret_cast<const char*>(&sizeSurfaceNodes), sizeof(sizeSurfaceNodes));
	out.write(reinterpret_cast<const char*>(&sizeLocalFacets), sizeof(sizeLocalFacets));
	out.write(reinterpret_cast<const char*>(&sizeLocalNodes), sizeof(sizeLocalNodes));
	// Write out surface nodes.
	for (auto it = begin (surfaceNodes); it != end (surfaceNodes); ++it) {
	    int nodeid = *it;
	    out.write(reinterpret_cast<const char*>(&nodeid), sizeof(nodeid));
	}
	// Write out all the local facets
	for (auto it = begin (localFacets); it != end (localFacets); ++it) {
	    int facetid = it->first;
	    int nid1 = it->second[0];
	    int nid2 = it->second[1];
	    int nid3 = it->second[2];
	    out.write(reinterpret_cast<const char*>(&facetid), sizeof(facetid));
	    out.write(reinterpret_cast<const char*>(&nid1), sizeof(nid1));
	    out.write(reinterpret_cast<const char*>(&nid2), sizeof(nid2));
	    out.write(reinterpret_cast<const char*>(&nid3), sizeof(nid3));
	}
	// Write out all the local nodes
	for (auto it = begin (localNodes); it != end (localNodes); ++it) {
	    int nodeid = it->first;
	    double x = it->second[0];
	    double y = it->second[1];
	    double z = it->second[2];
	    out.write(reinterpret_cast<const char*>(&nodeid), sizeof(nodeid));
	    out.write(reinterpret_cast<const char*>(&x), sizeof(x));
	    out.write(reinterpret_cast<const char*>(&y), sizeof(y));
	    out.write(reinterpret_cast<const char*>(&z), sizeof(z));
	}
	out.close();
	return 0;
}


int Grain::writeBackGroundSmsh(string folder, programSettings &settings) {
	string com;

	com = "tetgen -Y " + folder + "/" + filename + " -q"
			+ settings.backgroundQuality + " > /dev/null";

	// cout << com << endl;
	const char *command = com.c_str();
	int status = system(command);

	size_t lastindex = filename.find_last_of(".");
	string base_file = filename.substr(0, lastindex);

	bool failed = false;
	// First, check if the expected files even exist.
	if (!boost::filesystem::exists(folder + "/" + base_file + ".1.edge")) {
		std::cout << endl << "*** VOLUME MESHING ERROR: STL FILE " << filename
				<< " IS MISSING ITS EDGE FILE.";
		failed = true;
	}

	if (!boost::filesystem::exists(folder + "/" + base_file + ".1.ele")) {
		std::cout << endl << "*** VOLUME MESHING ERROR: STL FILE " << filename
				<< " IS MISSING ITS ELEMENT FILE.";
		failed = true;
	}

	if (!boost::filesystem::exists(folder + "/" + base_file + ".1.face")) {
		std::cout << endl << "*** VOLUME MESHING ERROR: STL FILE " << filename
				<< " IS MISSING ITS FACE FILE.";
		failed = true;
	}

	if (!boost::filesystem::exists(folder + "/" + base_file + ".1.node")) {
		std::cout << endl << "*** VOLUME MESHING ERROR: STL FILE " << filename
				<< " IS MISSING ITS NODE FILE.";
		failed = true;
	}

	if (!boost::filesystem::exists(folder + "/" + base_file + ".1.smesh")) {
		std::cout << endl << "*** VOLUME MESHING ERROR: STL FILE " << filename
				<< " IS MISSING ITS SMESH FILE.";
		failed = true;
	}

	if (failed == true) {
		cout << endl
				<< "*** Check your STL file to ensure it is volume meshable."
				<< endl;
		return -1;
	}
	return 0;
}

Grain::~Grain() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
