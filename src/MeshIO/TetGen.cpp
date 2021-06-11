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


#include "TetGen.h"

namespace std {

string cleanString(string str)
// https://www.geeksforgeeks.org/remove-extra-spaces-string/
		{
	// n is length of the original string
	int n = str.length();
	int i = 0, j = -1;
	bool spaceFound = false;
	while (++j < n && str[j] == ' ')
		;
	while (j < n) {
		// if current characters is non-space
		if (str[j] != ' ') {
			// remove preceding spaces before dot,
			// comma & question mark
			if ((str[j] == '.' || str[j] == ',' || str[j] == '?') && i - 1 >= 0
					&& str[i - 1] == ' ')
				str[i - 1] = str[j++];

			else
				// copy current character at index i
				// and increment both i and j
				str[i++] = str[j++];

			// set space flag to false when any
			// non-space character is found
			spaceFound = false;
		}
		// if current character is a space
		else if (str[j++] == ' ') {
			// If space is encountered for the first
			// time after a word, put one space in the
			// output and set space flag to true
			if (!spaceFound) {
				str[i++] = ' ';
				spaceFound = true;
			}
		}
	}

	// Remove trailing spaces
	if (i <= 1)
		str.erase(str.begin() + i, str.end());
	else
		str.erase(str.begin() + i - 1, str.end());
	return str;

}

TetGen::TetGen(string prefix, string folder) {
	// TODO Auto-generated constructor stub
	this->filename = prefix;
	this->folder = folder;
}

int TetGen::sMeshReader(vector<int>& surfaceNodes, map<int, vector<int> >& facetMap) {
	// cout << ">>> Reading in surface nodes from smesh file: "
	// 		<< folder + "/" + filename + ".smesh" << endl;
	ifstream infile(folder + "/" + filename + ".smesh");

	int nnodes = 0;
	int facetID = 1;

	std::string line;
	int mode = 0;
	bool set = false;
	vector<int> identifier;
	while (std::getline(infile, line)) {
		// cout << line << endl;
		if (boost::starts_with(line, "# part 2")) {
			//cout << ">>> Reading in nodes from smesh file: " << filename << endl;
			mode = 1;
			continue;
		} else if (boost::starts_with(line, "# part 3")) {
			//cout << ">>> Done reading surface nodes from smesh file: " << filename << endl;
			mode = 0;
			continue;
		} else if (line.compare("\n") == 0) {
			mode = 0;
			continue;
		} else if (line.length() < 3) {
			mode = 0;
			continue;
		}

		if (mode == 1) {
			if (set == false) {
				set = true;
			} else {
				parseSurfNodes(surfaceNodes,facetMap, facetID, line);
				facetID++;
				nnodes++;
			}
		}

	};

	//cout << ">>> Number of surface nodes read: " << surfaceNodes.size() << endl;
	if (surfaceNodes.size() == 0) {
		cout << endl << "*** SURFACE MESH READING ERROR: " << folder + "/" + filename + ".smesh" << " CONTAINS ZERO SURFACE MESH NODES." << endl;
		cout << "*** This may indicate that you either need to convert your STLs from binary to ascii (using the -b1 flag) or a error in your ascii files." << endl;
		exit(-1);
	}

	infile.close();

	return 0;
}


int TetGen::nodeReader(map<int, vector<double> >& nodes) {
	ifstream infile(folder + "/" + filename + ".node");

	int nnodes = 0;

	std::string line;
	bool set = false;
	vector<int> identifier;
	while (std::getline(infile, line)) {
		if (boost::starts_with(line, "#")) {
			continue;
		} else {
			if (set == false) {
				set = true;
			} else {
				line = line + ' ';
				line = cleanString(line);
				//cout << line << endl;
				std::stringstream linestream(line);

				std::string value;

				int nid, idx = 0;
				vector<double> xyz;
				while (getline(linestream, value, ' ')) {
					//cout << value << endl;
					if (idx == 0) {
						nid = stoi(value);
						idx++;
					} else {
						xyz.push_back(stod(value));
					}
				}
				nodes.insert(pair<int, vector<double> >(nid, xyz));
				nnodes++;
			}
		}
	}


	if (nnodes == 0) {
		cout << endl << "*** NODE READING ERROR: " << folder + "/" + filename + ".node" << " CONTAINS ZERO NODES." << endl;
		cout << "*** CHECK STL OR MESHING PARAMETERS." << endl;
		exit(-1);
	}
	infile.close();
	return 0;
}

int TetGen::elementReader(map<int, vector<int> >& elements) {
	ifstream infile(folder + "/" + filename + ".ele");

	int nels = 0;

	std::string line;
	bool set = false;
	vector<int> identifier;
	while (std::getline(infile, line)) {
		if (boost::starts_with(line, "#")) {
			continue;
		} else {
			if (set == false) {
				set = true;
			} else {
				line = line + ' ';
				line = cleanString(line);
				//cout << line << endl;
				std::stringstream linestream(line);

				std::string value;

				int eid, idx = 0;
				vector<int> nids;
				while (getline(linestream, value, ' ')) {
					//cout << value << endl;
					if (idx == 0) {
						eid = stoi(value);
						idx++;
					} else {
						nids.push_back(stoi(value));
					}
				}
				elements.insert(pair<int, vector<int> >(eid, nids));
				nels++;
				//throw;
			}
		}
	}
	infile.close();

	if (nels == 0) {
		cout << endl << "*** ELEMENT READING ERROR: " << folder + "/" + filename + ".node" << " CONTAINS ZERO ELEMENTS." << endl;
		cout << "*** CHECK STL OR MESHING PARAMETERS." << endl;
		exit(-1);
	}
	return 0;
}

int TetGen::parseSurfNodes(vector<int>& surfaceNodes, map<int, vector<int> >& facetMap, int facetID, string line) {
	line = cleanString(line + " ");
	std::stringstream linestream(line);
	std::string value;
	int idx = 0;
	int nid;
	vector<int> facetDef;

	while (getline(linestream, value, ' ')) {
		if (idx == 0) {
			idx++;
		} else {
			nid = stoi(value);
			if (std::find(surfaceNodes.begin(), surfaceNodes.end(), nid)
					!= surfaceNodes.end()) {
				// Do Nothing
			} else {
				surfaceNodes.push_back(nid);
			}

			facetDef.push_back(nid);
		}

	}

	facetMap.insert(pair<int, vector<int> >(facetID, facetDef));
	return 0;
}

TetGen::~TetGen() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
