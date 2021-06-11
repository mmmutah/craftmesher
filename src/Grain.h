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


#ifndef GRAIN_H_
#define GRAIN_H_

#include "debug.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "MeshIO/TetGen.h"
#include "programSettings.h"

namespace std {

class Grain {
public:

	// Class Data
	// TODO: move to private
	vector<int> surfaceNodes;
	map<int, vector<int> > localElements;
	map<int, vector<double> > localNodes;
	map<int, int> kdTreeMap;
	map<int, vector<int> > localFacets;
	map<int, int> grainToGlobalFacets;
	map<int, int> globalToGrainFacets;

	Grain(string fname);
	Grain(string stlFilename, int grainNumber, string stlPrefix);

	int meshGrain(string folder, string cont_str, string ref_str, int Mode, programSettings &settings );
	int meshGrain(string folder, string cont_str, string ref_str, int Mode, programSettings &settings, int &detectIsland, map<int, vector<vector<double> > > &hole_list );
	
	int convertStl(string folder);
	int writeSmshAndRead(string folder);
	int writeBackGroundSmsh(string folder, programSettings &settings);
	void writeStl();

	int importBinary(ifstream &in);
	int exportBinary(string filename);

	virtual ~Grain();

	int getKDMap(int localNID, int &outputNID);
	int inInternalNodeMap(int localNID, int &outputNID, int &global_nid_count);

	int grabSurfaceNodes(string folder, int Mode);

	int givePositiveGrainNumber() {
		if (grainNumber < 0) {
			return abs((grainNumber+1));
		} else {
			return abs(grainNumber);
		}

	}

	int giveBaseGrainNumber() {
		if (grainNumber < 0) {
			return abs(((grainNumber+1))) % 100000  ;
		} else {
			return abs(grainNumber) % 100000;
		}

	}

	int giveRawGrainNumber() {
		if (grainNumber < 0) {
			return grainNumber + 1 ;
		}
		return grainNumber;
	}


private:
	string filename;

	int grainNumber;
	string prefix;
	map<int, int> internalNodeMap;
	vector<vector<double> > localPoints;
	
	
	int facetCount = 1;


};

} /* namespace std */

#endif /* GRAIN_H_ */
