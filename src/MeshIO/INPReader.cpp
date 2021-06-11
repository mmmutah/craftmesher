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


#include "INPReader.h"

namespace std {

INPReader::INPReader() {
//	// Find the INP file in the directory
}

INPReader::~INPReader() {
	// TODO Auto-generated destructor stub
}

void INPReader::writeINP(GlobalMesh &mesh) {

	// Define our new INP name
	string inpout = "PPM_ABAQUS.inp";
	ofstream outfile(inpout);

	//std::string line;
	int mode = 0;

	outfile << "*Heading" << endl;
	outfile << "** nothing to see here" << endl;
	outfile << "*Part, name=Micro" << endl;
	outfile << "*Node" << endl;

	cout << ">>> Writing new nodes to INP file: " << inpout << endl;
	// For all the nodes,
	typedef map<int, vector<double> >::iterator it_type;
	for (it_type it = mesh.nodes.begin(); it != mesh.nodes.end(); it++) {
		int node_id = it->first;
		vector<double> xyz = it->second;
		//				cout << node_id << ", " << xyz[0] << ", " << xyz[1] << ", "
		//						<< xyz[2] << "\n";
		outfile << node_id << ", " << xyz[0] << ", " << xyz[1] << ", "
				<< xyz[2] << "\n";
	}
	cout << ">>> Completed writing new nodes to INP file: " << inpout
			<< endl;

	outfile << "*Element, type=C3D4" << endl;

	cout << ">>> Writing new elements to INP file: " << inpout << endl;
	// For all the elements
	typedef map<int, vector<int> >::iterator it_type2;
	for (it_type2 it = mesh.elements.begin(); it != mesh.elements.end(); it++) {
		int element_id = it->first;
		vector<int> n = it->second;
		outfile << element_id << ", " << n[0] << ", " << n[1] << ", "
				<< n[2] << ", " << n[3] << "\n";
	}
	cout << ">>> Completed writing new elements to INP file: " << inpout
			<< endl;

	cout << ">>> Writing elsets to INP file: " << inpout << endl;
	// For all the elements
	typedef map<int, vector<int> >::iterator it_type2;
	for (it_type2 it = mesh.elementSets.begin(); it != mesh.elementSets.end(); it++) {
		int elset = it->first;
		outfile << "*Elset, elset=Grain_" << elset << "\n";
		vector<int> els = it->second;
		for (uint e = 0; e < els.size(); e++) {
			outfile << els[e] << ",\n";
		}
	}
	cout << ">>> Completed writing new elsets to INP file: " << inpout
			<< endl;
	outfile << "*End Part" << endl;
	outfile << "* Assembly, name=Assem" << endl;
	outfile << "*Instance, name=Micro-1, part=Micro" << endl;
	outfile << "*End Instance" << endl;
	outfile << "*End Assembly" << endl;
}



} /* namespace std */
