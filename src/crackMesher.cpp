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


#include "debug.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp> // We need both -lboost_filesystem-mt -lboost_system-mt
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <unistd.h>

#include "GlobalMesh.h"
#include "Grain.h"

#include "MeshIO/asciiSTL.h"
#include "MeshIO/INPReader.h"

#include "refine/refineSurf.h"
#include "refine/coarsenSurf.h"

#include "programSettings.h"

using namespace std;
namespace fs = ::boost::filesystem;
namespace po = boost::program_options;

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v) {
	copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
	return os;
}

void enumerateStlFiles(fs::path& folder, vector<fs::path>& out) {
	string ext = ".stl";
	if (!fs::exists(folder) || !fs::is_directory(folder))
		return;

	fs::recursive_directory_iterator it(folder);
	fs::recursive_directory_iterator endit;

	while (it != endit) {
		//cout << it->path().extension() << endl;
		if (fs::is_regular_file(*it) && it->path().extension() == ext) {
			out.push_back(it->path().filename());
		}
		++it;

	}
	// Sort the input mesh filenames to ensure consistent order on different machines
	std::sort(out.begin(), out.end());
}

void enumerateStlFiles2(fs::path& folder, vector<fs::path>& out) {
	string ext = ".stl";
	if (!fs::exists(folder) || !fs::is_directory(folder))
		return;

	fs::recursive_directory_iterator it(folder);
	fs::recursive_directory_iterator endit;

	while (it != endit) {
		//cout << it->path().filename().string() << endl;
		if ((fs::is_regular_file(*it) && it->path().extension() == ext)
				&& (boost::starts_with(it->path().filename().string(),
						"indgrain_") == true)) {
			out.push_back(it->path().filename());
			//cout << "Detected :" << it->path().filename().string() << endl;
		}
		++it;

	}

}

int main(int ac, char* av[]) {
	string stlPath, cont_string, ref_string, advFilePath, settingsPath;

	int binary, smesh_ref, smesh_coarsen = 0;
	try {
		po::options_description desc("Allowed options");
		desc.add_options()("input-file", po::value<vector<string> >(),
				"input file");
		po::positional_options_description p;
		p.add("input-file", -1);

		po::variables_map vm;
		po::store(
				po::command_line_parser(ac, av).options(desc).positional(p).run(),
				vm);
		po::notify(vm);

		if (vm.count("input-file")) {
			settingsPath = vm["input-file"].as<vector<string> >()[0];
		} else {
			cerr << "Input file not specified." << endl;
			return 1;
		}

	} catch (std::exception& e) {
		cout << "Failed at processing input file: " << e.what() << "\n";
		return 1;
	}

	programSettings settings;
	try {
		settings.loadData(settingsPath);
		stlPath = settings.stlsPath;
		cout << "Stl-path is: " << stlPath << "\n";
		binary = settings.ascii2binaryConversion;
		if (binary == 1) {
			cout << "Binary to ascii enabled." << "\n";
		}
		smesh_ref = settings.refinementPasses;
		if (smesh_ref > 0) {
			cout << "Surface mesh refinement enabled: " << smesh_ref
					<< " passes.\n";
		}
		smesh_coarsen = settings.coarseningPasses;
		if (smesh_coarsen > 0) {
			cout << "Surface mesh coarsening enabled: " << smesh_coarsen
					<< " passes.\n";
		}

		cont_string = settings.continuumQuality;
		cout << "Tetgen quality argument: " << cont_string << "\n";

		advFilePath = settings.advancingCrackPath;
		ref_string = cont_string;
		cout << "File path for AdvancingCrack: " << advFilePath << "\n";

	} catch (std::exception& e) {
		cout << "Failed at reading input file: " << e.what() << "\n";
		cout << "Check that you have everything required to run crackMesher."
				<< endl;
		return 1;
	}


	vector<fs::path> stlFiles;
	fs::path folder = stlPath;

	enumerateStlFiles(folder, stlFiles);

	vector<Grain> grains;

	cout << ">>> Reading stl files... " << endl;
	double progress = 0.0;
	int barWidth = 30;
	// If we decide to write a summary binary file, write the number of grains
	bool binaryExists = false;
	if (settings.writeBinarySave.length() > 0) {
		// Check to see if this binary file is already written
		if ( !boost::filesystem::exists( settings.writeBinarySave ) )
		{
			cout << ">>> Could not find specified binary write file. Writing out new binary file." << endl;
			// Write the number of grains in the binary file
			std::ofstream out(settings.writeBinarySave, ios_base::out );
			int numberOfStlFiles = stlFiles.size();
			out.write(reinterpret_cast<char *>(&numberOfStlFiles), sizeof(numberOfStlFiles));
		} else {
			// If it's here, load the binary file instead:
			cout << ">>> Found binary file. Reading from binary file instead of stl files." << endl;
 			binaryExists = true;
		}

	}

	// Read in each individual grain's surface mesh and nodes
	if (binaryExists == false) {
		for (uint i = 0; i < stlFiles.size(); i++) {
			//cout << stlFiles[i] << endl;
			Grain test(stlFiles[i].string());
			if (binary == 1) {
				test.convertStl(folder.string());
				test.writeSmshAndRead(folder.string());
			}

			if (smesh_ref > 0 || smesh_coarsen > 0) {
				grains.push_back(test);
			}

			// If we want to write the binary file, do it now
			if (settings.writeBinarySave.length() > 0) {
				test.exportBinary(settings.writeBinarySave);
			}

			std::cout << "[";
			progress = double(i) / double(stlFiles.size());
			int pos = barWidth * progress;
			for (int p = 0; p < barWidth; ++p) {
				if (p < pos)
					std::cout << "=";
				else if (p == pos)
					std::cout << ">";
				else
					std::cout << " ";
			}
			std::cout << "] " << int(progress * 100.0) << " %\r";
			std::cout.flush();
		}
		cout << endl;
	} else {
		// Create a single binary reader
		std::ifstream in(settings.writeBinarySave, ios_base::binary | ios_base::in );
		// How many grains do we have in this file?
		int numberOfGrainsInBinaryFile;
		in.read(reinterpret_cast<char *>(&numberOfGrainsInBinaryFile), sizeof(numberOfGrainsInBinaryFile));
		// Start looping through that many files.
		for (uint i = 0; i < numberOfGrainsInBinaryFile; i++) {
			// Grab what we need to initialize Grain class
			int grainNumber;
			in.read(reinterpret_cast<char*>(&grainNumber), sizeof(grainNumber));
			int filenameSize;
			in.read(reinterpret_cast<char*>(&filenameSize), sizeof(filenameSize));
			string grainFilename; grainFilename.resize(filenameSize);
			in.read(&grainFilename[0], filenameSize);
			int prefixSize;
			in.read(reinterpret_cast<char*>(&prefixSize), sizeof(prefixSize));
			string grainPrefix; grainPrefix.resize(prefixSize);
			in.read(&grainPrefix[0], prefixSize);

			Grain test(grainFilename, grainNumber, grainPrefix);

			test.importBinary(in);

			if (smesh_ref > 0 || smesh_coarsen > 0) {
				grains.push_back(test);
			}

			std::cout << "[";
			progress = double(i) / double(stlFiles.size());
			int pos = barWidth * progress;
			for (int p = 0; p < barWidth; ++p) {
				if (p < pos)
					std::cout << "=";
				else if (p == pos)
					std::cout << ">";
				else
					std::cout << " ";
			}
			std::cout << "] " << int(progress * 100.0) << " %\r";
			std::cout.flush();
		}
		cout << endl;
		in.close();
	}
	cout << ">>> Read " << stlFiles.size() << " stl files." << endl;

	map<int, vector<vector<double> > > hole_list;
	if (settings.detectIslands == 1) {

		ifstream holefile("holes");
		std::string s, line;

		while (std::getline(holefile, line)) {
			std::istringstream iss(line); // string stream
			std::getline(iss, s, ',');
			int GID = stoi(s);

			vector<double> centroid;
			for (int zzz = 0; zzz < 3; zzz++) {
				std::getline(iss, s, ',');
				centroid.push_back(stod(s));
			}

			// This is hole in the associated grain ID
			std::map<int, vector<vector<double> > >::iterator it_assoc =
					hole_list.find(GID);
			if (it_assoc != hole_list.end()) {
				// It's in this set, just insert it
				hole_list.at(GID).push_back(centroid);
			} else {
				vector<vector<double> > temp;
				temp.push_back(centroid);
				hole_list.insert(
						pair<int, vector<vector<double> > >(GID, temp));
			}

		}

	}

	if (settings.defectType == 0) {
		cout << ">>> This model contains defects of type: Cracks" << endl;
	} else if (settings.defectType == 1) {
		cout << ">>> This model contains defects of type: Voids" << endl;
	} else {
		cout << ">>> Defect type not recognized!" << endl;
		exit(1);
	}

	vector<Grain> grains2;
	if ((smesh_coarsen > 0) || (smesh_ref > 0)) {
		GlobalMesh mesh_surface(&settings);
		mesh_surface.refine_iterations = smesh_ref;

		// Build KD tree for gradation function

		cout << ">>> Undoing shrinking with a factor of "
				<< mesh_surface.Undo_shrink_factorX << ","  << mesh_surface.Undo_shrink_factorY << "," << mesh_surface.Undo_shrink_factorZ << endl;

		mesh_surface.BuildCrackFrontKD(advFilePath);
	

		for (uint i = 0; i < grains.size(); i++) {
			cout << ">>> Current kd_node_id: ";
			int current_kd_id = mesh_surface.processGrainSurfaceNodes(
					grains[i]);
			cout << current_kd_id << " \r";
			cout.flush();
		}
		cout << endl;

		for (uint i = 0; i < grains.size(); i++) {
			cout << ">>> Current max_facet_id: ";
			int max_facet_id = mesh_surface.processGrainFacets(grains[i], i);
			cout << max_facet_id << " \r";
			cout.flush();
		}
		cout << endl;

		asciiSTL astl;
		// astl.writeGlobalFacets(mesh_surface, "original.stl");

		// First we do coarsening
		if (smesh_coarsen > 0) {
			for (int rr = 0; rr < smesh_coarsen; rr++) {
				cout << "\n\n&&& Coarsening surface meshes: pass " << rr + 1 << "."
						<< endl;
				mesh_surface.coarsen(rr);
				string astlname = "coarsen" + to_string(rr + 1) + ".stl";
				astl.writeGlobalFacets(mesh_surface, astlname);
			}
		}

		// Then we do refinement
		if (smesh_ref > 0) {
			for (int rr = 0; rr < smesh_ref; rr++) {
				cout << "\n\n&&& Refining surface meshes: pass " << rr + 1 << "."
						<< endl;
				mesh_surface.refine();
				string astlname = "split" + to_string(rr + 1) + ".stl";
				astl.writeGlobalFacets(mesh_surface, astlname);
			}
		}

		// Write the final, modified surface mesh 
		astl.writeGlobalFacets(mesh_surface, "modified.stl");

		// Now we have the surface mesh, let's separate the surface mesh such that we retain the grain ID
		asciiSTL writer;
		writer.writeIndividualFacets(folder.string(), mesh_surface);

		vector<fs::path> stlFiles2;
		enumerateStlFiles2(folder, stlFiles2);
		

		int Mode = 2;
		// Write a background volume mesh and calculate desired edge lengths
		progress = 0.0;
		cout << ">>> Writing background mesh and mtr files." << endl;
		for (uint i = 0; i < stlFiles2.size(); i++) {
			std::cout << "[";
			Grain test(stlFiles2[i].string());

			progress = double(i) / double(stlFiles2.size());
			int pos = barWidth * progress;
			for (int p = 0; p < barWidth; ++p) {
				if (p < pos)
					std::cout << "=";
				else if (p == pos)
					std::cout << ">";
				else
					std::cout << " ";
			}
			std::cout << "] " << int(progress * 100.0) << " % (GID: "
					<< test.giveRawGrainNumber() << ")               \r";
			std::cout.flush();

			int success = test.writeBackGroundSmsh(folder.string(), settings);
			if (success != 0) {
				cout << "ERROR: Background mesh generation failed! Please check the surface mesh to make sure it is volume-meshable!"
						<< endl;
				exit(1);
			}

			// Write mtr files
			mesh_surface.writeMTR(stlFiles2[i].string());

		}
		cout << endl;
		cout << ">>> Done writing background mesh and mtr files." << endl;

		// Can parallelize
		progress = 0.0;
		cout << ">>> Attempting to volume mesh..." << endl;
		for (uint i = 0; i < stlFiles2.size(); i++) {
			std::cout << "[";
			Grain test(stlFiles2[i].string());

			progress = double(i) / double(stlFiles2.size());
			int pos = barWidth * progress;
			for (int p = 0; p < barWidth; ++p) {
				if (p < pos)
					std::cout << "=";
				else if (p == pos)
					std::cout << ">";
				else
					std::cout << " ";
			}
			std::cout << "] " << int(progress * 100.0) << " % (GID: "
					<< test.giveRawGrainNumber() << ")               \r";
			std::cout.flush();

			bool success = -1;
			if (settings.detectIslands == 1) {
				success = test.meshGrain(folder.string(), cont_string,
						ref_string, Mode, settings, settings.detectIslands,
						hole_list);
			} else {
				success = test.meshGrain(folder.string(), cont_string,
						ref_string, Mode, settings);
			}

			if (success == 0) {
				test.grabSurfaceNodes(folder.string(), Mode);
				grains2.push_back(test);
			}

		}


	} 
	cout << ">>> Volume meshing complete!" << endl;
	GlobalMesh mesh(&settings);

	for (uint i = 0; i < grains2.size(); i++) {
		cout << ">>> Current kd_node_id: ";
		int current_kd_id = mesh.processGrainSurfaceNodes(grains2[i]);
		cout << current_kd_id << " \r";
		cout.flush();
	}
	cout << endl;
	// After Going through each Surface Mesh, run through each grain
	for (uint i = 0; i < grains2.size(); i++) {
		int current_nid, current_eid;
		mesh.processGrainElements(grains2[i], current_nid, current_eid);
		cout << ">>> Current max node id: " << current_nid
				<< ", current max_eid_id: " << current_eid << ".\r";
		cout.flush();
	}
	cout << endl;
	cout << ">>> Writing out Abaqus INP..." << endl;
	INPReader out;
	out.writeINP(mesh);

	cout << ">>> All DONE!" << endl;
	return 0;
}
