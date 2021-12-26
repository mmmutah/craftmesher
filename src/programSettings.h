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


#ifndef SRC_PROGRAMSETTINGS_H_
#define SRC_PROGRAMSETTINGS_H_

#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <iostream>
#include <vector>
#include <map>

namespace pt = boost::property_tree;

class programSettings {
public:
	std::string settingsPath, stlsPath, advancingCrackPath;
	int ascii2binaryConversion;
	std::string writeBinarySave;

	// Mesh Parameters
	double OriginalEdge;
	int gradationScheme;
	int defectType;

	std::string continuumQuality, refineQuality, backgroundQuality;

	// Refinement Parameters
	int refinementPasses;
	double ref_targetSize, ref_volTargetSize, R0, Rt, shrinkFactorX, shrinkFactorY, shrinkFactorZ;
	double powerExp;
	std::string refMethod;

	// Coarsening Parameters
	int coarseningPasses;
	double coar_targetSize, coar_volTargetSize;
	int RemoveHourglass;
	int SmoothBoundary;
	int Condense3;
	int EdgeSmoothing;

	// Detect Islands
	int detectIslands;


	// Load the Data
	int loadData(std::string configPath) {
		// Create empty property tree object
		pt::ptree tree;

		// Parse the XML into the property tree.
		settingsPath = configPath;
		pt::read_xml(configPath, tree);

		stlsPath = tree.get<std::string>("crackMesher.stls");

		ascii2binaryConversion = tree.get("crackMesher.binary", 0);

		writeBinarySave = tree.get<std::string>("crackMesher.writeBinarySave");

		// Mesh Parameters
		OriginalEdge = tree.get<double>("crackMesher.OriginalEdgeLength", 0);
		continuumQuality = tree.get<std::string>("crackMesher.TetGenQualityMeasure");
		backgroundQuality = tree.get<std::string>("crackMesher.BackgroundSMeshQuality");

		gradationScheme = tree.get("crackMesher.gradationScheme", 0);

		defectType = tree.get("crackMesher.defectType", 0);

		R0 = tree.get<double>("crackMesher.R0", 0);
		Rt = tree.get<double>("crackMesher.Rt", 0);
		powerExp = tree.get<double>("crackMesher.a", 3);

		advancingCrackPath = tree.get<std::string>("crackMesher.advancingCrackFile");
		shrinkFactorX = tree.get<double>("crackMesher.advCrackFileShrinkX");
		shrinkFactorY = tree.get<double>("crackMesher.advCrackFileShrinkY");
		shrinkFactorZ = tree.get<double>("crackMesher.advCrackFileShrinkZ");

		// Refinement Parameters
		refinementPasses = tree.get("crackMesher.refinement.refinement_passes", 0);
		ref_targetSize = tree.get<double>("crackMesher.refinement.SurfaceTargetEdgeLength", 0);
		ref_volTargetSize = tree.get<double>("crackMesher.refinement.VolumeTargetEdgeLength", 0);

		// Coarsening Parameters
		coarseningPasses = tree.get("crackMesher.coarsening.coarsening_passes", 0);
		coar_targetSize = tree.get<double>("crackMesher.coarsening.SurfaceTargetEdgeLength", 0);
		coar_volTargetSize = tree.get<double>("crackMesher.coarsening.VolumeTargetEdgeLength", 0);

		// Detect Islands?
		detectIslands = tree.get("crackMesher.detectIslands", 0);

		RemoveHourglass = tree.get("crackMesher.coarsening.RemoveHourglass", 1);
		SmoothBoundary = tree.get("crackMesher.coarsening.SmoothBoundary", 1);
		Condense3 = tree.get("crackMesher.coarsening.Condense3", 1);
		EdgeSmoothing = tree.get("crackMesher.coarsening.EdgeSmoothing", 2);

		return 0;
	}

private:

};

#endif /* SRC_PROGRAMSETTINGS_H_ */
