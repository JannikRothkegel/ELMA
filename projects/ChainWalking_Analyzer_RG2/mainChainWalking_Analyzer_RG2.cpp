/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2018,2021 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
    ooo                        | Ron Dockhorn
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#include <cstring>
#include <stdlib.h> //for atoi
#include <unistd.h> //for getopt

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/utility/DepthIteratorPredicates.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFileSubGroup.h>

#include "Analyzer_ChainWalking_RG2.h"


int main(int argc, char* argv[])
{
  try{
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	//typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureFixedMonomers, FeatureAttributes<>,FeatureExcludedVolumeSc<>) Features;
	//typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureFixedMonomers, FeatureAttributes<>) Features;
    typedef LOKI_TYPELIST_2(FeatureMoleculesIO, FeatureAttributes<>) Features;

	typedef ConfigureSystem<VectorInt3,Features, 8> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	std::string filename="test.bfm";
	//uint32_t number_of_monomers=100;
	//double probability = 1.0;

	long evalulation_time = 0;

	int skip = 0;

	int option_char(0);

	//read in options by getopt
	while ((option_char = getopt (argc, argv, "f:n:e:s:h"))  != EOF){
		switch (option_char)
		{
		case 'f':
			filename=optarg;
			break;
		case 'e': evalulation_time = atol(optarg);
				  break;
		case 's': skip = atol(optarg);
				  break;
		case 'h':
		default:
			std::cerr << "Usage: " << argv[0] << " [-f filename(=test.bfm)] [-e evaluation_time(=0)] [-s skipframe(=0)]\n";

			return 0;
		}
	}

	//seed the globally available random number generators
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    myIngredients.setName(filename);

    TaskManager taskmanager;
    taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(filename,myIngredients,UpdaterReadBfmFile<Ing>::READ_STEPWISE));

    taskmanager.addAnalyzer(new Analyzer_ChainWalking_RG2<Ing>(myIngredients, evalulation_time), (skip+1));

    taskmanager.initialize();
    taskmanager.run();
    taskmanager.cleanup();

	//TaskManager taskmanager;
	//taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
	//here you can choose to use MoveLocalBcc instead. Careful though: no real tests made yet
	//(other than for latticeOccupation, valid bonds, frozen monomers...)
	//taskmanager.addUpdater(new UpdaterSimpleSimulator<Ing,MoveLocalSc>(myIngredients,save_interval));

	//taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));
	
	//taskmanager.initialize();
	//taskmanager.run(max_mcs/save_interval);
	//taskmanager.cleanup();
	
	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

