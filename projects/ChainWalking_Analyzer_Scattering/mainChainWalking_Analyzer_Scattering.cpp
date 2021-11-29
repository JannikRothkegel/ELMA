

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

#include "Analyzer_ChainWalking_Scattering.h"


int main(int argc, char* argv[])
{
  try{
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	//typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureFixedMonomers, FeatureAttributes,FeatureExcludedVolumeSc<>) Features;
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureFixedMonomers, FeatureAttributes<>) Features;

	typedef ConfigureSystem<VectorInt3,Features, 8> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	std::string filename="test.bfm";
	//uint32_t number_of_monomers=100;
	//double probability = 1.0;

	long evalulation_time = 0;

	int option_char(0);

	//read in options by getopt
	while ((option_char = getopt (argc, argv, "f:n:e:h"))  != EOF){
		switch (option_char)
		{
		case 'f':
			filename=optarg;
			break;

		case 'e': evalulation_time = atol(optarg);
				  break;
		//case 'n':
		//	number_of_monomers = atoi(optarg);
		//	break;
		//case 'p':
		//	probability = atof(optarg);
		//	break;
		case 'h':
		default:
			//std::cerr << "Usage: " << argv[0] << " [-f filename] [-n number_of_monomers] [-p probability] \n";
			std::cerr << "Usage: " << argv[0] << " [-f filename(=test.bfm)] [-e evaluation_time(=0)] \n";

			return 0;
		}
	}

	//seed the globally available random number generators
    RandomNumberGenerators randomNumbers;
    randomNumbers.seedAll();

    myIngredients.setName(filename);

    TaskManager taskmanager;
    taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(filename,myIngredients,UpdaterReadBfmFile<Ing>::READ_STEPWISE));

    taskmanager.addAnalyzer(new Analyzer_ChainWalking_Scattering<Ing>(myIngredients, evalulation_time));

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

