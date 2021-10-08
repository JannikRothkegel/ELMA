

#include <cstring>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/updater/UpdaterAddLinearChains.h>

#include "catchorg/clara/clara.hpp"

int main(int argc, char* argv[])
{
  try{
	std::string infile  = "input.bfm";
	std::string outfile = "outfile.bfm";
	uint32_t newMolecules = 100;
	
	
    bool showHelp = false;
    
    auto parser
    = clara::Opt( infile, "input (=input.bfm)" )
        ["-i"]["--infile"]
        ("BFM-file to load.")
        .required()
    | clara::Opt( outfile, "output (=outfile.bfm)" )
        ["-o"]["--outfile"]
        ("BFM-file to save.")
        .required()
    | clara::Opt( [&newMolecules](int const m)
        {
         if (m < 0)
         {
            return clara::ParserResult::runtimeError("Number of new molecules must be greater or equal to 0.");
         }
         else
         {
            newMolecules = m;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }
        }, "newMolecules" )
        ["-m"]["--newMolecules"]
        ("(required) specifies the total number of new molecules.")
        .required()

    | clara::Help( showHelp );
        
    auto result = parser.parse( clara::Args( argc, argv ) );
    if( !result ) {
    std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
    exit(1);
    }
    else if(showHelp == true)
    {
        std::cout << "Simulator for the ScBFM with Ex.Vol and BondCheck" << std::endl
                  << "maximum number of connections per monomer is 6" << std::endl
                  << "Features used: FeatureBondset, FeatureAttributes, FeatureExcludedVolumeSc<FeatureLattice<bool> >" << std::endl
		          << "Updaters used: ReadFullBFMFile, SimpleSimulator" << std::endl
		          << "Analyzers used: WriteBfmFile" << std::endl;
        
        parser.writeToStream(std::cout);
        exit(0);
    }
    else
    {
        std::cout << "infile:        " << infile << std::endl
                  << "outfile:       " << outfile << std::endl
                  << "newMolecules:       " << newMolecules << std::endl;
                  
    }
       
	//seed the globally available random number generators
	RandomNumberGenerators rng;
	rng.seedAll();
	
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes< >,FeatureExcludedVolumeSc<>) Features;
	
	typedef ConfigureSystem<VectorInt3,Features, 6> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	TaskManager taskmanager;
	taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
	

    // here new Molecules are added by using UpdaterAddLinearChains with a lengh of 1 
	taskmanager.addUpdater(new UpdaterAddLinearChains<Ing>(myIngredients, newMolecules, 1,1,2, true));
    


	taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));
	
	taskmanager.initialize();
	
	taskmanager.cleanup();
	
	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

