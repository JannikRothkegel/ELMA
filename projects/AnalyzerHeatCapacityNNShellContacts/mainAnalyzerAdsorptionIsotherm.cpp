

#include <cstring>

#include <iostream>
#include <iomanip>

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

#include "catchorg/clara/clara.hpp"

#include "AnalyzerAdsorptionIsotherm.h"



int main(int argc, char* argv[])
{
	try{
		std::string infile  = "input.bfm";
		std::string outfile = "outfile.bfm";
		uint64_t max_mcs=100;
		uint32_t save_interval=100;
		uint32_t numCoSolvent = 0;
		double nn_interation = -0.8;
		
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
		| clara::Opt( [&max_mcs](uint64_t const m)
					   {
			if (m <= 0)
			{
				return clara::ParserResult::runtimeError("Simulation time must be greater than 0");
			}
			else
			{
				max_mcs = m;
				return clara::ParserResult::ok(clara::ParseResultType::Matched);
			}
					   }, "max MCS(=100)" )
		["-m"]["--max-mcs"]
			   ("(required) specifies the total Monte-Carlo steps to simulate.")
			   .required()
			   | clara::Opt( [&save_interval](int const s)
					   {
			if (s < 0)
			{
				return clara::ParserResult::runtimeError("Save intervall must be greater than 0");
			}
			else
			{
				save_interval = s;
				return clara::ParserResult::ok(clara::ParseResultType::Matched);
			}
					   }, "save MCS(=100)")
		["-s"]["--save-mcs"]
			   ("(required) Save after every <integer> Monte-Carlo steps to the output file." )
			   .required()
		| clara::Opt( [&numCoSolvent](int const n)
					   {
			if (n < 0)
			{
				return clara::ParserResult::runtimeError("Number of CoSolvent must be greater than 0");
			}
			else
			{
				numCoSolvent = n;
				return clara::ParserResult::ok(clara::ParseResultType::Matched);
			}
					   }, "(=0)")
		["-n"]["--number-cosolvent"]
			   ("(required) Number of CoSolvent for interacting with polymer." )
			   .required()
		| clara::Opt(  nn_interation, "(=-0.8)" )
		["-e"]["--epsilon"]
		 ("NNShell interaction (=-0.8)")
		 .required()
		 | clara::Help( showHelp );

		auto result = parser.parse( clara::Args( argc, argv ) );
		if( !result ) {
			std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
			exit(1);
		}
		else if(showHelp == true)
		{
			std::cout << "Simulator for the ScBFM with Ex.Vol and BondCheck in NN-shell" << std::endl
					<< "maximum number of connections per monomer is 2" << std::endl;

			parser.writeToStream(std::cout);
			exit(0);
		}
		else
		{
			std::cout << "infile:        " << infile << std::endl
					<< "outfile:       " << outfile << std::endl
					<< "max_mcs:       " << max_mcs << std::endl
					<< "save_interval: " << save_interval << std::endl
					<< "numCoSolvent: " << numCoSolvent << std::endl
					<< "nn-interation: " << nn_interation << std::endl

					;
		}
	
	//seed the globally available random number generators
	RandomNumberGenerators rng;
	rng.seedAll();
	
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	//typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureFixedMonomers,FeatureAttributes,FeatureExcludedVolumeSc<>) Features;
	//typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes, FeatureNNInteractionSc< FeatureLattice >) Features;
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes< >, FeatureNNInteractionSc< FeatureLattice >) Features;

	typedef ConfigureSystem<VectorInt3,Features, 2> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	TaskManager taskmanager;
	taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile, myIngredients,UpdaterReadBfmFile<Ing>::READ_STEPWISE));
	//taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
	//here you can choose to use MoveLocalBcc instead. Careful though: no real tests made yet
	//(other than for latticeOccupation, valid bonds, frozen monomers...)
	//taskmanager.addUpdater(new UpdaterSimpleSimulator<Ing,MoveLocalSc>(myIngredients,save_interval));
    
    taskmanager.addAnalyzer(new AnalyzerAdsorptionIsotherm<Ing>(myIngredients, 50000,  "./", 1024));

	//taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));
	
	taskmanager.initialize();
	taskmanager.run();
	//taskmanager.run(max_mcs/save_interval);
	taskmanager.cleanup();
	
    /*
	TaskManager taskmanager;
	taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile, ingredients,UpdaterReadBfmFile<Ing>::READ_STEPWISE));

	taskmanager.addAnalyzer(new AnalyzerAdsorptionIsotherm<Ing>(ingredients, 0,  "./", 1024));

	taskmanager.initialize();
	taskmanager.run();
	taskmanager.cleanup();
	*/
	
	// Export everything to a new file (overwrites maybe an old one)
	//AnalyzerWriteBfmFile<Ing> AW_BFM( ingredients.getName(),ingredients, AnalyzerWriteBfmFile<Ing>::OVERWRITE);
	//AW_BFM.initialize();
	//AW_BFM.execute();
	//AW_BFM.cleanup();

	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

