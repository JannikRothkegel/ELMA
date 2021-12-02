

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

//#include "Updater_RelabelingLinearChain.h"
#include <LeMonADE/feature/FeatureNNInteractionReadWrite.h>

int main(int argc, char* argv[])
{
  try{
	//std::string infile  = "input.bfm";
	std::string outfile = "outfile.bfm";
	uint32_t CoSolvent_1 = 0;
	uint32_t CoSolvent_2 = 0;

    uint32_t ChainLength = 128;

    int32_t AttributCoSolvent_1 = 2;
	int32_t AttributCoSolvent_2 = 3;

    bool mark_CS1_asSolvent = false;
    bool mark_CS2_asSolvent = false;

    int32_t StandartBoxSize_X = 32;
    int32_t StandartBoxSize_Y = 32;
    int32_t StandartBoxSize_Z = 32;

    

    // Interaction Energy in kT between the 3 types of monomers
    double Interaction_1_2 = 0.0;
    double Interaction_1_3 = 0.0;
    double Interaction_2_3 = 0.0;

    bool is_negative_1_2 = false;
    bool is_negative_1_3 = false;
    bool is_negative_2_3 = false;

    bool showHelp = false;
    
    auto parser
    = clara::Opt( outfile, "output (=outfile.bfm)" )
        ["-o"]["--outfile"]
        ("BFM-file to save.")
        .required()

    | clara::Opt( [&CoSolvent_1](int const cs_1)
        {
         if (cs_1 < 0)
         {
            return clara::ParserResult::runtimeError("Number of CoSolvent_1 must be greater or equal to 0.");
         }
         else
         {
            CoSolvent_1 = cs_1;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }
        }, "CoSolvent_1" )
        ["-c"]["--CoSolvent_1"]
        ("(required) specifies the total number of CoSolvent_1 molecules.")
        .required()

    | clara::Opt( [&ChainLength](int const cl)
        {
         if (cl < 0)
         {
            return clara::ParserResult::runtimeError("Chain length must not be greater than 0.");
         }
         else
         {
            ChainLength = cl;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }
        }, "ChainLength" )
        ["-l"]["--ChainLength"]
        ("(required) specifies the total number of Chain elements. ")
        .required()

    | clara::Opt( [&CoSolvent_2](int const cs_2)
        {
         if (cs_2 < 0)
         {
            return clara::ParserResult::runtimeError("Number of CoSolvent_2 must be greater or equal to 0.");
         }
         else
         {
            CoSolvent_2 = cs_2;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }
        }, "CoSolvent_2" )
        ["-s"]["--CoSolvent_2"]
        ("(required) specifies the total number of CoSolvent_2 molecules.")
        .required()

    | clara::Opt( [&Interaction_1_2](double const int_1_2)
        {
            Interaction_1_2 = int_1_2;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }, "Interaction1_2" )
        ["-a"]["--Interation1_2"]
        ("Specifies the energy bewteen Polymere (1) and Solvent (2).")
        .required()

    | clara::Opt( [&Interaction_1_3](double const int_1_3)
        {
            Interaction_1_3 = int_1_3;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }, "Interaction1_3" )
        ["-b"]["--Interation1_3"]
        ("Specifies the energy bewteen Polymere (1) and Cosolvent (3).")
        .required()

    | clara::Opt( [&Interaction_2_3](double const int_2_3)
        {
            Interaction_2_3 = int_2_3;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }, "Interaction2_3" )
        ["-c"]["--Interation2_3"]
        ("Specifies the energy bewteen Solvent (2) and Cosolvent (3).")
        .required()

    | clara::Opt( [&is_negative_1_2](bool const is_n_1_2)
        {
            is_negative_1_2 = is_n_1_2;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }, "is_negative_1_2" )
        ["-x"]["--is_negative_1_2"]
        ("Specifies wether the energy bewteen Polymere (1) and Solvent (2) is negative.")
        .required()

    | clara::Opt( [&is_negative_1_3](bool const is_n_1_3)
        {
            is_negative_1_3 = is_n_1_3;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }, "is_negative_1_3" )
        ["-y"]["--is_negative_1_3"]
        ("Specifies wether the energy bewteen Polymere (1) and Cosolvent (3) is negative.")
        .required()

    | clara::Opt( [&is_negative_2_3](bool const is_n_2_3)
        {
            is_negative_2_3 = is_n_2_3;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }, "is_negative_2_3" )
        ["-z"]["--is_negative_2_3"]
        ("Specifies wether the energy bewteen Solvent (2) and Cosolvent (3) is negative.")
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
                  << "Features used: FeatureMoleculesIO, FeatureAttributes< >,FeatureNNInteractionSc< FeatureLattice >" << std::endl
		          << "Updaters used: ReadFullBFMFile, SimpleSimulator,AddLinearChains" << std::endl
		          << "Analyzers used: WriteBfmFile" << std::endl;
        
        parser.writeToStream(std::cout);
        exit(0);
    }
    else
    {
         
        std::cout << "outfile:       " << outfile << std::endl
                  << "Number of CS 1:       " << CoSolvent_1 << std::endl
                  << "Number of CS 2:       " << CoSolvent_2 << std::endl
                  << "Chainlength:       " << ChainLength << std::endl
                  << "Interaction between Polymere and Solvent in kT:       " << Interaction_1_2 << std::endl
                  << "Interaction between Polymere and CoSolvent in kT:       " << Interaction_1_3 << std::endl
                  << "Interaction between Solvent and CoSolvent in kT:       " << Interaction_2_3 << std::endl;
                  
    }

    

	//seed the globally available random number generators
	RandomNumberGenerators rng;
	rng.seedAll();
	
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	// typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes< >,FeatureExcludedVolumeSc<>) Features;
    typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes< >,FeatureNNInteractionSc< FeatureLattice >) Features;

	
	typedef ConfigureSystem<VectorInt3,Features, 6> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

    myIngredients.setBoxX(StandartBoxSize_X);
    myIngredients.setBoxY(StandartBoxSize_Y);
    myIngredients.setBoxZ(StandartBoxSize_Z);

    if(is_negative_1_2){
        Interaction_1_2 *= -1;
    }

    if(is_negative_1_3){
        Interaction_1_3 *= -1;
    }

    if(is_negative_2_3){
        Interaction_2_3 *= -1;
    }


    myIngredients.setPeriodicX(true);
    myIngredients.setPeriodicY(true);
    myIngredients.setPeriodicZ(true);

    myIngredients.modifyBondset().addBFMclassicBondset();
    
    myIngredients.synchronize();

    myIngredients.setNNInteraction(1,2, Interaction_1_2);
    myIngredients.setNNInteraction(2,3, Interaction_2_3);
    myIngredients.setNNInteraction(1,3, Interaction_1_3);
    myIngredients.synchronize();

	TaskManager taskmanager;
	//taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
    taskmanager.addUpdater(new UpdaterAddLinearChains<Ing>(myIngredients, 1, ChainLength, 1, 1, false));

    // here CoSelvent_1 is added as a chain of lengh 1 with attribute 2
	taskmanager.addUpdater(new UpdaterAddLinearChains<Ing>(myIngredients, CoSolvent_1, 1, AttributCoSolvent_1, AttributCoSolvent_1, mark_CS1_asSolvent));

    // here CoSelvent_2 is added as a chain of lengh 1 with attribute 3
    taskmanager.addUpdater(new UpdaterAddLinearChains<Ing>(myIngredients, CoSolvent_2, 1, AttributCoSolvent_2, AttributCoSolvent_2, mark_CS2_asSolvent));

	taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));
	
	taskmanager.initialize();

    taskmanager.run(1);
	
	taskmanager.cleanup();
	
	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

