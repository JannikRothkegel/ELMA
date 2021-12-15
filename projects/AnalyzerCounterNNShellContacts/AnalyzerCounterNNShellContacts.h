/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2018,2021 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (Ron Dockhorn)
    ooo                        |
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

#ifndef AnalyzerCounterNNShellContacts_H
#define AnalyzerCounterNNShellContacts_H

#include <vector>
#include <string>
#include <utility> // std::pair
#include <map>
#include <vector>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/utility/DistanceCalculation.h>

#include "StatisticMoment.h"

// polymere attribute tag = 1
// solvent attribute tag = 2
// cosolvemt attribute tag = 3

template <class IngredientsType>
class AnalyzerCounterNNShellContacts : public AbstractAnalyzer
{
public:
	AnalyzerCounterNNShellContacts(const IngredientsType &ing, uint64_t startTime_, std::string dstDir_);

	virtual ~AnalyzerCounterNNShellContacts(){

	};

	//typedef typename IngredientsType::molecules_type molecules_type;
	const typename IngredientsType::molecules_type &molecules;

	const IngredientsType &getIngredients() const { return ingredients; }

	virtual void initialize();
	virtual bool execute();
	virtual void cleanup();

	

	void getNumberCoSolventInNNShell();

private:
	const IngredientsType &ingredients;

	StatisticMoment Statistic_numCosolventInShell;
	StatisticMoment Statistic_numCosolventAsBridge;

	StatisticMoment Statistic_NumCosolventMonomers;

	uint64_t startTime;

	std::string filename;
	std::string dstdir;

	
};

/////////////////////////////////////////////////////////////////////////////

template <class IngredientsType>
AnalyzerCounterNNShellContacts<IngredientsType>::AnalyzerCounterNNShellContacts(const IngredientsType &ing, uint64_t startTime_, std::string dstDir_)
	: ingredients(ing), molecules(ing.getMolecules()), startTime(startTime_), dstdir(dstDir_)
{

	Statistic_numCosolventInShell.clear();
	Statistic_numCosolventAsBridge.clear();
	Statistic_NumCosolventMonomers.clear();
}

template <class IngredientsType>
void AnalyzerCounterNNShellContacts<IngredientsType>::initialize()
{

	//execute();
}

template <class IngredientsType>
bool AnalyzerCounterNNShellContacts<IngredientsType>::execute()
{
	// time of the conformations
	//uint64_t timeInSim =  ingredients.getMolecules().getAge();
	//molecules <-> ingredients.getMolecules()

	if (ingredients.getMolecules().getAge() >= startTime)
	{
		std::cout << "AnalyzerCounterNNShellContacts.execute() at MCS:" << ingredients.getMolecules().getAge() << std::endl;

		

		getNumberCoSolventInNNShell();

		
	}

	return true;
}

template <class IngredientsType>
void AnalyzerCounterNNShellContacts<IngredientsType>::getNumberCoSolventInNNShell() 
{

	int32_t numberOfCosolventInShell = 0;
	int32_t numberOfCosolventAsBridges = 0;

	//loop through all polymer
	for (size_t n = 0; n < ingredients.getMolecules().size(); n++)
	{
		
		int32_t monoType = ingredients.getMolecules()[n].getAttributeTag();

		VectorInt3 posOfCosolvent = ingredients.getMolecules()[n];

		// only for cosolvent
		if (monoType == 3)
		{
			int32_t counterContacts = 0;

			for (size_t m = 0; m < ingredients.getMolecules().size(); m++)
			{
				int32_t monoType = ingredients.getMolecules()[m].getAttributeTag();

				// only for polymerechain
				if (monoType == 1)
				{
					VectorInt3 posOfChain = ingredients.getMolecules()[m];
					VectorInt3 diff = posOfCosolvent - posOfChain;

					diff.setX(LemonadeDistCalcs::fold(diff.getX(), ingredients.getBoxX()));
					diff.setY(LemonadeDistCalcs::fold(diff.getY(), ingredients.getBoxY()));
					diff.setZ(LemonadeDistCalcs::fold(diff.getZ(), ingredients.getBoxZ()));

					if (diff * diff <= 6)
					{
						counterContacts++;
					}
				}
			}
			// add one to statistic if at least one interaction is detected
			if (counterContacts > 0)
			{
				numberOfCosolventInShell++;
			}
			// add one to bridge statistic if more thanone interactioni detected
			if (counterContacts > 1)
			{
				numberOfCosolventAsBridges++;
			}



		}

	}
	Statistic_numCosolventInShell.AddValue(numberOfCosolventInShell);
	Statistic_numCosolventAsBridge.AddValue(numberOfCosolventAsBridges);


	return;
}

struct PathSeparator
{
	bool operator()(char ch) const
	{
		return ch == '\\' || ch == '/';
	}
};

template <class IngredientsType>
void AnalyzerCounterNNShellContacts<IngredientsType>::cleanup()
{	

	int32_t counterCosolvent = 0;

	for (size_t n = 0; n < ingredients.getMolecules().size(); n++)
	{
		int32_t monoType = ingredients.getMolecules()[n].getAttributeTag();
		
		// only for cosolvent
		if (monoType == 3)
		{
			counterCosolvent++;
		}
	}

	std::cout << "File output" << std::endl;
	

	std::vector<std::vector<double>> tmpResults;

	tmpResults.resize(5);

	for (int i = 0; i < 5; i++)
		tmpResults[i].resize(1);

	tmpResults[0][0] = counterCosolvent;
	tmpResults[1][0] = Statistic_numCosolventInShell.ReturnM1();
	tmpResults[2][0] = Statistic_numCosolventInShell.ReturnM2();
	tmpResults[3][0] = Statistic_numCosolventAsBridge.ReturnM1() ;
	tmpResults[4][0] = Statistic_numCosolventAsBridge.ReturnM2();

	std::stringstream comment;
	comment << "File produced by analyzer AnalyzerCounterNNShellContacts\n"
			<< "Analyze CoSolventPolyereBridges\n"
			<< "eta: Number of Cosolvent in NNShell in vicinity of polymere\n"
			<< "gamma: Number of Bridge building Cosolvent with Polymere\n"
			<< "\n"
			<< "numCoSolvent\t<eta>\t<eta²>\t<gamma>\t<gamma^2>\n";

	// find the filename without path and extensions
	std::string filenameGeneral = std::string(std::find_if(ingredients.getName().rbegin(), ingredients.getName().rend(), PathSeparator()).base(), ingredients.getName().end());

	std::string::size_type const p(filenameGeneral.find_last_of('.'));
	filenameGeneral = filenameGeneral.substr(0, p);

	std::string filenameRg2_Ree_b2 = filenameGeneral + "_AnalyzerCounterNNShellContacts.dat";

	std::cout << " Write output to: " << dstdir << "/" << filenameRg2_Ree_b2 << std::endl;

	ResultFormattingTools::writeResultFile(dstdir + "/" + filenameRg2_Ree_b2, this->ingredients, tmpResults, comment.str());
}

#endif /*AnalyzerCounterNNShellContacts_H*/
