/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2021 by
  o/.|.\o    E   nvironment    | Jannik
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

#ifndef LEMONADE_UPDATER_RELABEL_LINEAR_CHAINS
#define LEMONADE_UPDATER_RELABEL_LINEAR_CHAINS
/**
 * @file
 *
 * @class Updater_RelabelingLinearChain
 *
 * @brief Updater to create a solution of monomdisperse linear chains.
 *
 * @details This is a simple implementation of a system setup starting from an empty ingredients
 * or a system with some monomers inside. This updater requires FeatureAttributes.
 * Two tags are added to the monomers in alternating manner, usually needed for GPU computing.
 *
 * @tparam IngredientsType
 *
 * @param ingredients_ The system, holding eigther an empty simulation box for system setup
 * or a prefilled ingredients where the linear chains shall be added
 * @param NChain_ number of chains that are added to ingredients
 * @param NMonoPerChain_ number of monomer is each chain
 * @param type1_ attribute tag of "even" monomers
 * @param type2_ attribute tag of "odd" monomers
 **/

#include <LeMonADE/updater/UpdaterAbstractCreate.h>
#include <LeMonADE/utility/Vector3D.h>
#include <cmath>

template<class IngredientsType>
class Updater_RelabelingLinearChain: public UpdaterAbstractCreate<IngredientsType>
{
  typedef UpdaterAbstractCreate<IngredientsType> BaseClass;

public:
  Updater_RelabelingLinearChain(IngredientsType& ingredients_, uint32_t Size_, uint32_t Tag_);

  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();


private:
  // provide access to functions of UpdaterAbstractCreate used in this updater
  using BaseClass::ingredients;
  using BaseClass::addMonomerToParent;
  using BaseClass::addSingleMonomer;
  using BaseClass::linearizeSystem;

  //! bool for execution
  bool wasExecuted;


  //! number of linear chains in the box
  uint32_t Size;
  uint32_t Tag;

};

/**
* @brief Constructor handling the new systems paramters
*
* @param ingredients_ a reference to the IngredientsType - mainly the system
* @param NChain_ number of chains to be added in the system instead of solvent
* @param NMonoPerChain_ number of monomers in one chain
*/
template < class IngredientsType >
Updater_RelabelingLinearChain<IngredientsType>::Updater_RelabelingLinearChain(IngredientsType& ingredients_, uint32_t Size_, uint32_t Tag_):
BaseClass(ingredients_), Size(Size_), Tag(Tag_), wasExecuted(false)
{

}

/**
* @brief initialise function, calculate the target density to compare with at the end.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void Updater_RelabelingLinearChain<IngredientsType>::initialize(){
  

  execute();
}

/**
* @brief Execution of the system creation
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
bool Updater_RelabelingLinearChain<IngredientsType>::execute(){
    
    if(wasExecuted)
        return true;

    //ingredients.setBoxX(Size);
    //ingredients.setBoxY(Size);
    //ingredients.setBoxZ(Size);

    int32_t i = 0;
    while (i < (ingredients.getMolecules().size())){
        ingredients.modifyMolecules()[i].setAttributeTag(Tag);

        std::cout<<"Tag: " << i << "  " << ingredients.getMolecules()[i].getAttributeTag() <<std::endl;
        i++;
       
    }




    ingredients.synchronize();

    std::cout<<"Labling done"<<std::endl;




  
    return true;
  
}

/**
* @brief Standard clean up.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void Updater_RelabelingLinearChain<IngredientsType>::cleanup(){

}


#endif /* LEMONADE_UPDATER_RELABEL_LINEAR_CHAINS */