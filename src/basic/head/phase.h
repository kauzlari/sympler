/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2015, 
 * David Kauzlaric <david.kauzlaric@frias.uni-freiburg.de>,
 * and others authors stated in the AUTHORS file in the top-level 
 * source directory.
 *
 * SYMPLER is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SYMPLER is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SYMPLER.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Please cite the research papers on SYMPLER in your own publications. 
 * Check out the PUBLICATIONS file in the top-level source directory.
 *
 * You are very welcome to contribute extensions to the code. Please do 
 * so by making a pull request on https://github.com/kauzlari/sympler
 * 
 */



#ifndef __PHASE_H
#define __PHASE_H 

#include <set>
#include <string>

#include "general.h"
#include "boundary.h"
#include "particle_list.h"
#include "smart_pointer.h"
#include "node_many_children.h"
#include "random.h"
// the following gives bad include loops 
/* #include "triplet.h" */
/* #include "triplet_calc_angular_f.h" */


using namespace std;

//---- Macros ---- 

/* If the group list is empty we have to loop over all particles. If someone
   knows how to do it without an "if (groups.empty())" and therefore writing
   the same code twice let me know. */
	 
	 
	 
#define ALL_COLOURS ((size_t) 0xFFFFFFFF)
		


#define FOR_EACH_FREE_PARTICLE_C(phase, col, code)    \
if (col == ALL_COLOURS) {                             \
	for(size_t c = 0; c < phase->nColours(); ++c)       \
		SL_FOR_EACH(Particle, phase->particles(c), code); \
} else {                                              \
  size_t c = col;				      \
  SL_FOR_EACH(Particle, phase->particles(c), code); \
} while(0)

#define FOR_EACH_FREE_PARTICLE(phase, code)           \
{                                                     \
	for(size_t c = 0; c < phase->nColours(); ++c)       \
		SL_FOR_EACH(Particle, phase->particles(c), code); \
} while(0)

 
#define FOR_EACH_FROZEN_PARTICLE_ALL_C(phase, code)        \
{                                                     \
	for(size_t c = 0; c < phase->nColours(); ++c)       \
		SL_FOR_EACH(Particle, phase->frozenParticles(c), code); \
} while(0)

#define FOR_EACH_PARTICLE_C(phase, col, code)                           \
{                                                                     \
  if (col == ALL_COLOURS) {                                             \
    for(size_t c = 0; c < phase->nColours(); ++c) {                     \
      SL_FOR_EACH(Particle, phase->particles(c), code);                 \
      SL_FOR_EACH(Particle, phase->frozenParticles(c), code);           \
    }                                                                   \
  } else {                                                              \
    size_t c = col;                                                     \
    SL_FOR_EACH(Particle, phase->particles(c), code);                 \
    SL_FOR_EACH(Particle, phase->frozenParticles(c), code);           \
  }                                                                     \
} while(0)

#define FOR_EACH_PARTICLE(phase, code)                           \
{                                                                \
	for(size_t c = 0; c < phase->nColours(); ++c) {                \
    SL_FOR_EACH(Particle, phase->particles(c), code);            \
    SL_FOR_EACH(Particle, phase->frozenParticles(c), code);      \
  }                                                              \
} while(0)



/*
#define FOR_EACH_PARTICLE_SLOT_C(phase, c, code)                        \
{                                                                       \
  SL_FOR_EACH_IN_MEMORY(Particle, phase->particles(c), code);           \
  SL_FOR_EACH_IN_MEMORY(Particle, phase->frozenParticles(c), code);     \
} while(0)
*/
  

#define FOR_EACH_FROZEN_PARTICLE(phase, c, code)        \
SL_FOR_EACH(Particle, phase->frozenParticles(c), code)


#define FOR_EACH_FREE_PARTICLE_IN_GROUP(phase, col, groups, code)    \
if (col == ALL_COLOURS) {                                            \
  for (size_t c = 0; c < phase->nColours(); ++c) {                   \
    if (groups.empty()) {                                            \
      SL_FOR_EACH                                                    \
        (Particle,                                                   \
         phase->particles(c),                                        \
         code                                                        \
        );                                                           \
    } else {                                                         \
      SL_FOR_EACH                                                    \
        (Particle,                                                   \
         phase->particles(c),                                        \
         if (groups.find(__iSLFE->g) != groups.end()) {                    \
           code                                                      \
         }                                                           \
        );                                                           \
    }                                                                \
  }                                                                  \
} else {                                                             \
  if (groups.empty()) {                                              \
    SL_FOR_EACH                                                      \
      (Particle,                                                     \
       phase->particles(col),                                        \
       code                                                          \
      );                                                             \
  } else {                                                           \
    SL_FOR_EACH                                                      \
      (Particle,                                                     \
       phase->particles(col),                                        \
       if (groups.find(__iSLFE->g) != groups.end()) {                      \
         code                                                        \
       }                                                             \
      );                                                             \
  }                                                                  \
} while(0)



#define FOR_EACH_GROUP(groups, code)                                    \
{                                                                       \
  for (group_t::const_iterator group = groups.begin(); group != groups.end(); group++) { \
    code                                                                \
  }                                                                 \
} while(0)




//---- Phase ----

class Simulation;
class ManagerCell;
class PairCreator;
class IntegratorPosition;
class TripletCalculator;
class QuintetCalculator;
struct triplet_t;
struct quintet_t;

typedef list<triplet_t> tripletList;
typedef list<triplet_t>::iterator tripletListItr;  

typedef list<quintet_t> quintetList;
typedef list<quintet_t>::iterator quintetListItr;  

/*!
 * The \a Phase is the main object in charge of the whole simulation domain.
 * It keeps the main list of all particles, as well as instances of a \a Boundary
 * \a ManagerCell and \a PairCreator.
 */
class Phase: public NodeManyChildren
{
 protected:

  /*!
   * The center of mass velocity
   */
  point_t velCM;

  /*!
   * Does the center of mass velocity need to be recalculated
   */
  bool velCMold;

  /*!
   * The center of mass velocity per group
   */
  map<size_t, point_t> velCMPerGroup;
	
  /*!
   * Filename of the VTK to which cell subdivision information is stored
   */
  string m_cell_filename;

  /*!
   * List of all particles, sorted by color
   */
  vector<ParticleList> m_particles;
		
  /*!
   * List of all frozen particles, sorted by color
   *
   * frozen particles; should not need a reference list like the free
   * particles, since their number should not change during simulation
   * (at least for now 2006-02-22); this is CONVENTION 3
   */
  vector<ParticleList> m_frozen_particles;

  /*!
   * Stores lists of bonded triplets of \a Particle s
   */
  vector<tripletList*> m_tripletLists;

  /*!
   * String identifiers of the triplet lists in \a m_tripletLists 
   */
  vector<string> m_tripletListNames;

  /*!
   * The biggest stage occuring in \a m_bondedTripletCalculators . This is determined during runtime in Phase::sortStages()
   */
  vector<size_t> m_maxBondedStage_triplet;

  /*!
   * The biggest stage occuring in \a m_bondedTripletCalculators_0 . This is determined during runtime in Phase::sortStages_0()
   */
  vector<size_t> m_maxBondedStage_triplet_0;

  /*!
   * The biggest stage occuring in \a m_bondedQuintetCalculators . This is determined during runtime in Phase::sortStages()
   */
  vector<size_t> m_maxBondedStage_quintet;

  /*!
   * The biggest stage occuring in \a m_bondedQuintetCalculators_0 . This is determined during runtime in Phase::sortStages_0()
   */
  vector<size_t> m_maxBondedStage_quintet_0;


  /*!
   * Total number of particle
   */
  size_t nOfParticles;

  /*!
   * Total number of frozen particles
   */
  size_t nOfFrozenP;

  /*
  vector<size_t> nOfParticlesPerColour;
  vector<size_t> nOfFrozenPerColour;    
  */

  /*!
   * Total number of free particle sorted by group
   */
  map<size_t, size_t> nOfFreeParticlesPerGroup;
		
  /*!
   * Only create cells (for cell subdivision) that lie within the simulation
   * region
   */
  bool m_smartCells;

  /*!
   * Should the order of the pairs be randomised?
   */
  bool m_randomPairs;

  /*!
   * Should half cutoff width cells be created?
   */
  bool m_halfCutoff;
  /*!
   * This Calculators are for 
   * bonded triplets. An additional vector-layer is used for distinguishing 
   * between the different available bonded lists.
   */
  vector<vector<vector<TripletCalculator*> > > m_bondedTripletCalculators;

  /*!
   * As \a m_bondedTripletCalculators but called at another instant during one timestep by the \a Controller
  */
  vector<vector<vector<TripletCalculator*> > > m_bondedTripletCalculators_0;

  /*!
   * This is a helper for building \a m_bondedTripletCalculators . First, the \a ValCalculator s
   * are registered here and afterwards sorted by stages into \a m_bondedTripletCalculators .
   */
  vector<TripletCalculator*> m_bondedTripletCalculators_flat;

  /*!
   * This is a helper for building \a m_bondedTripletCalculators_0 . First, the \a ValCalculator s
   * are registered here and afterwards sorted by stages into \a m_bondedTripletCalculators .
   */
  vector<TripletCalculator*> m_bondedTripletCalculators_flat_0;


//------------------------------------------------------------------------------------------------------------------------
  /*!
   * Stores lists of bonded quintets of \a Particle s
   */
  vector<quintetList*> m_quintetLists;

  /*!
   * String identifiers of the quintet lists in \a m_quintetLists
   */
  vector<string> m_quintetListNames;

  /*!
   * This Calculators are for 
   * bonded quintet. An additional vector-layer is used for distinguishing 
   * between the different available bonded lists.
   */
  vector<vector<vector<QuintetCalculator*> > > m_bondedQuintetCalculators;

  /*!
   * As \a m_bondedQuintetCalculators but called at another instant during one timestep by the \a Controller
  */
  vector<vector<vector<QuintetCalculator*> > > m_bondedQuintetCalculators_0;

  /*!
   * This is a helper for building \a m_bondedQuintetCalculators . First, the \a ValCalculator s
   * are registered here and afterwards sorted by stages into \a m_bondedQuintetCalculators .
   */
  vector<QuintetCalculator*> m_bondedQuintetCalculators_flat;

  /*!
   * This is a helper for building \a m_bondedQuintetCalculators_0 . First, the \a ValCalculator s
   * are registered here and afterwards sorted by stages into \a m_bondedQuintetCalculators .
   */
  vector<QuintetCalculator*> m_bondedQuintetCalculators_flat_0;


//------------------------------------------------------------------------------------------------------------------------

  /*!
   * The cell subdivsion manager
   */
  ManagerCell* m_manager;

  /*!
   * Pointer to the \a PairCreator object
   */
  PairCreator *m_pairCreator;

  /*!
   * Pointer to the \a Boundary object
   */
  Boundary *m_boundary;
  
  /*!
   * My random number generator
   */
  RandomNumberGenerator m_rng;

  /*!
   * Initialize the property list
   */
  void init();

  /*!
   * Compute the sum of the squares of the velocities
   */
  void computeVVsum();
		
  /*!
   * Create a new particle of color \a colour
   * @param colour The color of the new particle
   */
  Particle* newParticle(size_t colour);

  /*!
   * Create the node corresponding to name \a name. In this case, this
   * Can be a \a Boundary, \a Phase, \a WeightingFunction, \a PairCreator etc.
   */
  virtual Node *instantiateChild(const string &name);
  
  /*!
   * Add a new colour (type of particle) to the simulation
   */
  void addColour(size_t colour);


public:
  /*!
   * Constructor
   * @param simulation Pointer to the parent simulation object
   */
  Phase(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~Phase();

  /*!
   * As the name says
   */
  const double& returnVelCMPerGroup(size_t particle, size_t group);

  /*!
   * As the name says
   */
  size_t returnNOfFreeParticlesPerGroup(size_t group);

  /*!
   * Invalidate tells phase that the coordinates have been updated.
   */
  void invalidate();

  /*!
   * Clear the tags for new calculation of cached quantities
   * currently, only additional doubles of particles are reset to zero
   * Fixme!!! what else could be done here?
   */
  void setForNewIntegration();

  /*!
   * Compute the center of mass velocity
   * Fixme!!! I don't like this function to be public
   */
  void computeVelCM();

  /*!
   * Add particle \a particle to the simulation and return a
   * pointer to its position in memory.
   * @param particle Add this particle
   */
  Particle *addParticle(const Particle &particle);

  /*!
   * Add a particle to the simulation and return a pointer to its position
   * in memory.
   * @param pos Position of the particle
   * @param vel Velocity of the particle
   * @param group Group of this particle
   * @param colour Color of this particle
   */ 
  Particle *addParticle
    (const point_t &pos, const point_t &vel, size_t group, size_t colour);

/*!
   * Add a frozen particle \a particle to the simulation and return a
   * pointer to its position in memory.
   * @param particle Add this particle
   */
  Particle *addFrozenParticle(const Particle &particle);

  /*!
   * Add a frozen particle to the simulation and return a pointer to its position
   * in memory.
   * @param pos Position of the particle
   * @param vel Velocity of the particle
   * @param group Group of this particle
   * @param colour Color of this particle
   */ 
  Particle *addFrozenParticle
    (const point_t &pos, const point_t &vel, size_t group, size_t colour);

  /*!
   * Change the group of particle \a p to \a newg
   * @param p Change group of this particle
   * @param newg New group of this particle
   */
  void groupChanged(Particle *p, size_t newg) {
    --nOfFreeParticlesPerGroup[p->g];
    ++nOfFreeParticlesPerGroup[newg];
  }

  /*!
   * Count all particles
   */
  void countParticles();

  /*!
   * Advances the position (including collision detection)
   * @param integrator Integrato to use for advancing positions
   */
  void invalidatePositions(IntegratorPosition *integrator);

  /*!
   * Is the center of mass velocity still valid. Obsolete?
   */
  bool velCMIsOld();

  /*!
   * Is there more than one group in the system?
   */
  bool groupsDefined() const;

  /*!
   * Return the square of the velocities
   */
  double returnVVsum();

  /*!
   * Return the square of the velocities for group \a group
   * @param group Group to return the sum for
   */
  double returnVVsum(size_t group);

  /*!
   * Return the square of the velocities for groups \a groups
   * @param groups Set of groups of the particles
   */
  double returnVVsum(const group_t &groups);

  /*!
   * Return the number of particles
   */
  size_t returnNofPart() const;

  /*!
   * Return whether assigneParticlesToCells has been called
   */
  bool m_particles_assigned;

  /*!
   * Return the number of particles for color \a colour
   * @param colour Color of the particles
   */
  size_t returnNofPartC(size_t colour);

  /*!
   * Return the number of particles for group \a group
   * @param group Group of the particles
   */
  size_t returnNofPart(size_t group);

  /*!
   * Return the number of particles for groups \a groups
   * @param groups Set of groups of the particles
   */
  size_t returnNofPart(const group_t &groups);

  /*!
   * Return the number of frozen particles
   */
  size_t returnNofFrozenP() const;

  /*!
   * Return the number of frozen particle for color \a colour
   * @param colour Color of the particles
   */
  size_t returnNofFrozenPC(size_t colour) const;

  /*!
   * Return the average density in the simulation. Fixme!!! Obsolete?
   */
  double returnDensity() const;

  /*!
   * Return the volume of the bounding box of the simulation. Fixme!!! Obsolete?
   */
  double cuboidVolume();

  /*!
   * Check if there is any boundary defined
   */
  virtual void setup();

  /*!
   * Remove particle \a p
   * @param p Particle to remove
   */
  void removeParticle(Particle *p);

  /*!
   * Remove particle of color \a colour and slot \a slot.
   * @param colour Color of the particle to be removed
   * @param slot List slot of the particle to be removed
   */
  void removeParticle(size_t colour, size_t slot);

  /* "New-style" interface to particles. */
  /*!
   * Return the free particle with color \a c and slot \a i
   * @param c Color of the particle
   * @param i List slot of the particle
   */
  Particle& freeP(size_t c, int i) {
    return m_particles[c][i];
  }

  /*!
   * Return the free particle with color \a c and slot \a i
   * @param c Color of the particle
   * @param i List slot of the particle
   */
  const Particle& freeP(size_t c, int i) const {
    return m_particles[c][i];
  }

  /*!
   * Return the frozen particle with color \a c and slot \a i
   * @param c Color of the particle
   * @param i List slot of the particle
   */
  Particle& frozenP(size_t c, int i) {
    return m_frozen_particles[c][i];
  }

  /*!
   * Return the frozen particle with color \a c and slot \a i
   * @param c Color of the particle
   * @param i List slot of the particle
   */
  const Particle& frozenP(size_t c, int i) const {
    return m_frozen_particles[c][i];
  }
		
  /*!
   * Return a list of all free particles
   */
  ParticleList &particles(size_t c) {
    return m_particles[c];
  }

  /*!
   * Return a list of all frozen particles
   */
  ParticleList &frozenParticles(size_t c) {
    return m_frozen_particles[c];
  }

  /*!
   * Return the number of colors present in the system
   */
  size_t nColours() const;

  /*!
   * Return the center of mass velocity for the system
   */
  const point_t &centerOfMassVelocity() const {
    return velCM;
  }

  /*!
   * Return the center of mass velocity for particles in group \a g
   * @param g Group
   */
  const point_t &centerOfMassVelocity(int g) {
    return velCMPerGroup[g];
  }

  /*!
   * Return a pointer to the \a Boundary object
   */
  Boundary *boundary();

  /*!
   * Return a pointer to the cell subdivision manager
   */
  inline ManagerCell *manager() {
    return m_manager;
  }

  /*!
   * Return a pointer to the cell subdivision manager
   */
//   PairCreator *pairCreator() {
//     return m_pairCreator;
//   }
  PairCreator *pairCreator();
  
  /*!
   * Returns the \a Phase 's \a RandomNumberGenerator
   */
  RandomNumberGenerator& rng()
  {
    return m_rng; 
  }
  
  /*!
   * Will the cells be only created in the simulation domain?
   */
  bool smartCells() {
    return m_smartCells;
  }

  /*!
   * Should the order of the pairs be randomised?
   */
  bool randomPairs() {
    return m_randomPairs;
  }

  virtual bool& particlesAssigned() {
    return m_particles_assigned;
  }

  /*!
   * Initially, sort particles to the cell subdivision
   */
  virtual void assignParticlesToCells();

  /*!
   * Read from XML file and initialize the cell subdivision manager
   */
  virtual void read(const xmlNode *xmln);

  /*!
   * Write a restart file that can be used by \a ParticleCreatorFile
   */
  virtual void writeRestartFile(string name);


  /*!
   * Return the capacity of the particle array
   */
  virtual size_t returnArrayCapacity(size_t c) {
    return m_particles[c].capacity();
  }


  vector<TripletCalculator*>* bondedTripletCalculatorsFlat() {
    return &m_bondedTripletCalculators_flat;
  }

  vector<TripletCalculator*>* bondedTripletCalculatorsFlat_0() {
    return &m_bondedTripletCalculators_flat_0;
  }

  /*!
   * Return a pointer to \a m_tripletLists
   */
  vector<tripletList*>* tripletLists() {
    return &m_tripletLists;
  }

  /*!
   * Return the connected list from \a m_connectedLists at slot \a slot.
   */
  tripletList* returnTripletList(size_t slot) {
    return m_tripletLists[slot];
  }
  
  /*!
   * Add the three \a Particle s to the \a tripletList of \a m_tripletLists indexed by \a listIndex
   */
  void addTriplet(Particle *p1, Particle *p2, Particle *p3, size_t listIndex);

  /*!
   * Add the three \a Particle s to the \a tripletList of \a m_tripletLists indexed by \a listIndex and write into a file
   */
  void addTripletAndWrite(Particle *p1, Particle *p2, Particle *p3, size_t listIndex);

  /*!
   * Returns the string identifier of the \a tripletList from \a m_tripletLists with index \a listIndex
   */
  string tripletListName(size_t listIndex);

  /*!
   * Returns the index of the \a tripletList from \a m_tripletLists with string identifier \a name
   */
  size_t tripletListIndex(string name);

  /*!
   * Returns the index of the \a tripletList from \a m_tripletLists with string identifier \a name . If it does not exist yet it is created and writte to a file.
   */
  size_t createTripletListIndexAndWrite(string name);

  /*!
   * Returns the index of the \a tripletList from \a m_tripletLists with string identifier \a name . If it does not exist yet it is created.
   */
  size_t createTripletListIndex(string name);

  /*!
   * Returns the index of the \a tripletList from \a m_tripletLists with string identifier \a name
   */
  int searchTripletListIndex(string name);

  /*!
   * Returns the \a tripletList from \a m_tripletLists with string identifier \a name
   */
  tripletList* returnTripletList(string name);

    
  /*!
   * Register the bonded triplet calculators
   */
  void registerBondedCalc(TripletCalculator* vC);

  /*!
   * Register the bonded triplet calculators for instant "_0"
   */
  void registerBondedCalc_0(TripletCalculator* vC);

 //----------------------------------------------------------------------------------------------------------------------------
  /*!
   * Quintet Calculator Definition of lists an itterators
   */

  vector<QuintetCalculator*>* bondedQuintetCalculatorsFlat() {
    return &m_bondedQuintetCalculators_flat;
  }

  vector<QuintetCalculator*>* bondedQuintetCalculatorsFlat_0() {
    return &m_bondedQuintetCalculators_flat_0;
  }

  /*!
   * Return a pointer to \a m_quintetLists
   */
  vector<quintetList*>* quintetLists() {
    return &m_quintetLists;
  }

  /*!
   * Return the connected list from \a m_connectedLists at slot \a slot.
   */
  quintetList* returnQuintetList(size_t slot) {
    return m_quintetLists[slot];
  }
  
  /*!
   * Add the three \a Particle s to the \a quintetList of \a m_quintetLists indexed by \a listIndex
   */
  void addQuintet(Particle *p1, Particle *p2, Particle *p3, Particle *p4, Particle *pc, size_t listIndex);

  /*!
   * Add the three \a Particle s to the \a quintetList of \a m_quintetLists indexed by \a listIndex and write into a file
   */
  void addQuintetAndWrite(Particle *p1, Particle *p2, Particle *p3, Particle *p4, Particle *pc, size_t listIndex);

  /*!
   * Returns the string identifier of the \a quintetList from \a m_quintetLists with index \a listIndex
   */
  string quintetListName(size_t listIndex);

  /*!
   * Returns the index of the \a quintetList from \a m_quintetLists with string identifier \a name
   */
  size_t quintetListIndex(string name);

  /*!
   * Returns the index of the \a quintetList from \a m_quintetLists with string identifier \a name . If it does not exist yet it is created and writte to a file.
   */
  size_t createQuintetListIndexAndWrite(string name);

  /*!
   * Returns the index of the \a quintetList from \a m_quintetLists with string identifier \a name . If it does not exist yet it is created.
   */
  size_t createQuintetListIndex(string name);

  /*!
   * Returns the index of the \a quintetList from \a m_quintetLists with string identifier \a name
   */
  int searchQuintetListIndex(string name);

  /*!
   * Returns the \a quintetList from \a m_quintetLists with string identifier \a name
   */
  quintetList* returnQuintetList(string name);

    
  /*!
   * Register the bonded quintet calculators
   */
  void registerBondedQuinCalc(QuintetCalculator* vC);

  /*!
   * Register the bonded quintet calculators for instant "_0"
   */
  void registerBondedQuinCalc_0(QuintetCalculator* vC);


  /*!
   * Return the list of calculators with index \a listIndex 
   * from \a m_bondedQuintetCalculators for stage \a stage. 
   */
  vector<QuintetCalculator*>& bondedQuintetCalculators(size_t stage, size_t listIndex) {
    return (m_bondedQuintetCalculators[listIndex][stage]);
  }

  /*!
   * Return the list of calculators with index \a listIndex 
   * from \a m_bondedQuintetCalculators_0 for stage \a stage. 
   */
  vector<QuintetCalculator*>& bondedQuintetCalculators_0(size_t stage, size_t listIndex) {
    return (m_bondedQuintetCalculators_0[listIndex][stage]);
  }



//----------------------------------------------------------------------------------------------------------------------------


  /*!
   * Sort the stages of the Calculators
   */
  void sortStages();

  /*!
   * Sort the stages of the Calculators called at instant "_0"
   */
  void sortStages_0();
      
  /*!
   * return the \a m_maxBondedStage_triplet[icl] 
   */
  size_t maxBondedStage_triplet(size_t icl) const
  { 
    return m_maxBondedStage_triplet[icl];
  }
      
  /*!
   * return the \a m_maxBondedStage_triplet_0[icl] 
   */
  size_t maxBondedStage_triplet_0(size_t icl) const
  { 
    return m_maxBondedStage_triplet_0[icl];
  }

  /*!
   * return the \a m_maxBondedStage_quintet[icl] 
   */
  size_t maxBondedStage_quintet(size_t icl) const
  { 
    return m_maxBondedStage_quintet[icl];
  }
      
  /*!
   * return the \a m_maxBondedStage_quintet_0[icl] 
   */
  size_t maxBondedStage_quintet_0(size_t icl) const
  { 
    return m_maxBondedStage_quintet_0[icl];
  }
      
  /*!
   * Return the list of calculators with index \a listIndex 
   * from \a m_bondedTripletCalculators for stage \a stage. 
   */
  vector<TripletCalculator*>& bondedTripletCalculators(size_t stage, size_t listIndex) {
    return (m_bondedTripletCalculators[listIndex][stage]);
  }

  /*!
   * Return the list of calculators with index \a listIndex 
   * from \a m_bondedTripletCalculators_0 for stage \a stage. 
   */
  vector<TripletCalculator*>& bondedTripletCalculators_0(size_t stage, size_t listIndex) {
    return (m_bondedTripletCalculators_0[listIndex][stage]);
  }

  /*!
   * Find the stages of the calculators stored in the \a Phase
   */
  bool findStages();

  /*!
   * Find the stages of the calculators stored in the \a Phase
   */
  bool findStages_0();

  /*!
   * This variable states whether the declaration of the lists containing connections is already completely written in the output file containing particle connections.
   */
  static bool connectionDeclarationFinished;

  friend class ManagerCell;

};



//---- inline functions ----------------------------------------------------

inline bool Phase::velCMIsOld()
{
  return velCMold;
}

inline size_t Phase::returnNOfFreeParticlesPerGroup(size_t group)
{
  return nOfFreeParticlesPerGroup[group];
}


// number of FREE particles only
inline size_t Phase::returnNofPart() const
{
	return nOfParticles;
}

inline size_t Phase::returnNofFrozenP() const
{
	return nOfFrozenP;
}

inline size_t Phase::returnNofFrozenPC(size_t colour) const
{
	return m_frozen_particles[colour].size();
}

// number of FREE particles only
inline size_t Phase::returnNofPart(size_t group)
{
		return nOfFreeParticlesPerGroup[group];
}

inline size_t Phase::returnNofPartC(size_t colour)
{
// 		MSG_DEBUG("Phase::returnNofPartC", "in");
// 		return m_particles[colour]->size();
//	return nOfParticlesPerColour[colour];
  return m_particles[colour].size();
}

// number of FREE particles only
inline size_t Phase::returnNofPart(const group_t &groups)
{
    size_t nofp = 0;

    if (groups.empty()) {
        return nOfParticles;
    } else {
        FOR_EACH_GROUP(groups,
            nofp += nOfFreeParticlesPerGroup[*group];
        );
    }

    return nofp;
}

// we try to get rid of the following commented out functions

/* //CONVENTION 4: following function only needed for free particles */
/* inline double Phase::returnVVsum() */
/* { */
/* 		if (vvSumOld) */
/* 				computeVVsum(); */
/* 		return vvSum; */
/* } */

/* //CONVEBTION 4: following function only needed for free particles */
/* inline double Phase::returnVVsum(size_t group) */
/* { */
/* 		if (vvSumOld) */
/* 				computeVVsum(); */
/* 		return vvSumPerGroup[group]; */
/* } */

/* //CONVENTION 4: following function only needed for free particles */
/* inline double Phase::returnVVsum(const group_t &groups) */
/* { */
/*     double vvsum = 0; */

/*     if (vvSumOld) */
/*         computeVVsum(); */

/*     if (groups.empty()) { */
/*         return vvSum; */
/*     } else { */
/*         FOR_EACH_GROUP(groups, */
/*             vvsum += vvSumPerGroup[*group]; */
/*         ); */
/*     } */

/*     return vvsum; */
/* } */


inline double Phase::cuboidVolume()
{
		return boundary() -> cuboidVolume();
}


#endif
