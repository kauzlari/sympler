/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2013, 
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



#ifndef __MANAGER_CELL_H
#define __MANAGER_CELL_H

using namespace std;


#include <list>
#include "gen_f.h"
#include "pairdist.h"
#include "pair_list.h"
#include "integrator.h"
/*#include "colour_pair.h"*/
#include "wall_container.h"
// #include "integrator_position.h"
// #include "integrator_velocity_verlet.h"
// #include "integrator_velocity_verlet_disp.h"
// #include "integrator_static.h"


#define OUTLET 0
#define NEIGHBOR 1

#define P_CREATE 0
#define P_DELETE 1

#//---- Constants ----

#define NUM_DIRECT_NEIGHBORS 26
#define NUM_HALF_DIRECT_NEIGHBORS 13
#define OFFSET2DIRECTNEIGHBOR(off, n) \
    n = (off[0]+1)*9 + (off[1]+1)*3 + (off[2]+1); \
    if (n > NUM_DIRECT_NEIGHBORS/2) \
    n--; \
    while(0)


#define INV_DIRECT_NEIGHBOR(n)   NUM_DIRECT_NEIGHBORS-n-1
#define INV_NEIGHBOR(n)   m_manager->num_neighbors()-n-1
//TODO: make pairs to point to 2x
const int_point_t c_direct_offsets[NUM_DIRECT_NEIGHBORS] = {
    {-1,-1,-1}, {-1,-1, 0}, {-1,-1, 1}, {-1, 0,-1}, {-1, 0, 0}, {-1, 0, 1},
    {-1, 1,-1}, {-1, 1, 0}, {-1, 1, 1},
    { 0,-1,-1}, { 0,-1, 0}, { 0,-1, 1}, { 0, 0,-1}, /*{0, 0, 0},*/ { 0, 0, 1},
    { 0, 1,-1}, { 0, 1, 0}, { 0, 1, 1},
    { 1,-1,-1}, { 1,-1, 0}, { 1,-1, 1}, { 1, 0,-1}, { 1, 0, 0}, { 1, 0, 1},
    { 1, 1,-1}, { 1, 1, 0}, { 1, 1, 1}
};

const int c_1x_direct_neighbors[NUM_DIRECT_NEIGHBORS] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, };

#define NUM_2X_NEIGHBORS 124
#define NUM_HALF_2X_INDIRECT_NEIGHBORS 49 
//TODO::make pairs to point to direct neighbor
const pair<int, int_point_t> c_2x_offsets[NUM_2X_NEIGHBORS] = {
    {-1, {-2,-2,-2}}, {-1, {-2,-2,-1}}, {-1, {-2,-2, 0}}, {-1, {-2,-2, 1}}, {-1, {-2,-2, 2}},
    {-1, {-2,-1,-2}}, {-1, {-2,-1,-1}}, {-1, {-2,-1, 0}}, {-1, {-2,-1, 1}}, {-1, {-2,-1, 2}},
    {-1, {-2, 0,-2}}, {-1, {-2, 0,-1}}, {-1, {-2, 0, 0}}, {-1, {-2, 0, 1}}, {-1, {-2, 0, 2}},
    {-1, {-2, 1,-2}}, {-1, {-2, 1,-1}}, {-1, {-2, 1, 0}}, {-1, {-2, 1, 1}}, {-1, {-2, 1, 2}},
    {-1, {-2, 2,-2}}, {-1, {-2, 2,-1}}, {-1, {-2, 2, 0}}, {-1, {-2, 2, 1}}, {-1, {-2, 2, 2}},

    {-1, {-1,-2,-2}}, {-1, {-1,-2,-1}}, {-1, {-1,-2, 0}}, {-1, {-1,-2, 1}}, {-1, {-1,-2, 2}},
    {-1, {-1,-1,-2}}, {0, {-1,-1,-1}}, {1, {-1,-1, 0}}, {2, {-1,-1, 1}}, {-1, {-1,-1, 2}},
    {-1, {-1, 0,-2}}, {3, {-1, 0,-1}}, {4, {-1, 0, 0}}, {5, {-1, 0, 1}}, {-1, {-1, 0, 2}},
    {-1, {-1, 1,-2}}, {6, {-1, 1,-1}}, {7, {-1, 1, 0}}, {8, {-1, 1, 1}}, {-1, {-1, 1, 2}},
    {-1, {-1, 2,-2}}, {-1, {-1, 2,-1}}, {-1, {-1, 2, 0}}, {-1, {-1, 2, 1}}, {-1, {-1, 2, 2}},

    {-1, { 0,-2,-2}}, {-1, { 0,-2,-1}}, {-1, { 0,-2, 0}}, {-1, { 0,-2, 1}}, {-1, { 0,-2, 2}},
    {-1, { 0,-1,-2}}, {9, { 0,-1,-1}}, {10, { 0,-1, 0}}, {11, { 0,-1, 1}}, {-1, { 0,-1, 2}},
    {-1, { 0, 0,-2}}, {12, { 0, 0,-1}},/*{0, 0, 0},*/{13, {0, 0, 1}}, {-1, { 0, 0, 2}},
    {-1, { 0, 1,-2}}, {14, { 0, 1,-1}}, {15, { 0, 1, 0}}, {16, { 0, 1, 1}}, {-1, { 0, 1, 2}},
    {-1, { 0, 2,-2}}, {-1, { 0, 2,-1}}, {-1, { 0, 2, 0}}, {-1, { 0, 2, 1}}, {-1, { 0, 2, 2}},

    {-1, { 1,-2,-2}}, {-1, { 1,-2,-1}}, {-1, { 1,-2, 0}}, {-1, { 1,-2, 1}}, {-1, { 1,-2, 2}},
    {-1, { 1,-1,-2}}, {17, { 1,-1,-1}}, {18, { 1,-1, 0}}, {19, { 1,-1, 1}}, {-1, { 1,-1, 2}},
    {-1, { 1, 0,-2}}, {20, { 1, 0,-1}}, {21, { 1, 0, 0}}, {22, { 1, 0, 1}}, {-1, { 1, 0, 2}},
    {-1, { 1, 1,-2}}, {23, { 1, 1,-1}}, {24, { 1, 1, 0}}, {25, { 1, 1, 1}}, {-1, { 1, 1, 2}},
    {-1, { 1, 2,-2}}, {-1, { 1, 2,-1}}, {-1, { 1, 2, 0}}, {-1, { 1, 2, 1}}, {-1, { 1, 2, 2}},

    {-1, { 2,-2,-2}}, {-1, { 2,-2,-1}}, {-1, { 2,-2, 0}}, {-1, { 2,-2, 1}}, {-1, { 2,-2, 2}},
    {-1, { 2,-1,-2}}, {-1, { 2,-1,-1}}, {-1, { 2,-1, 0}}, {-1, { 2,-1, 1}}, {-1, { 2,-1, 2}},
    {-1, { 2, 0,-2}}, {-1, { 2, 0,-1}}, {-1, { 2, 0, 0}}, {-1, { 2, 0, 1}}, {-1, { 2, 0, 2}},
    {-1, { 2, 1,-2}}, {-1, { 2, 1,-1}}, {-1, { 2, 1, 0}}, {-1, { 2, 1, 1}}, {-1, { 2, 1, 2}},
    {-1, { 2, 2,-2}}, {-1, { 2, 2,-1}}, {-1, { 2, 2, 0}}, {-1, { 2, 2, 1}}, {-1, { 2, 2, 2}},
};

const vector<pair<int, int> > interMap = {
    {0, 0}, {0, 1}, {1, 1}, {2, 1}, {2, 2},
    {0, 3}, {0, 4}, {1, 4}, {2, 4}, {2, 5},
    {3, 3}, {3, 4}, {4, 4}, {5, 4}, {5, 5},
    {6, 3}, {6, 4}, {7, 4}, {8, 4}, {8, 5},
    {6, 6}, {6, 7}, {7, 7}, {8, 7}, {8, 8},

    {0,  9}, {0, 10}, {1, 10}, {2, 10}, {2, 11},
    {0, 12}, {-0, -13}, {-1, -13}, {-2, -13}, {2, 13},
    {3, 12}, {-3, -13}, {-4, -13}, {-5, -13}, {5, 13},
    {6, 12}, {-6, -13}, {-7, -13}, {-8, -13}, {8, 13},
    {6, 14}, {6, 15}, {7, 15}, {8, 15}, {8, 16},

    {9, 9}, {9, 10}, {10, 10}, {11, 10}, {11, 11},
    {9, 12},{-9, -13}, {-10, -13}, {-11, -13}, {11, 13}, 
    {12, 12}, {-12, -13}, /*{13,13}*/{-13, -13}, {13, 13},
    {14, 12}, {-14, -13}, {-15, -13}, {-16, -13}, {16, 13},
    {14, 14}, {14, 15}, {15, 15}, {16, 15}, {16, 16},

    {17, 9}, {17, 10}, {18, 10}, {19, 10}, {19, 11},
    {17, 12}, {-17, -13}, {-18, -13}, {-19, -13}, {19, 13},
    {20, 12}, {-20, -13}, {-21, -13}, {-22, -13}, {22, 13},
    {23, 12}, {-23, -13}, {-24, -13}, {-25, -13}, {25, 13},
    {23, 14}, {23, 15}, {24, 15}, {25, 15}, {25, 16},

    {17, 17}, {17, 18}, {18, 18}, {19, 18}, {19, 19},
    {17, 20}, {17, 21}, {18, 21}, {19, 21}, {19, 22},
    {20, 20}, {20, 21}, {21, 21}, {22, 21}, {22, 22},
    {23, 20}, {23, 21}, {24, 21}, {25, 21}, {25, 22},
    {23, 23}, {23, 24}, {24, 24}, {25, 24}, {25, 25},
};        

const int c_indirect_neighbors[NUM_HALF_2X_INDIRECT_NEIGHBORS] = {
     0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
     25,26,27,28,29,30,34,35,39,40,44,45,46,47,48,49,
     50,51,52,53,54,55,59,60,};/*63,64,68,69,70,71,72,73,
     74,75,76,77,78,79,83,84,88,89,93,94,95,96,97,98,
     99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,
       115,116,117,118,119,120,121,122,123,};*/

const int c_2x_direct_neighbors[NUM_DIRECT_NEIGHBORS] = {
     31, 32, 33, 36, 37, 38, 41, 42, 43, 
     56, 57, 58, 61, 62, 65, 66, 67,
     80, 81, 82, 85, 86, 87, 90, 91, 92, };
//---- Macros ----

#define FOR_EACH_COLOUR_PAIR(manager, code)                             \
{                                                                       \
  vector<ColourPair*>::iterator __end = manager->colourPairs().end();   \
  for(vector<ColourPair*>::iterator __cp = manager->colourPairs().begin(); \
      __cp != __end; ++__cp)                                            \
	{                                                                     \
    ColourPair *cp = *__cp;                                             \
    code                                                                \
  }                                                                     \
} while(0)



#define TOCELLINDEX(p, n)  ((p.z * n.y) + p.y) * n.x + p.x
#define TOCELLINDEXD(x, y, z, n)  ((z * n.y) + y) * n.x + x

class Cell;
/*!
 * This structure holds information about a region in space. A region
 * is a cuboid which is in return subdivided by cells.
 */
struct region_t {
  bool_point_t m_oneCellPeriodicDims; 
  
  /*!
   * First corner of this region
   */
  point_t corner1;

  /*!
   * Second corner of this region
   */
  point_t corner2;

  /*!
   * Number of cells in x, y and z direction
   */
  int_point_t n_cells;

  /*!
   * The inverse of a cell width in x, y and z direction
   */
  point_t inv_width;

  /*!
   * A list of all the cells in this region
   */
  vector<Cell*> cells;

  /*!
   * A map giving the index in the \a cells vector from
   * the direct cell index. This is needed because not all
   * cells within a region need to exist.
   */
  vector<int> cells_by_pos;

  /*!
   * Check whether the point \a pos lies inside this region
   * @param pos Point to check whether inside or outside
   */
  bool isInside(const point_t &pos) {
    bool match = true;

    for (int i = 0; (i < SPACE_DIMS) && match; i++) {
      if ((pos[i] < corner1[i]) || (pos[i] >= corner2[i]))
        match = false;
    }
/*    if(!match) MSG_DEBUG("region_t::isInside", "match = " << match << endl << "corner1 = " << corner1 << endl
   << "corner2 = "<< corner2 << endl << "pos = " << pos);       */
    return match;
  }

  /*!
   * Return the located at cell position \a p
   * @param p Position of the cell
   */
  Cell *cellByPos(const int_point_t &p) {
    int index = TOCELLINDEX(p, n_cells);
    int i = cells_by_pos[index];

    if (i == -1)
      return NULL;
    else
      return cells[i];
  }

  /*!
   * Return the cell located at space position \a pos
   * @param pos Space position of the cell
   */
  Cell *cellAtPos(const point_t &pos) {
    int_point_t p;

    assert(isInside(pos));

    p.x = (int) (inv_width.x * (pos.x - corner1.x));
    p.y = (int) (inv_width.y * (pos.y - corner1.y));
    p.z = (int) (inv_width.z * (pos.z - corner1.z));

    return cellByPos(p);
  }
};



//---- ManagerCell ----

class ParticleCreator;
class ColourPair;
class AbstractCellLink;
/*!
 * The \a ManagerCell is in charge of cell subdivision for the whole simulation
 * region.
 */
class ManagerCell
{
protected:
  /*!
   * common constructor code
   */
  void init();

  int m_divby;
  const int *m_direct_neighbors = c_1x_direct_neighbors;

  /*!
   * List of all cells in the simulation region
   */
  vector<Cell*> m_cells;

  /*!
   * First active cell, active cells are given by a linked list.
   */
  Cell* m_first_cell;

  /*!
   * Number of active cells, i.e., cells with particle, in the
   * simulation region.
   */
  size_t m_n_active_cells;

#ifdef ENABLE_PTHREADS
  /*!
   * Mutex for access to the cell list
   */
  pthread_mutex_t m_cells__mutex;
#endif

  /*!
   * List of all cell links in the system
   */
  vector<AbstractCellLink*> m_links;

  /*!
   * Number of active cell links, i.e., cell links where both of the connected
   * cells contain at least one particle
   */
#ifdef _OPENMP
  vector<size_t> m_n_active_links;
#else   
  size_t m_n_active_links;
#endif

  /*!
   * First active cell link, active cell links are given by a linked list.
   */
#ifdef _OPENMP   
  vector<AbstractCellLink*> m_first_link;
#else
  AbstractCellLink* m_first_link;
#endif  
#ifdef ENABLE_PTHREADS
  /*!
   * Mutex for access to the cell link list
   */
  pthread_mutex_t m_links__mutex;
#endif

  /*!
   * List of all regions
   */
  list<region_t*> m_regions;

  /*!
   * Pointer to the parent phase
   */
  Phase *m_phase;

  /*!
   * Name <-> Colour correspondence table.
   */
  vector<string> m_species;

  /*!
   * The objects managing the different pairs of colours.
   */
  vector<ColourPair*> m_colourPairs;
  
  /*!
   * (SameColourPair) 
   * Check if the colour pairs in the manager are in the same order as in Pairdist
   * Needed for "Forces"->computeForces()
   */
  bool m_scp; 

  /*!
   * Activate cell \a c, this means that the first particle entered this cell
   * @param c Cell to activate
   */
  void activateCell(Cell *c);

  /*!
   * Deactivate cell \a c, this means that the last particle left this cell
   * @param c Cell to deactivate
   */
  void deactivateCell(Cell *c);

  /*!
   * Activate cell link \a c, this means both cells have been activated.
   * @param c Cell link to activate
   */
  void activateCellLink(AbstractCellLink *c);

  /*!
   * Deactivate cell link \a c, this means one of both cells has been deactivated.
   * @param c Cell link to deactivate
   */
  void deactivateCellLink(AbstractCellLink *c);

    /*!
   * Find the species with name \a species and return its color. Adds the species
   * as new if it is not found.
   * @param species Name of the species
     */
  size_t getColourAndAdd(string species) {
    int r = findSpecies(species);

    if (r == -1)
      r = addColour(species);

    return r;
  }

  /*!
   * Add the species \a species and return its color
   * @param species Name of the new species
   */
  virtual size_t addColour(string species);

  void addBoundaryCell(region_t *r, cuboid_t cuboid, int_point_t cell_pos, int cellIndex, int group);
  /*!
   * This function creates the colours for firstS and secondS if they do not exist.
   * @param firstS Name of the first species
   * @param secondS Name of the second species
   */
  virtual ColourPair *cp(string firstS, string secondS);

  /*!
   * This function creates the colours firstS and secondS if they do not exist.
   * @param species Pair of the names of the first and the second species
   */
  virtual ColourPair *cp(pair<string, string> species) {
    return cp(species.first, species.second);
  }

  /*!
   * Connect the cells that belong to two regions but only one plane pair in both of them
   * @param a First region
   * @param b Second region
   * @param periodic Are the two region periodic in x, y or z direction?
   * @param t Create outlet or link? Outlet means, particle can escape, link means particles interact as well
   * @param dir Direction in which to connect the two regions
   */
  void connectPlanePair(region_t *a, region_t *b, bool_point_t periodic, AbstractCellLink *(Cell::*add)(Cell *, int, bool), void (ManagerCell::*dist)(Cell *, Cell*, int, int, int), int klbound, int_point_t p, int_point_t q, int_point_t off, int_point_t dir, int orderInter = 0);

public:
  /*!
   * Constructor
   * @param p Pointer to the parent phase object
   */
  ManagerCell(Phase* p);
  /*!
   * Constructor
   * @param p Pointer to the parent phase object
   * @param divby Cell width will be >=cutoff/divby, only 1 and 2 is supported, due to performance restrictions.
   */
  ManagerCell(Phase* p, int divby);

  /*!
   * Destructor
   */
  virtual ~ManagerCell();

  /*!
   * Returns total number of neighbours (direct + indirect)
   */
  virtual inline int num_neighbors() { return NUM_DIRECT_NEIGHBORS; };

  /*!
   * Returns the 3-coord offset vector of a specified neighbour \e a.
   * @param a Index of neighbour
   */
  virtual int_point_t c_offsets(int a) { return c_direct_offsets[a]; };
  virtual int c_2x_1x(int a) { return a; };

  /*!
   * Returns offset in specified direction \e b of 3-coord offset vector of specified neighbour \e a.
   * @param a Index of neighbour
   * @param b Index of direction
   */
  virtual int c_offsets (int a, int b) { return c_direct_offsets[a][b]; }
  
  /*!
   * Finds and returns the index of neighbour specified by offset vector \e off
   * @param off 3-coord int offset vector
   */
  virtual inline int offset2neighbor(int_point_t off){
    int n = (off[0]+1)*9 + (off[1]+1)*3 + (off[2]+1);
    if (n > NUM_DIRECT_NEIGHBORS/2)
      n--;
    return n;
  };
  
  /*!
   * Determines the cell pitch in real coordinates
   * @param first Pointer to the first Cell object
   * @param second Pointer to the second Cell object
   * @param alignment Index of neighbor, which \e second is for \e first. 
   * @param[out] 3 real element vector. For each coordinate, width of cell which is on the left.
   */
  virtual void cellDist(Cell *first, Cell *second, int alignment, int orderTarget, int orderInter = 0);
  virtual void cellDist2x(Cell *first, Cell *second, int alignment, int orderTarget, int orderInter) {};

  /*!
   * Advance the position of the particles, taking care of walls in the system
   * and set \a m_distances_valid to false 
   * @param integrator Position integrator to use
   */
  virtual void invalidatePositions(IntegratorPosition *integrator);
  
  /*!
   * Thread-number counter
   */
  static int thread_counter;

  /*!
   * Create the pair list
   */
  //virtual void createDistances();
  void findRegionCellDirectNeighbors(region_t *r, bool_point_t periodic, bool oneCellPeriodicDims);
  virtual void findRegionCellIndirectNeighbors(region_t *r, bool_point_t periodic, bool oneCellPeriodicDims) {};
  /*!
   * Create a new cell subdivided region
   * @param cutoff Cut-off radius (i.e., cell size)
   * @param corner1 First corner of the region
   * @param corner2 Second corner of the region
   * @param periodic Periodicity in x, y and z direction of the region
   * @param group The group all cells in this region should belong to
   * @param pc Particle creator for inlet/outlet definition. If = NULL => no inlet/outlet
   * @param axis Axis in which the inlet/outlet can be found
   * @param dir Direction in which the inlet/outlet can be found
   * @param action Create an inlet or an outlet?
   */
  region_t *cellSubdivide
    (double cutoff, point_t corner1, point_t corner2, bool_point_t periodic,
     int group = 0, ParticleCreator *pc = NULL, int axis = 0, int dir = 0, int action = P_CREATE);

  /*!
   * Connect the cells that belong to two regions
   * @param dir Direction in which to connect the two regions
   * @param a First region
   * @param b Second region
   * @param periodic Are the two region periodic in x, y or z direction?
   * @param t Create outlet or link? Outlet means, particle can escape, link means particles interact as well
   */
  void connect(int dir, region_t *a, region_t *b, bool_point_t periodic, int t);

  /*!
   * This methods tells the manager that the boundary is done
   * preparing the cell subdivision. The boundary will sort the
   * newly created cells in space.
   */
  void cellSubdivisionFinished();

  /*!
   * This method tells the manager the exact form of the boundayr
   * @param container Container object holding the triangulated surface
   */
  void assignContainer(WallContainer *container);

  /*!
   * Find a cell within the simulation domain
   * @param pos Find cell that contains this position
   * @param region If not NULL, the function return the region this cell belongs to
   */
  virtual Cell* findCell(const point_t &pos, region_t **region = NULL);

  /*!
   * Clear all pair distances
   */
  virtual void clearAll();

  /*!
   * Clear all tags of all particle
   */
  virtual void clearTags();

  /*!
   * Return the list of all cells
   */
  virtual vector<Cell*> &cells() {
    return m_cells;
  }

  /*!
   * Return cell with index \a index
   * @param index The index of the cell to return
   */
  virtual Cell* cell(int index) {
    return m_cells[index];
  }

  /*!
   * Return the list of all cell links
   */
  virtual vector<AbstractCellLink*> &links() {
    return m_links;
  }
  
  virtual bool& sameCP() {
    return m_scp;
  }

  /*!
   * Return the number of all active links
   */
#ifdef _OPENMP
  virtual vector<size_t> activeLinks() {
    return m_n_active_links;
  }
#else   
  virtual size_t activeLinks() {
    return m_n_active_links;
  }
#endif  

  /*!
   * Return the first cell link
   * Each thread has its own firstLink (for openMP version)
   */
#ifdef _OPENMP   
  virtual vector<AbstractCellLink*> firstLink() {
    return m_first_link;
  }
#else
  virtual AbstractCellLink* firstLink() {
    return m_first_link;
  }
#endif  

  /*!
   * Return the first active cell
   */
  virtual Cell* firstCell() {
    return m_first_cell;
  }

  /*!
   * Return the number of all active cells
   */
  virtual size_t activeCells() {
    return m_n_active_cells;
  }

  /*!
   * Return the list of all color pairs
   */
  inline virtual vector<ColourPair*> &colourPairs()
  {
    return m_colourPairs;
  }

  /*!
   * Return a pointer to the parent phase object
   */
  virtual Phase *phase() {
    return m_phase;
  }

  friend class Cell;
  friend class BoundaryCell;
  friend class AbstractCellLink;
  friend class LinkedListCreator;

  /*!
   * Export the cell subdivison information to the VTK file with filename \a filename
   * @param filename File name of the VTK
   */
  virtual void toVTK(string filename);

  /*!
   * Export the cell subdivision information as VTK to the stream \a s
   * @param s Stream to write VTK information to
   */
  virtual void toVTK(ostream &s);

  /*!
   * Return the colour pair for colors \a c1 and \a c2
   * @param c1 First color
   * @param c2 Second color
   */
  virtual ColourPair *cp(size_t c1, size_t c2)
  {
    /* We have the following assigment table (for the index used for m_colourPairs):

        0 1 2 3 4 ... < c2
      0 0 1 3 6 .
      1   2 4 7 .
      2     5 8 .
      3       9 .
      4         .
      ...
      ^
      c1

      With index = Sum[i, {i, 1, c2}] + c1 = c2*(c2+1)/2+c1
    */

/*     MSG_DEBUG("ColourPair::cp", "START"); */

    size_t index;

    if (!(c1 < m_species.size() && c2 < m_species.size())) {
      throw gError("ColourPair::cp", "(c1 < m_species.size() && c2 < m_species.size()) failed: c1 = " + ObjToString(c1) + ", c2 = " + ObjToString(c2) + ", m_colourPairs.size() = " + ObjToString( m_colourPairs.size()) + ", m_species.size() = " + ObjToString(m_species.size()) );
    }

    if(c2 < c1) {
      size_t temp = c1;
      c1 = c2;
      c2 = temp;
    }

    index = c2*(c2+1)/2+c1;

    if (!(index < m_colourPairs.size())) {
      throw gError("ColourPair::cp", "index < m_colourPairs.size() failed: c1 = " + ObjToString(c1) + ", c2 = " + ObjToString(c2) + ", index = " + ObjToString(index) + ", m_colourPairs.size() = " + ObjToString( m_colourPairs.size()));
    }

/*     MSG_DEBUG("ColourPair::cp", "END"); */

    return m_colourPairs[index];
  }

  /*!
   * Return the number of registered colors
   */
  size_t nColours() const {
    return m_species.size();
  }

  /*!
   * Find the species with name \a species and return its color
   * @param species Name of the species
   */
  int findSpecies(string species) {
    int r = -1;

    for (size_t i = 0; i < m_species.size(); i++) {
      if (m_species[i] == species) {
        r = i;
        assert(r >= 0);
      }
    }

    return r;
  }

  /*!
   * Find the species with name \a species and return its color. Throws an exception
   * if the species is not found.
   * @param species Name of the species
   */
  size_t getColour(string species) {
    int r = findSpecies(species);

    if (r == -1)
      throw gError
        ("ManagerCell::getColour", "Unknown species: '" + species + "'");

    return r;
  }


  /*!
   * Return the name of the species for a certain color
   * @param c Color for which to find the name
   */
  const string &species(size_t c) const {
    return m_species[c];
  }

  friend class Integrator;
};

class ManagerCell2x: public ManagerCell
{
  public:

    /*!
     * Constructor
     * @param p Pointer to the parent phase object
     */
    ManagerCell2x(Phase* p) : ManagerCell(p, 2) {m_direct_neighbors = c_2x_direct_neighbors;};
    
    inline void cellDist2x(Cell *first, Cell *second, int alignment, int orderTarget, int orderInter);
    /*!
     * Returns total number of neighbours (direct + indirect)
     */
    virtual inline int num_neighbors() { return NUM_2X_NEIGHBORS; };

    virtual void findRegionCellIndirectNeighbors(region_t *r, bool_point_t periodic, bool oneCellPeriodicDims);
    /*!
     * Returns the 3-coord offset vector of a specified neighbour \e a.
     * @param a Index of neighbour
     */
    virtual int_point_t c_offsets(int a) { return c_2x_offsets[a].second; };
    virtual int c_2x_1x(int a) { return c_2x_offsets[a].first; };

    /*!
     * Returns offset in specified direction \e b of 3-coord offset vector of specified neighbour \e a.
     * @param a Index of neighbour
     * @param b Index of direction
     */
    virtual int c_offsets (int a, int b) { return c_2x_offsets[a].second[b]; };
    
    /*!
     * Finds and returns the index of neighbour specified by offset vector \e off
     * @param off 3-coord int offset vector
     */
    virtual inline int offset2neighbor(int_point_t off){
      MSG_DEBUG("ManagerCell2x::offset2neighbor", "");
      const int BY = m_divby;
      int n = (off[0]+BY)*(pow(BY*2+1,2)) + (off[1]+BY)*(BY*2+1) + (off[2]+BY);
      if (n > num_neighbors()/2)
        n--;
      return n;
    };
};

#endif

