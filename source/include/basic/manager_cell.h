/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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

#include "cell.h"
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



/*!
 * This structure holds information about a region in space. A region
 * is a cuboid which is in return subdivided by cells.
 */
struct region_t {
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

/*!
 * The \a ManagerCell is in charge of cell subdivision for the whole simulation
 * region.
 */
class ManagerCell
{
protected:
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

  /*!
   * List of all cell links in the system
   */
  vector<CellLink*> m_links;

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
  vector<CellLink*> m_first_link;
#else
  CellLink* m_first_link;
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
  void activateCellLink(CellLink *c);

  /*!
   * Deactivate cell link \a c, this means one of both cells has been deactivated.
   * @param c Cell link to deactivate
   */
  void deactivateCellLink(CellLink *c);

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

public:
  /*!
   * Constructor
   * @param p Pointer to the parent phase object
   */
  ManagerCell(Phase* p);

  /*!
   * Destructor
   */
  virtual ~ManagerCell();

  /*!
   * Advance the position of the particles, taking care of walls in the system
   * and set \a m_distances_valid to false 
   * @param integrator Position integrator to use
   */
  virtual void invalidatePositions(IntegratorPosition *integrator);
  
  /*!
   * Notifies \a Cell s of changed \a Particle positions which are 
   * not caused by an \a IntegratorPosition
   * @param colour Colour of the \a Particle s with changed positions
   */
  virtual void invalidatePositions(size_t colour);
  
  /*!
   * Thread-number counter
   */
  static size_t thread_counter;

  /*!
   * Create the pair list
   */
  //virtual void createDistances();

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
  static void connect(int dir, region_t *a, region_t *b, bool_point_t periodic, int t);

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
  virtual vector<CellLink*> &links() {
    return m_links;
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
  virtual vector<CellLink*> firstLink() {
    return m_first_link;
  }
#else
  virtual CellLink* firstLink() {
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
  friend class CellLink;
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


#endif

