/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
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



using namespace std;

#include <algorithm>

#include "cell.h"
#include "phase.h"
#include "pairdist.h"
#include "pair_list.h"
#include "controller.h"
#include "manager_cell.h"
#include "wall_triangle.h"
#include "colour_pair.h"


// const int_point_t Cell::c_offsets[NUM_NEIGHBORS] = {
//   {{-1},{-1},{-1}}, {{-1},{-1}, {0}}, {{-1},{-1}, {1}}, {{-1}, {0},{-1}}, {{-1}, {0}, {0}}, {{-1}, {0}, {1}},
//   {{-1}, {1},{-1}}, {{-1}, {1}, {0}}, {{-1}, {1}, {1}},
//   { {0},{-1},{-1}}, { {0},{-1}, {0}}, { {0},{-1}, {1}}, { {0}, {0},{-1}}, /*{{0}, {0}, {0}},*/ { {0}, {0}, {1}},
//   { {0}, {1},{-1}}, { {0}, {1}, {0}}, { {0}, {1}, {1}},
//   { {1},{-1},{-1}}, { {1},{-1}, {0}}, { {1},{-1}, {1}}, { {1}, {0},{-1}}, { {1}, {0}, {0}}, { {1}, {0}, {1}},
//   { {1}, {1},{-1}}, { {1}, {1}, {0}}, { {1}, {1}, {1}}
// };

const int_point_t Cell::c_offsets[NUM_NEIGHBORS] = {
    {-1,-1,-1}, {-1,-1, 0}, {-1,-1, 1}, {-1, 0,-1}, {-1, 0, 0}, {-1, 0, 1},
    {-1, 1,-1}, {-1, 1, 0}, {-1, 1, 1},
    { 0,-1,-1}, { 0,-1, 0}, { 0,-1, 1}, { 0, 0,-1}, /*{0, 0, 0},*/ { 0, 0, 1},
    { 0, 1,-1}, { 0, 1, 0}, { 0, 1, 1},
    { 1,-1,-1}, { 1,-1, 0}, { 1,-1, 1}, { 1, 0,-1}, { 1, 0, 0}, { 1, 0, 1},
    { 1, 1,-1}, { 1, 1, 0}, { 1, 1, 1}
};


// const int_point_t Cell::c_half_offsets[NUM_HALF_NEIGHBORS] = {
//     {-1,-1, 1}, {-1, 0, 1},
//     {-1, 1, 0}, {-1, 1, 1},
//     {0,-1, 1}, {0, 0, 1},
//     {0, 1, 0}, {0, 1, 1},
//     {1,-1, 1}, {1, 0, 0}, {1, 0, 1},
//     {1, 1, 0}, {1, 1, 1}
// };


/*---- Class CellLink ----*/

CellLink::CellLink()
{
  throw gError("CellLink::CellLink(default)", "Should not be called. Contact the programmer.");

  set(NULL, NULL, -1, true, true);
}


CellLink::CellLink(Cell *first, Cell *second, int alignment, bool acts_on_first, bool acts_on_second)
  : m_n_active_cells(0), next(NULL), prev(NULL)/*, m_linkUsed(false)*/
{
  set(first, second, alignment, acts_on_first, acts_on_second);
}


CellLink::CellLink(const CellLink &copy)
  : m_n_active_cells(copy.m_n_active_cells), next(NULL), prev(NULL)/*, m_linkUsed(false)*/
{
  throw gError("CellLink::CellLink(copy)", "Should not be called. Contact the programmer.");

  set(copy.m_first, copy.m_second, copy.m_alignment, copy.m_acts_on.first, copy.m_acts_on.second);
}


CellLink::~CellLink()
{
}



//---- Methods ----

/* cellDist determines the distance between two cells in real coordinates.
   i.e.
   if first is left of second it takes
            second.corner1 - first.corner2
   and otherwise
            first.corner1 - second.corner2
*/
inline void cellDist(Cell *first, Cell *second, int alignment, point_t &dist)
{
  int_point_t off;
  point_t width1, width2;

  off = Cell::c_offsets[alignment];
  width1 = first->corner2 - first->corner1;
  width2 = second->corner2 - second->corner1;

  for (int s = 0; s < SPACE_DIMS; s++) {
    if (off[s] == 1)
      dist[s] = width1[s];
    else if (off[s] == -1)
      dist[s] = -width2[s];
    else
      dist[s] = 0;
  }

}


void CellLink::set(Cell *first, Cell *second, int alignment, bool acts_on_first, bool acts_on_second)
{
  m_first = first;
  m_second = second;
  m_alignment = alignment;
  m_acts_on.first = acts_on_first;
  m_acts_on.second = acts_on_second;

  assert(first && second && ((alignment != -1) || (first == second)) );

  if (first == second) {
    for (int s = 0; s < SPACE_DIMS; s++)
      m_cell_dist[s] = 0;
  } else {
    cellDist(m_first, m_second, m_alignment, m_cell_dist);
  }
}


/* One of the two cells connected by a cell link is activated.
   If both links are activated, activate this link
*/
void CellLink::cellActivated() {

  ++m_n_active_cells;

  if (!(m_n_active_cells >= 0 && m_n_active_cells <= 2)) {
    MSG_DEBUG
      ("CellLink::cellActivated",
       "FATAL INTERNAL ERROR: m_n_active_cells = " << m_n_active_cells << " is an invalid value at this point.");
    abort();
  }

// Assigning the activated CellLink in the proper thread list
  if (m_n_active_cells == 2) 
    m_first->manager()->activateCellLink(this);

}


/* One of the two cells connected by a cell link is deactivated.
   If link is active, deactivate link.
*/
void CellLink::cellDeactivated() {

  --m_n_active_cells;

  if (!(m_n_active_cells >= 0 && m_n_active_cells <= 2)) {
    MSG_DEBUG
      ("CellLink::cellDeactivated",
       "FATAL INTERNAL ERROR: m_n_active_cells = " << m_n_active_cells << " ! Aborting.");
    abort();
  }

  if (m_n_active_cells == 1)
    m_first->manager()->deactivateCellLink(this);

}


/* ---- Inline helper functions ---- */

#ifdef _OPENMP
inline void createDistancesForSame
  (vector<PairList> &distances,
   double cutoff_sq,
   int dir,
   Cell *c,
   list<Particle*> &p,
   point_t &cell_dist,
   int thread_no)
{
#else
inline void createDistancesForSame
  (vector<PairList> &distances,
   double cutoff_sq,
   int dir,
   Cell *c,
   list<Particle*> &p,
   point_t &cell_dist)
{
#endif
  list<Particle*>::iterator p_end = p.end();

  for (list<Particle*>::iterator i = p.begin(); i != p_end; ++i) {
    list<Particle*>::iterator j = i;
    for (++j; j != p_end; ++j) {

#ifdef _OPENMP
      addPair
        (distances, cutoff_sq,
         dir, c, c,
         *i, *j,
         true, true,
         cell_dist,
         thread_no);
#else
      addPair
        (distances, cutoff_sq,
         dir, c, c,
         *i, *j,
         true, true,
         cell_dist);
#endif         
    }
  }
}

#ifdef _OPENMP
inline void createDistancesForDifferent
  (vector<PairList> &distances,
   double cutoff_sq,
   int dir,
   Cell *first_c,
   Cell *second_c,
   list<Particle*> &first_p,
   list<Particle*> &second_p,
   bool ao_f,
   bool ao_s,
   point_t &cell_dist, 
   int thread_no)
{
#else
inline void createDistancesForDifferent
  (vector<PairList> &distances,
   double cutoff_sq,
   int dir,
   Cell *first_c,
   Cell *second_c,
   list<Particle*> &first_p,
   list<Particle*> &second_p,
   bool ao_f,
   bool ao_s,
   point_t &cell_dist)
{	
#endif


  list<Particle*>::iterator p1_end = first_p.end();
  list<Particle*>::iterator p2_end = second_p.end();

  for (list<Particle*>::iterator i = first_p.begin(); i != p1_end; ++i) {
    for (list<Particle*>::iterator j = second_p.begin(); j != p2_end; ++j) {

#ifdef _OPENMP       
      addPair
        (distances, cutoff_sq,
         dir, first_c, second_c,
         *i, *j,
         ao_f, ao_s,
         cell_dist,
         thread_no);
#else
      addPair
        (distances, cutoff_sq,
         dir, first_c, second_c,
         *i, *j,
         ao_f, ao_s,
         cell_dist);
#endif         
    }
  }
}


/* Creating the pair distances, using the cell linked-list method */

#ifdef _OPENMP
void CellLink::createDistances(int thread_no)
{
// MSG_DEBUG("CellLink::createDistances", "start");
 ManagerCell *manager = m_first->manager();
  size_t n_colours = manager->nColours();
  /* Are we calculating pairs within the same cell? */
  if (m_first == m_second) {
    /* Loop over the first colour */
    for (size_t c1 = 0; c1 < n_colours; ++c1) {
      /* --- Pairs of the same colour -------------------------------------------------------------- */
      ColourPair *cp = manager->cp(c1, c1);

      if (cp->needPairs()) {
	double cutoff_sq = cp->cutoff() * cp->cutoff();

	createDistancesForSame
	  (cp->freePairs(),
	   cutoff_sq,
	   0,
	   m_first,
	   m_first->particles(c1),
	   m_cell_dist,
	   thread_no);

	createDistancesForDifferent
	  (cp->frozenPairs(),
	   cutoff_sq,
	   0,
	   m_first, m_first,
	   m_first->particles(c1), m_first->frozenParticles(c1),
	   true, false,
	   m_cell_dist,
	   thread_no);
      }

      for (size_t c2 = c1+1; c2 < n_colours; ++c2) {
	/* --- Pairs of different colour -------------------------------------------------------------- */
	cp = manager->cp(c1, c2);

	if (cp->needPairs()) {
	  double cutoff_sq = cp->cutoff() * cp->cutoff();

	  createDistancesForDifferent
	    (cp->freePairs(),
	     cutoff_sq,
	     0,
	     m_first, m_first,
	     m_first->particles(c1), m_first->particles(c2),
	     true, true,
	     m_cell_dist,
	     thread_no);

	  createDistancesForDifferent
	    (cp->frozenPairs(),
	     cutoff_sq,
	     0,
	     m_first, m_first,
	     m_first->particles(c1), m_first->frozenParticles(c2),
	     true, false,
	     m_cell_dist,
	     thread_no);

	  createDistancesForDifferent
	    (cp->frozenPairs(),
	     cutoff_sq,
	     0,
	     m_first, m_first,
	     m_first->frozenParticles(c1), m_first->particles(c2),
	     false, true,
	     m_cell_dist,
	     thread_no);
	}
      } /* Loop over c2 */
    } /* Loop over c1 */
  } else { /* We have different cells. */

    for (size_t c1 = 0; c1 < n_colours; ++c1) {
      for (size_t c2 = 0; c2 < n_colours; ++c2) {
	ColourPair *cp = manager->cp(c1, c2);

	if (cp->needPairs()) {
	  double cutoff_sq = cp->cutoff() * cp->cutoff();

	  if (c1 < c2) {

	    createDistancesForDifferent
	      (cp->freePairs(),
	       cutoff_sq,
	       1,
	       m_first, m_second,
	       m_first->particles(c1), m_second->particles(c2),
	       m_acts_on.first, m_acts_on.second,
	       m_cell_dist,
	       thread_no);
	    if (m_acts_on.first)
	      createDistancesForDifferent
		(cp->frozenPairs(),
		 cutoff_sq,
		 1,
		 m_first, m_second,
		 m_first->particles(c1), m_second->frozenParticles(c2),
		 true, false,
		 m_cell_dist,
		 thread_no);
	    if (m_acts_on.second)
	      createDistancesForDifferent
		(cp->frozenPairs(),
		 cutoff_sq,
		 1,
		 m_first, m_second,
		 m_first->frozenParticles(c1), m_second->particles(c2),
		 false, true,
		 m_cell_dist,
		 thread_no);
	  } else {

     createDistancesForDifferent
	      (cp->freePairs(),
	       cutoff_sq,
	       -1,
	       m_second, m_first,
	       m_second->particles(c2), m_first->particles(c1),
	       m_acts_on.second, m_acts_on.first,
	       m_cell_dist,
	       thread_no);

	    if (m_acts_on.first)
	      createDistancesForDifferent
		(cp->frozenPairs(),
		 cutoff_sq,
		 -1,
		 m_second, m_first,
		 m_second->frozenParticles(c2), m_first->particles(c1),
		 false, true,
		 m_cell_dist,
		 thread_no);

	    if (m_acts_on.second)
	      createDistancesForDifferent
		(cp->frozenPairs(),
		 cutoff_sq,
		 -1,
		 m_second, m_first,
		 m_second->particles(c2), m_first->frozenParticles(c1),
		 true, false,
		 m_cell_dist,
		 thread_no);
	  }
	}
      }
    }
  }
}
#else
void CellLink::createDistances()
{

  ManagerCell *manager = m_first->manager();
  size_t n_colours = manager->nColours();
  /* Are we calculating pairs within the same cell? */
  if (m_first == m_second) {
    /* Loop over the first colour */
    for (size_t c1 = 0; c1 < n_colours; ++c1) {
      /* --- Pairs of the same colour -------------------------------------------------------------- */
      ColourPair *cp = manager->cp(c1, c1);

      if (cp->needPairs()) {
	double cutoff_sq = cp->cutoff() * cp->cutoff();

	createDistancesForSame
	  (cp->freePairs(),
	   cutoff_sq,
	   0,
	   m_first,
	   m_first->particles(c1),
	   m_cell_dist);

	createDistancesForDifferent
	  (cp->frozenPairs(),
	   cutoff_sq,
	   0,
	   m_first, m_first,
	   m_first->particles(c1), m_first->frozenParticles(c1),
	   true, false,
	   m_cell_dist);
      }

      for (size_t c2 = c1+1; c2 < n_colours; ++c2) {
	/* --- Pairs of different colour -------------------------------------------------------------- */
	cp = manager->cp(c1, c2);

	if (cp->needPairs()) {
	  double cutoff_sq = cp->cutoff() * cp->cutoff();

	  createDistancesForDifferent
	    (cp->freePairs(),
	     cutoff_sq,
	     0,
	     m_first, m_first,
	     m_first->particles(c1), m_first->particles(c2),
	     true, true,
	     m_cell_dist);

	  createDistancesForDifferent
	    (cp->frozenPairs(),
	     cutoff_sq,
	     0,
	     m_first, m_first,
	     m_first->particles(c1), m_first->frozenParticles(c2),
	     true, false,
	     m_cell_dist);

	  createDistancesForDifferent
	    (cp->frozenPairs(),
	     cutoff_sq,
	     0,
	     m_first, m_first,
	     m_first->frozenParticles(c1), m_first->particles(c2),
	     false, true,
	     m_cell_dist);
	}
      } /* Loop over c2 */
    } /* Loop over c1 */
  } else { /* We have different cells. */

    for (size_t c1 = 0; c1 < n_colours; ++c1) {
      for (size_t c2 = 0; c2 < n_colours; ++c2) {
	ColourPair *cp = manager->cp(c1, c2);

	if (cp->needPairs()) {
	  double cutoff_sq = cp->cutoff() * cp->cutoff();

	  if (c1 < c2) {

	    createDistancesForDifferent
	      (cp->freePairs(),
	       cutoff_sq,
	       1,
	       m_first, m_second,
	       m_first->particles(c1), m_second->particles(c2),
	       m_acts_on.first, m_acts_on.second,
	       m_cell_dist);
	    if (m_acts_on.first)
	      createDistancesForDifferent
		(cp->frozenPairs(),
		 cutoff_sq,
		 1,
		 m_first, m_second,
		 m_first->particles(c1), m_second->frozenParticles(c2),
		 true, false,
		 m_cell_dist);
	    if (m_acts_on.second)
	      createDistancesForDifferent
		(cp->frozenPairs(),
		 cutoff_sq,
		 1,
		 m_first, m_second,
		 m_first->frozenParticles(c1), m_second->particles(c2),
		 false, true,
		 m_cell_dist);
	  } else {

     createDistancesForDifferent
	      (cp->freePairs(),
	       cutoff_sq,
	       -1,
	       m_second, m_first,
	       m_second->particles(c2), m_first->particles(c1),
	       m_acts_on.second, m_acts_on.first,
	       m_cell_dist);

	    if (m_acts_on.first)
	      createDistancesForDifferent
		(cp->frozenPairs(),
		 cutoff_sq,
		 -1,
		 m_second, m_first,
		 m_second->frozenParticles(c2), m_first->particles(c1),
		 false, true,
		 m_cell_dist);

	    if (m_acts_on.second)
	      createDistancesForDifferent
		(cp->frozenPairs(),
		 cutoff_sq,
		 -1,
		 m_second, m_first,
		 m_second->particles(c2), m_first->frozenParticles(c1),
		 true, false,
		 m_cell_dist);
	  }
	}
      }
    }
  }
}
#endif

/*---- Class Cell ----*/

//---- Constructors/Destructor ----

Cell::Cell(ManagerCell *mgr, int group)
  : m_particles(mgr->nColours()), m_frozen_particles(mgr->nColours()),
    m_injected_particles(mgr->nColours()), m_n_particles(0), m_manager(mgr), m_group(group), m_cellUsed(false),
    next(NULL), prev(NULL)
{
}


Cell::Cell(ManagerCell *mgr, const point_t &c1, const point_t &c2, int group)
  : m_particles(mgr->nColours()), m_frozen_particles(mgr->nColours()),
    m_injected_particles(mgr->nColours()),
    m_n_particles(0),
    m_manager(mgr), m_group(group), m_cellUsed(false),
    next(NULL), prev(NULL)
{
  corner1 = c1;
  corner2 = c2;
}


Cell::~Cell()
{
}



//--- Methods ---

void Cell::establishLink(Cell *neighbor, int where, bool first, bool second)
{
  list<CellLink*>::iterator cur = m_neighbors[where].begin();

  while (cur != m_neighbors[where].end() && (*cur)->other(this) != neighbor)
    ++cur;

  if (cur == m_neighbors[where].end()) {
    CellLink *link;

    link = new CellLink(this, neighbor, where, first, second);

    m_neighbors[where].push_back(link);
    neighbor->m_neighbors[INV_NEIGHBOR(where)].push_back(link);

    /* Links are being deleted when the manager gets out of scope. */
    m_manager->m_links.push_back(link);
  }
}

void Cell::addNeighbor(Cell *neighbor, int where)
{
    establishLink(neighbor, where, true, true);

    if (find(m_outlets[where].begin(), m_outlets[where].end(), neighbor) == m_outlets[where].end())
        m_outlets[where].push_back(neighbor);

    where = INV_NEIGHBOR(where);

    if (find(neighbor->m_outlets[where].begin(), neighbor->m_outlets[where].end(), this)
        == neighbor->m_outlets[where].end())
        neighbor->m_outlets[where].push_back(this);
}


void Cell::addOutlet(Cell *neighbor, int where)
{
    establishLink(neighbor, where, false, true);
    //    establishLink(neighbor, where, true, false);

    if (find(m_outlets[where].begin(), m_outlets[where].end(), neighbor) == m_outlets[where].end())
        m_outlets[where].push_back(neighbor);
}

void Cell::addPeriodic(Cell *neighbor, int where)
{
    if (find(m_outlets[where].begin(), m_outlets[where].end(), neighbor) == m_outlets[where].end())
        m_outlets[where].push_back(neighbor);
}

// #define TRACK_PARTICLE 399

void Cell::doCollision(Particle *p, point_t& r, point_t& v, const point_t &force, IntegratorPosition *integratorP) {
  double t_travelled;
  point_t hit_pos, old_r;
  Wall *wall;
  bool hit = true;
  
#ifdef TRACK_PARTICLE

  if (p->mySlot == TRACK_PARTICLE)
    cout << "Cell::doCollision: tracking particle " << TRACK_PARTICLE << " !!! " << m_all_walls.size() << " walls in cell." << endl;

#endif
  size_t iterations = 0;
  while (hit) {      
    ++iterations;
    if(iterations > 100)
      throw gError("Cell::doCollision", "More than 100 wall collisions for particle " + ObjToString(p->mySlot) + ", " + ObjToString(p->r) + ", colour = " + ObjToString(p->c) + "!!!\nCheck your settings (forces, size of timestep, ...)!");
    hit = false;
    t_travelled = HUGE_VAL;
    old_r = r;

    /* First, find hit that occurs first. */
    // this cell knows about all walls, intersecting this and its neighbouring cells
    checkForHit(p, force, hit, t_travelled, hit_pos, wall, integratorP);

#ifdef TRACK_PARTICLE

    if (p->mySlot == TRACK_PARTICLE && hit) {
      cout << "--> hit @ " << hit_pos << " after time " << t_travelled << " on wall:" << endl;
      cout << wall->toString() << endl;
    }

#endif
    // this checks for the indirect outlets = "far neighbours" due to, e.g., periodic BCs
    for (list< pair<int, Cell*> >::iterator c = m_indirect_outlets.begin();
         c != m_indirect_outlets.end(); c++) {
      point_t dist;
      bool new_hit = false;

      cellDist(this, c->second, c->first, dist);


      /* Update position in case this is periodic, etc. */
      r = old_r + c->second->corner1 - corner1 - dist;

      c->second->checkForHit(p, force, new_hit, t_travelled, hit_pos, wall, integratorP);

#ifdef TRACK_PARTICLE

      if (p->mySlot == TRACK_PARTICLE && hit) {
        cout << "--> new hit @ " << hit_pos << " after time " << t_travelled << " on wall:" << endl;
        cout << wall->toString() << endl;
      }

#endif

      if (new_hit) {

        hit = true;
        hit_pos = hit_pos - c->second->corner1 + corner1 + dist;
//         MSG_DEBUG("Cell::doCollision", "new_hit = true for " << p -> r << ", " << p->mySlot << ", " << p->c);
      }
    }

    if (hit) {

//       MSG_DEBUG("Cell::doCollision", "hit = true for " << p->r << ", " << p->mySlot << ", " << p->c << "wall = " << wall->toString());
#ifdef TRACK_PARTICLE

      if (p->mySlot == TRACK_PARTICLE) {
        cout << "--> reflecting" << endl;
        cout << "before: p->r = " << p->r << ", p->v = " << /*p->*/v << ", p->dt = " << p->dt << endl;
      }

#endif

// MSG_DEBUG("Cell::doCollision", "HITCASE: wall normal = " << wall->normal());
      /* Hit happened inside. Reflect particle. */
      wall->reflector()->reflect(p, r, v, hit_pos, wall->normal(), wall->inPlane());
      p->dt -= t_travelled;
      // next is for the case that the reflector has aborted the collisions
      // by setting p->dt = 0
      if(p->dt < 0) p->dt = 0;

#ifdef TRACK_PARTICLE

      if (p->mySlot == TRACK_PARTICLE) {
        cout << "after: p->r = " << p->r << ", p->v = " << /*p->*/v << ", p->dt = " << p->dt << endl;
      }

#endif

    } else
      r = old_r;

  } // end of while(hit)
}

// Check of position updates not caused by integrators
void Cell::updatePositions(size_t colour)
{
  list<Particle*>::iterator m_particles_end = m_particles[colour].end();
  for (list<Particle*>::iterator i = m_particles[colour].begin(); i != m_particles_end; ) {

    i = checkNewPosition(i);
  }
}

/* Main integration and collision detection logic. */
void Cell::updatePositions(IntegratorPosition *integrator)
{
  Phase *phase = m_manager->phase();

  double dt;
  int force_index;
  size_t colour;

  dt = integrator->dt();
  colour = integrator->colour();

  list<Particle*>::iterator m_particles_end = m_particles[colour].end();
  for (list<Particle*>::iterator i = m_particles[colour].begin(); i != m_particles_end; ) {

    Particle *p = *i;

    p->dt = dt;

    integrator->integratePosition(p, this);
    integrator->integrateVelocity(p);

    i = checkNewPosition(i);
  }
}

// Checks if new position of given Particle requires finding a new Cell
// Return value: valid iterator to next particle
   list<Particle*>::iterator Cell::checkNewPosition(list<Particle*>::iterator i) {

  Phase *phase = m_manager->phase();
     
    Particle *p = *i;

    size_t colour = p->c;

    bool erase = false;
    bool total_erase = false;

    int_point_t off = {0, 0, 0};
    int n;
    
    // Last modified: 2007-12-27: changed from isInside to isInsideEps due to problem with geometrical epsilons. There was an inconsistency due to a small epsilon between the isInside of the cell the particle is leaving and isInside of the cell the particle should enter (see below)
    if (!isInsideEps(p->r, g_geom_eps)) {

      
      // if((*i)->mySlot == 399) MSG_DEBUG("Cell::updatePositions", "particle 399 left cell");

      for (int j = 0; j < SPACE_DIMS; j++)
        if (p->r[j] < corner1[j])
          off[j] = -1;
        else if (p->r[j] >= corner2[j])
          off[j] = 1;

      OFFSET2NEIGHBOR(off, n);

      if (m_outlets[n].empty()) {
        erase = true;
        total_erase = true;
      } else {
        point_t old_r = p->r;
        Particle *new_p = p;

        for (list<Cell*>::iterator c = m_outlets[n].begin(); c != m_outlets[n].end(); c++) {
          point_t dist;

          cellDist(this, *c, n, dist);

          if (!new_p) {
            new_p = phase->addParticle(*p);
          }

          /* Update position in case this is periodic, etc. */
          new_p->r = old_r + (*c)->corner1 - corner1 - dist;
          if ((*c)->isInside(new_p->r)) {
            (*c)->injectFree(colour, new_p);
            erase = true;

            /* Notify any subclass that a particle has left the cell in a certain direction. */
            particleLeftCell(off, n);
          } else {

	    Controller* controller = ((Simulation*) phase->parent())->controller();
	    size_t force_index = controller->forceIndex();
	    
            MSG_DEBUG("Cell::updatePositions", "Particle did not enter neighboring cell:");

	    cout << "=== DEBUGING INFORMATION: PARTICLE ===" << endl;
            cout << "particle id (slot)               = " << p->mySlot << endl;
            cout << "position before update/collision = " << old_r << endl;
            cout << "current position                 = " << p->r << endl;
	    cout << "current velocity                 = " << p->v << endl;
            cout << "force                            = " << p->force[force_index] << endl;
	    cout << "=== tag ===" << endl;
	    cout << p->tag.toString();

	    cout << "=== DEBUGING INFORMATION: CELLS ===" << endl;
            cout << "current cell: corner 1                   = " << corner1 << endl;
	    cout << "current cell: corner 2                   = " << corner2 << endl;
            cout << "entering cell: corner 1                  = " << (*c)->corner1 << endl;
	    cout << "entering cell: corner 2                  = " << (*c)->corner2 << endl;
            cout << "offset between cells (current->entering) = " << off << endl;
	    cout << "neigbor index                            = " << n << endl;
            cout << "distance between cells                   = " << dist << endl;

            throw gError("Cell::updatePositions", "Particle flew farther than "
                         "a cell in one timestep. Please check the correctness of your forces, decrease dt or "
                         "increase the cut-off radius.", gError::PARTICLEFLEWTOOFAR);
          }

          new_p = NULL;
        }
      }
    }

//     if((*i)->mySlot == 52)
//     {
//       if(phase->boundary()->isInside((*i)->r))
//         MSG_DEBUG("Cell::updatePositionsEND", "52 INSIDE");
//       else MSG_DEBUG("Cell::updatePositionsEND", "52 NOT INSIDE");
//     }

    if (erase) {
/*      if((*i)->mySlot == 52) MSG_DEBUG("Cell::updatePositions", "52: erase = TRUE");*/

      Particle *p = *i;
      list<Particle*>::iterator e = i;
      i++;
      m_particles[colour].erase(e);

      assert(m_n_particles != 0);

      m_n_particles--;

      if (m_n_particles == 0)
        deactivate();

      if (total_erase) {
	
                // cout << ">>> DELETING PARTICLE " << p->mySlot << " >>>" << endl;
                // cout << "off = " << off << endl;
                // cout << "corner1 = " << corner1 << endl;
                // cout << "corner2 = " << corner2 << endl;
                // cout << "p->r = " << p->r << endl;
                // cout << "p->v = " << p->v << endl;
                // cout << "p->dt = " << p->dt << endl;
                // cout << "p->force[0] = " << p->force[0] << endl;
                // cout << "p->force[1] = " << p->force[1] << endl;

        m_manager->phase()->removeParticle(p);
      }
    } else
      i++;

    return i;
    }
 
void Cell::init()
{
  m_local_link = new CellLink(this, this, -1);

  m_manager->m_links.push_back(m_local_link);
}

void Cell::injectFree(size_t colour, Particle *p)
{
  if (p->g != m_group)  {
    m_manager->phase()->groupChanged(p, m_group);
    p->g = m_group;
  }

  m_injected_particles[colour].push_back(p);
}

void Cell::injectFrozen(size_t colour, Particle *p)
{
  p->g = m_group;
  p->isFrozen = 1;
  m_frozen_particles[colour].push_back(p);

  if (!m_n_particles) {
    activate();
  }
  ++m_n_particles;
}


void Cell::commitInjections()
{
  size_t np = 0;

  for(size_t colour = 0; colour < m_particles.size(); ++colour)
  {
    np += m_injected_particles[colour].size();

    FOR_EACH
      (list<Particle*>,
       m_injected_particles[colour],
       m_particles[colour].push_back(*__iFE);
      );

    m_injected_particles[colour].clear();

  }

  if (np && !m_n_particles)
    activate();
  m_n_particles += np;
}


void Cell::clearTags()
{
  // Phase *phase = m_manager->phase();
  for (size_t c = 0; c < m_particles.size(); ++c) {
    FOR_EACH
      (list<Particle*>,
       m_particles[c],
       (*__iFE)->tag.clear();
       );
    
  }
}


void Cell::clearOutlets(int where)
{
  m_outlets[where].clear();
}

#if 0
void Cell::clearNeighbors(int where)
{
  for (list<CellLink*>::iterator i = m_neighbors[where].begin();
       i != m_neighbors[where].end(); i++) {
    vector<CellLink*>::iterator help;

    // Fixme!!! Slow as crap.
    help = find(m_manager->m_links.begin(), m_manager->m_links.end(), *i);

    if (help != m_manager->m_links.end()) {
      m_manager->m_links.erase(help);
    }

    //    m_links.erase(i->first);
  }

  m_neighbors[where].clear();
}
#endif

void Cell::eraseParticle(list<Particle*>::iterator e)
{
  Particle *p = *e;

  m_particles[p->c].erase(e);

  assert(m_n_particles != 0);

  m_n_particles--;

  if (m_n_particles == 0)
    deactivate();

  m_manager->phase()->removeParticle(p);
}


void Cell::assignContainer(WallContainer *container)
{
  for (list<Wall*>::iterator i = container->walls().begin(); i != container->walls().end(); i++) {
    if ((*i)->intersects(*this)) {
      addWall(*i);

      // WallTriangle *wt = (WallTriangle*) *i;
      // if(
      // 	 corner1.x > -1. && corner1.x < 1. &&
      // 	 corner1.y > 2. && corner1.y < 3. &&
      // 	 corner1.z > 103. && corner1.z < 104. 
      // 	 )
      // 	MSG_DEBUG("Cell::assignContainer", "wall ACCEPTED: " << "corner0=" << wt->corner(0) << ", corner1=" << wt->corner(1) << ", corner2=" << wt->corner(2) << ", cell:  " << corner1 << " " << corner2);

    }
    else {
      
      // WallTriangle *wt = (WallTriangle*) *i;
      // if(
      // 	 corner1.x > -1. && corner1.x < 1. &&
      // 	 corner1.y > 2. && corner1.y < 3. &&
      // 	 corner1.z > 103. && corner1.z < 104. 
      // 	 )
      // 	MSG_DEBUG("Cell::assignContainer", "wall REJECTED: " << "corner0=" << wt->corner(0) << ", corner1=" << wt->corner(1) << ", corner2=" << wt->corner(2) << ", cell:  " << corner1 << " " << corner2);
            
    }
  }
}


void Cell::setupWalls()
{
//   MSG_DEBUG("Cell::setupWalls", "START");

  copy(m_walls.begin(), m_walls.end(), inserter(m_all_walls, m_all_walls.begin()));

   // size_t counter = 0;
   // for (list<Wall*>::iterator ii = m_all_walls.begin(); ii != m_all_walls.end(); ii++) ++counter;
   // if(corner1.x > -1. && corner1.x < 1. &&
   //    corner1.y > 2. && corner1.y < 3. &&
   //    corner1.z > 103. && corner1.z < 104. 
   //    )
   //  MSG_DEBUG("Cell::setupWalls", "elements in m_all_walls before outlet-loop = " << counter << ", NUM_NEIGHBORS= " << NUM_NEIGHBORS << ", m_walls.size()=" << m_walls.size());

  /* Loop over all outlets. */
  for (int n = 0; n < NUM_NEIGHBORS; n++) {
    for (list<Cell*>::iterator i = m_outlets[n].begin(); i != m_outlets[n].end(); i++) {

      point_t dist;

      /* Check whether cell is a direct or an indirect outlet. */
      cellDist(this, *i, n, dist);

      if ((corner1+dist - (*i)->corner1).abs() < g_geom_eps) {

//       if(
// 	 corner1.x > 19.5 && corner1.x < 19.52 &&
// 	 corner1.y > -0.1 && corner1.y < 0.1 &&
// 	corner1.z > 6.59 && corner1.z < 6.6 &&
// 	corner2.x > 20.5 && corner2.x < 20.6 &&
// 	corner2.y > 1.2 && corner2.y < 1.3 &&
// 	corner2.z > 7.6 && corner2.z < 7.8
// 		       )
//       MSG_DEBUG("Cell::setupWalls", "direct outlet: " << (*i)->corner1 << ", " << (*i)->corner2);

        for (list<Wall*>::iterator j = (*i)->m_walls.begin(); j != (*i)->m_walls.end(); j++) {
          if (find(m_all_walls.begin(), m_all_walls.end(), *j) == m_all_walls.end())
            m_all_walls.push_back(*j);
        }
      } 
      else {

// 	if(
// 	   corner1.x > 19.5 && corner1.x < 19.52 &&
// 	   corner1.y > -0.1 && corner1.y < 0.1 &&
// 	  corner1.z > 6.59 && corner1.z < 6.6 &&
// 	  corner2.x > 20.5 && corner2.x < 20.6 &&
// 	  corner2.y > 1.2 && corner2.y < 1.3 &&
// 	  corner2.z > 7.6 && corner2.z < 7.8
// 			 )
// 	MSG_DEBUG("Cell::setupWalls", "INdirect outlet: " << (*i)->corner1 << ", " << (*i)->corner2);

        m_indirect_outlets.push_back(pair<int, Cell*>(n, *i));
      }
    }
  }
}


bool Cell::hasParticles() const
{
  bool has = false;

  FOR_EACH_CONST
    (vector< list<Particle*> >,
     m_particles,
     if (!i->empty())
       has = true;
    );

  return has;
}


bool Cell::hasFrozenParticles() const
{
  bool has = false;

  FOR_EACH_CONST
    (vector<list<Particle*> >,
     m_frozen_particles,
     if (!i->empty())
       has = true;
    );

  return has;
}



typedef map<Cell*, CellLink*> my_map_t;
void Cell::activate()
{
  /* Update cell list. */
  m_manager->activateCell(this);

  /* Notify local cell link. */
  m_local_link->cellActivated();
  m_local_link->cellActivated();

  for (int j = 0; j < NUM_NEIGHBORS; ++j) {
    FOR_EACH
      (list<CellLink*>,
       m_neighbors[j],
       (*__iFE)->cellActivated();
       );
  }
}


void Cell::deactivate()
{
  /* Update cell list. */
  m_manager->deactivateCell(this);

  /* Notify cell links. */
  m_local_link->cellDeactivated();
  m_local_link->cellDeactivated();

  for (int j = 0; j < NUM_NEIGHBORS; ++j) {
    FOR_EACH
      (list<CellLink*>,
       m_neighbors[j],
       (*__iFE)->cellDeactivated();
       );
  }
}
