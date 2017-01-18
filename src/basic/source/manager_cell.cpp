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



#include <algorithm>

#include "threads.h"
#include "inlet_cell.h"
#include "simulation.h"
#include "vertex_list.h"
#include "manager_cell.h"
#include "particle_cache.h"
#include "particle_creator.h"
// #include "colour_pair.h"
#include "pair_creator.h"
//#include "cell.h"
// necessary for arguments that are pointers to Phase
#include "phase.h"

using namespace std;

int ManagerCell::thread_counter = 0;

// class ColourPair;
//---- Constructors/Destructor ----

#ifdef _OPENMP
void ManagerCell::init()
{
//    MSG_DEBUG("ManagerCell::ManagerCell", "Constructor.");
  m_first_link.resize(global::n_threads);
  m_n_active_links.resize(global::n_threads);
  m_checkForFarNeighbors = m_phase->pairCreator()->checkForFarNeighbors();

  for (int t = 0; t < global::n_threads; ++t) {
    m_first_link[t]= NULL;
    m_n_active_links[t] = 0;
  }

#ifdef ENABLE_PTHREADS
  pthread_mutex_init(&m_cells__mutex, &g_mutex_attr);
  pthread_mutex_init(&m_links__mutex, &g_mutex_attr);
#endif
}

ManagerCell::ManagerCell(Phase* p)
  : m_first_cell(NULL), m_n_active_cells(0),
    /*m_n_active_links(0), m_first_link(NULL),*/
    m_phase(p)/* m_distances_valid(false), m_pairsValid(false) */,
    m_divby(1)
{ 
  init();
}

ManagerCell::ManagerCell(Phase* p, int divby)
  : m_first_cell(NULL), m_n_active_cells(0),
    /*m_n_active_links(0), m_first_link(NULL),*/
    m_phase(p), /* m_distances_valid(false), m_pairsValid(false) */,
    m_divby(divby)
{
  init();
}

#else
void ManagerCell::init()
{
#ifdef ENABLE_PTHREADS
  pthread_mutex_init(&m_cells__mutex, &g_mutex_attr);
  pthread_mutex_init(&m_links__mutex, &g_mutex_attr);
#endif
  m_checkForFarNeighbors = m_phase->pairCreator()->checkForFarNeighbors();
}

ManagerCell::ManagerCell(Phase* p)
  : m_first_cell(NULL), m_n_active_cells(0),
    m_n_active_links(0), m_first_link(NULL),
    m_phase(p),m_divby(1)
{
//    MSG_DEBUG("ManagerCell::ManagerCell", "Constructor.");
  init();
}

ManagerCell::ManagerCell(Phase* p, int divby)
  : m_first_cell(NULL), m_n_active_cells(0),
    m_n_active_links(0), m_first_link(NULL),
    m_phase(p),m_divby(divby)
{
//  MSG_DEBUG("ManagerCell::ManagerCell", "Constructor.");
  init();
}

#endif


ManagerCell::~ManagerCell()
{
//    MSG_DEBUG("ManagerCell::~ManagerCell", "Destructor.");

#ifdef ENABLE_PTHREADS
  pthread_mutex_destroy(&m_cells__mutex);
  pthread_mutex_destroy(&m_links__mutex);
#endif

  for (vector<AbstractCellLink*>::iterator i = m_links.begin(); i != m_links.end(); i++)
    delete *i;

  for (vector<Cell*>::iterator i = m_cells.begin(); i != m_cells.end(); i++)
    delete *i;

  for (list<region_t*>::iterator i = m_regions.begin(); i != m_regions.end(); i++)
    delete *i;

  for (vector<ColourPair*>::iterator i = m_colourPairs.begin(); i != m_colourPairs.end(); i++)
    delete *i;
}

inline void ManagerCell::cellDist(Cell *first, Cell *second, int alignment, int orderTarget, int orderInter)
{
  point_t dist;
  assert(orderInter == 0);
  int s = 0;
  for (int_point_t off = c_offsets(alignment); s < SPACE_DIMS; s++)
    if (off[s] == 1)
      dist[s] = first->corner2[s] - first->corner1[s];
    else if (off[s] == -1)
      dist[s] = second->corner1[s] - second->corner2[s]; //dist[s] = -(second->corner2[s] - second->corner1[s]);
    else
      dist[s] = 0;

  first->m_cell_dist[alignment][orderTarget] = dist;
  second->m_cell_dist[num_neighbors()-alignment-1][orderTarget] = -dist;
}

inline void ManagerCell2x::cellDist2x(Cell *first, Cell *second, int alignment, int orderTarget, int orderInter)
{
  pair<int, int> interAlStep = interMap[alignment];

  Cell *inter = first->m_outlets[interAlStep.first][orderInter];
  //TODO either prev or next
  first->m_cell_dist[alignment][orderTarget] = first->m_cell_dist[m_direct_neighbors[interAlStep.first]][orderInter] + inter->m_cell_dist[m_direct_neighbors[interAlStep.second]][(orderTarget)?1-orderInter:0];
  second->m_cell_dist[num_neighbors()-alignment-1][orderTarget] = -first->m_cell_dist[alignment][orderTarget];
}

void ManagerCell::invalidatePositions(IntegratorPosition *integrator)
{
/* --- Checking if adaptive cell subdivision works! DEBUG CODE! --- */

/*
  MSG_DEBUG("ManagerCell::invalidatePositons", "Checking if adaptive cell subdivision works.");

  for (Cell *cur = m_first_cell; cur; cur = cur->next) {
    if (!(cur->hasParticles() || cur->hasFrozenParticles())) {
      throw gError("ManagerCell::invalidatePositions", "Cell without particles in cell list!");
    }
  }

  FOR_EACH
    (vector<Cell*>,
     m_cells,

     if ((*i)->hasParticles() || (*i)->hasFrozenParticles()) {
       bool found = false;
       for (Cell *cur = m_first_cell; cur; cur = cur->next) {

         if (cur->identify() == (*i)->identify())
           found = true;
       }

       if (!found) {
         throw gError("ManagerCell::invalidatePositions", "Cell with particles not in cell list!");
       }
     }
    );

  MSG_DEBUG("ManagerCell::invalidatePositons", "Seems to work.");
*/

/* --- Check done. ------------------------------------------------ */


//  int size = m_cells.size();


  /* Advance positions of all particles. This is part of the
     integration algorithm of the respective Integrator algorithm. */
  LL_FOR_EACH__PARALLEL
    (Cell,
     m_first_cell,
     m_n_active_cells,
     integrator,
//     i->updatePositions((IntegratorPosition*) data);
     ((BoundaryCell *) i)->updatePositions((IntegratorPosition*) data);
    );

  /*Tell the ParticleCreators to create additional particles if they plan to*/
  phase()->boundary()->createMoreParticles();

  /* All particles that left cell in the above step are in an
     additional commit list in the new cells to avoid updating
     them twice. Insert them into the real particle lists now. */
  FOR_EACH__PARALLEL
    (vector<Cell*>,
     m_cells,
     NULL,
     (*__iFE)->commitInjections();
    );

}


void ManagerCell::cellSubdivisionFinished()
{
  MSG_DEBUG("ManagerCell::cellSubdivisionFinished", m_links.size() << " cell links established.");
  MSG_DEBUG("ManagerCell::cellSubdivisionFinished", "APPROX. MEM USED. cells: " << m_cells.size() * (sizeof(Cell) + num_neighbors() * 20) << ", links: " << m_links.size() * sizeof(AbstractCellLink));
}


Cell *ManagerCell::findCell(const point_t &pos, region_t **region)
{
  /* Loop over all regions and see whether pos is in region or not.
     When within a region the cell can be found directly because all
     cells in a region have the same size.
  */
  for (list<region_t*>::iterator r = m_regions.begin(); r != m_regions.end(); r++) {
    if ((*r)->isInside(pos)) {
      int_point_t p;

      p.x = (int) ((*r)->inv_width.x * (pos.x - (*r)->corner1.x));
      p.y = (int) ((*r)->inv_width.y * (pos.y - (*r)->corner1.y));
      p.z = (int) ((*r)->inv_width.z * (pos.z - (*r)->corner1.z));

      Cell* cell = (*r)->cellByPos(p);

      if (cell) {
        if (!cell->isInsideEps(pos, g_geom_eps)) {
          cout << "pos = " << pos << ", p = " << p << ", " << TOCELLINDEX(p, (*r)->n_cells) << ", inv_width = " << (*r)->inv_width << endl;
          cout << cell->tag << " " << TOCELLINDEX(cell->tag, (*r)->n_cells) << " c1 = " << cell->corner1 << ", c2 = " << cell->corner2 << endl;
          cout << "r->corner1 = " << (*r)->corner1 << endl;
          throw gError("ManagerCell::findCell", "FATAL: Point is not inside. ");
        }

        if (region)
          *region = *r;
      }

      return cell;
    }
  }

  return NULL;
}


#if 0
#define ODD(x)  (x & 1)

bool ManagerCell::isInside(const point_t &pos)
{
  /*    Cell *c = findCell(pos);

  if (c)
  return m_phase->boundary()->isInside(pos);
  else
  return false;*/


  line_t tester;
  int c_up, c_down, direction;
  Cell *c;
  region_t *r;
  int_point_t cpos, cur_pos;
  cuboid_t bounding_box;

  bounding_box = m_phase->boundary()->boundingBox();

  /* Is point in box bounding the total geometry?
   */
  if (!bounding_box.isInside(pos)) {
    return false;
  }


  /* Send line from current position to a position
     outside the box and count how many walls are cut.
     Odd -> Particle is inside
     Even -> Particle is outside
  */
  tester.from = pos;
  tester.to = pos;

  /* Don't check into a direction that is periodic.
   */
  if(m_periodic.x) {
    if(m_periodic.y) {
      // all periodic or z NOT periodic => take z
      direction = 2;
      //			tester.to.z = m_bounding_box.corner2.z+10;
    } else {
      direction = 1;
      //			tester.to.y = m_bounding_box.corner2.y+10;
    }
  } else {
    direction = 0;
    //		tester.to.x = m_bounding_box.corner2.x+10;
  }

  /* Count intersections. Loop over cells because they know by which
     walls they are intersected. Speeds up test (hopefully).
   */
  c = findCell(pos, &r);
  cpos = c->tag;
  cur_pos = cpos;

  c_up = 0;
  c_down = 0;

  tester.to[direction] = bounding_box.corner2[direction] + 10;
  tester.from[direction] = pos[direction] - g_geom_eps;

  for (int i = cpos[direction]; i < r->n_cells[direction]; ++i) {
    cur_pos[direction] = i;

    c = r->cellByPos(cur_pos);

    if (c)
      c_up += c->intersections(tester, direction);
  }

  tester.to[direction] = bounding_box.corner1[direction] - 10;
  tester.from[direction] = pos[direction] + g_geom_eps;

  for (int i = cpos[direction]; i >= 0; --i) {
    cur_pos[direction] = i;

    c = r->cellByPos(cur_pos);

    if (c)
      c_down += c->intersections(tester, direction);
  }

  /* Up and down check for consistency. */


  /* Note: The assertion fails when the point sits on a surface,
     so we're just checking if one of the test got a positive
     result. fixme!!! */
  //    assert(ODD(c_up/*+add*/) == ODD(c_down/*+add*/));

  if (ODD(c_up/*+add*/) || ODD(c_down/*+add*/)) {
    return true;
  }

  return false;
}

#undef ODD
#endif

void ManagerCell::clearAll()
{
  for (vector<Cell*>::iterator i = m_cells.begin(); i != m_cells.end(); i++)
  {
    for(size_t colour = 0; colour < nColours(); ++colour)
    {
      (*i)->particles(colour).clear();
      (*i)->frozenParticles(colour).clear();
    }
  }
}


void ManagerCell::clearTags()
{
  FOR_EACH__PARALLEL
    (vector<Cell*>,
     m_cells,
     NULL,
     (*__iFE)->clearTags();
    );
}


/**
 *
 * @param container
 */
void ManagerCell::assignContainer(WallContainer *container)
{
  MSG_DEBUG("ManagerCell::assignContainer", "Sorting walls to cells.");

  /* Assign this container to every cell. The cell is then looking for
     the walls crossing them and keeps track of these walls.
  */
  FOR_EACH__PARALLEL
    (vector<Cell*>,
     m_cells,
     container,
     (*__iFE)->assignContainer((WallContainer*) data);
  );

//   MSG_DEBUG("ManagerCell::assignContainer", "setupWalls");

  /* setupWalls will determine which walls cross the *neighboring*
     cells and keep track of them as well. Needed for the test if
     a particle that leaves a cell actually did hit a wall.
  */
  FOR_EACH__PARALLEL
    (vector<Cell*>,
     m_cells,
     NULL,
     (*__iFE)->setupWalls();
  );
}

/* cell subdivision tools */

/* Look for neighbors (which means periodicities) and connect
   cells appropriately
*/
void ManagerCell::findRegionCellDirectNeighbors(region_t *r, bool_point_t periodic, AbstractCellLink *(Cell::*ap_addNeighbor)(Cell*, int, bool))
{ 
  for (vector<Cell*>::iterator i = r->cells.begin(); i != r->cells.end(); i++)  {
    for (int j = 0, n; j < NUM_HALF_DIRECT_NEIGHBORS; j++) {
      n = m_direct_neighbors[j];
//      MSG_DEBUG("ManagerCell::findRegionCellDirectNeighbors", j << " -> " << n);
      if ((*i)->neighbors(n).empty()) {
      //if (!((*i)->neighbors(n)[0])) {
	/* No neighbor on this side? Find one! */

	int_point_t neighbor;

//        MSG_DEBUG("ManagerCell::", "c_offsets("<<n<<") = " << this->c_offsets(n));
	bool neighbourFound = true;

	for (int dim = 0; dim < SPACE_DIMS; dim++) {
	  neighbor[dim] = (*i)->tag[dim] + c_offsets(n, dim);

	  if (periodic[dim])
	    neighbor[dim] = (neighbor[dim] + r->n_cells[dim]) % r->n_cells[dim];

	  if (r->n_cells[dim] == 1 && c_offsets(n,dim) != 0) neighbourFound = false;
        }
//
//	MSG_DEBUG("ManagerCell::findRegionCellDirectNeighbors", "this = " << (*i)->tag);
	if (neighbor.x >= 0 && neighbor.x < r->n_cells.x &&
	    neighbor.y >= 0 && neighbor.y < r->n_cells.y &&
	    neighbor.z >= 0 && neighbor.z < r->n_cells.z) {// && !(neighbor == (*i)->tag)) {

//	  MSG_DEBUG("ManagerCell::findRegionCellDirectNeighbors", "neighbor = " << neighbor );
	  if (Cell *c = r->cellByPos(neighbor)){
            cellDist(*i, c, n, 0);
            if (neighbourFound){
              m_links.push_back(((*i)->*ap_addNeighbor)(c, n, false));
            }else{ 
              (*i)->addPeriodic(c, n);
              int inv_direct_neighbor = num_neighbors() - n - 1;
              (*i)->m_cell_dist[inv_direct_neighbor][0] = -(*i)->m_cell_dist[n][0];
              (*i)->addPeriodic(c, inv_direct_neighbor);
            }
          }
	}
      }
    }
  }
}

/* Look for neighbors (which means periodicities) and connect
   cells appropriately
*/
void ManagerCell2x::findRegionCellIndirectNeighbors(region_t *r, bool_point_t periodic, AbstractCellLink *(Cell::*p_establishLink)(Cell*, int, bool, bool, bool))
{
  for (vector<Cell*>::iterator i = r->cells.begin(); i != r->cells.end(); i++)  {
    for (int j = 0, n=0; j < NUM_HALF_2X_INDIRECT_NEIGHBORS;j++) {
      n= c_indirect_neighbors[j];
//      MSG_DEBUG("ManagerCell::findRegionCellInDirectNeighbors", j << " -> " << n);
      int_point_t off = c_offsets(n);
//      if ((*i)->neighbors(n).empty()) {
//      No neighbor on this side? Find one!

      int_point_t neighbor;

      //MSG_DEBUG("ManagerCell::cellSubdivide", "in neighbour loop");
      bool neighbourFound = true;
//      bool directNeighbour = true;

      for (int dim = 0; dim < SPACE_DIMS; dim++) {
        neighbor[dim] = (*i)->tag[dim] + off[dim];

        if (periodic[dim])
          neighbor[dim] = (neighbor[dim] + r->n_cells[dim]) % r->n_cells[dim];

        if(r->n_cells[dim] == 1 && off[dim] != 0) neighbourFound = false;
//        if(c_offsets(n,dim) > 1 || c_offsets(n,dim) < -1) directNeighbour = false;
      }

//	MSG_DEBUG("ManagerCell::findRegionCellIndirectNeighbors", "this = " << (*i)->tag);
      if (neighbourFound &&
          neighbor.x >= 0 && neighbor.x < r->n_cells.x &&
          neighbor.y >= 0 && neighbor.y < r->n_cells.y &&
          neighbor.z >= 0 && neighbor.z < r->n_cells.z) {

        if (Cell *c = r->cellByPos(neighbor)){
//	  MSG_DEBUG("ManagerCell::findRegionCellIndirectNeighbors", "neighbor = " << neighbor );
          cellDist2x(*i, c, n, 0, 0);
            m_links.push_back(((*i)->*p_establishLink)(c, n, true, true, false));
        }
      }
    }
  }
}

void ManagerCell::addBoundaryCell(region_t *r, cuboid_t cuboid, int_point_t cell_pos, int cellIndex, int group)
{
  BoundaryCell *bc = new BoundaryCell(this, cuboid.corner1, cuboid.corner2, cell_pos, r, group);
  r->cells.push_back(bc);
  r->cells_by_pos[cellIndex] += r->cells.size();
}

region_t *ManagerCell::cellSubdivide
(double cutoff, point_t corner1, point_t corner2, bool_point_t periodic,
 int group, ParticleCreator *pc, int axis, int dir, int action)
{
  point_t d, width;
  region_t *r = new region_t;

  MSG_DEBUG
    ("ManagerCell::cellSubdivide",
     "corner1 = " << corner1 << ", corner2 = " << corner2);

  MSG_DEBUG("ManagerCell::cellSubdivide", "cutoff = " << cutoff);
  MSG_DEBUG("ManagerCell::cellSubdivide", "m_divby = " << m_divby);
  MSG_DEBUG("ManagerCell::cellSubdivide", "periodic = " << periodic);
  MSG_DEBUG("ManagerCell::cellSubdivide", "pc = " << (pc != NULL));

  /* Setup the region information */
  r->corner1 = corner1;
  r->corner2 = corner2;

  d = corner2 - corner1;

  //    m_cutoff *= 2;

  /* Determine the number of cells that fit in x, y and z direction
   */
  for (int i = 0; i < SPACE_DIMS; i++) {
    // 2010-05-04: the if-else attempts to avoid the abort when no cutoff is needed
    if(cutoff > 0)
      r->n_cells[i] = (int) (d[i] / (cutoff / m_divby));
    if (cutoff == 0 || r->n_cells[i] <= (2*m_divby + 1))
      r->n_cells[i] = 1;

    r->inv_width[i] = r->n_cells[i] / d[i];
    width[i] = d[i] / r->n_cells[i];
  }

  /* Currently only works for 3 space dimensions. */
  const int xyPlane = r->n_cells.y*r->n_cells.x;
  int intervalSize = r->n_cells.z * r->n_cells.y * r->n_cells.x;
  int intervalCounter = 1;

  //r->boundaryCells.reserve(m_divby * 2 * (r->n_cells.z*(r->n_cells.y-1) + (r->n_cells.z-1) * r->n_cells.x + (r->n_cells.x-1)*r->n_cells.y));
  r->cells.reserve(intervalSize);
  r->cells_by_pos.resize(intervalSize,-2);

  bool oneCellPeriodicDims = false;
  for (int dim = 0; dim < SPACE_DIMS; dim++)//{
    oneCellPeriodicDims |= (r->m_oneCellPeriodicDims[dim] = (periodic[dim] && r->n_cells[dim] == 1));
  MSG_DEBUG("ManagerCell::cellSubdivide", "r->n_cells = " << r->n_cells);
  if (pc) {
    int_point_t axes = {(axis+1)%3, (axis+2)%3, axis};
    int_point_t cell_pos = {0,0,0}; cell_pos[axis] = (dir==1)?r->n_cells[axis]-1:0;//Note: assuming that dir==1 is always for InletCell_Create and dir==-1 is always for InletCell_Delete
    cuboid_t cuboid; cuboid.corner1[axis] = cell_pos[axis] * width[axis];

    for (cell_pos[axes[1]] = 0; cell_pos[axes[1]] < r->n_cells[axes[1]]; cell_pos[axes[1]]++, cuboid.corner1[axes[1]] += width[axes[1]])
      for (cell_pos[axes[0]] = 0; cell_pos[axes[0]] < r->n_cells[axes[0]]; cell_pos[axes[0]]++, cuboid.corner1[axes[0]] += width[axes[0]]){
        Cell *c;
        const int cellIndex = cell_pos[2]*xyPlane+cell_pos[1]*r->n_cells.y+cell_pos[0];
        cuboid.corner2 = cuboid.corner1 + width;
        if (!m_phase->smartCells() || m_phase->boundary()->isInside \
            (cuboid, m_phase -> pairCreator() -> interactionCutoff())){
          if (action == P_CREATE)//Note: assuming that dir==1 is always for InletCell_Create and dir==-1 is always for InletCell_Delete
              c = new InletCell_Create(this, pc, axis, dir, cuboid, cell_pos, r, group);
          else 
              c = new InletCell_Delete(this, pc, axis, dir, cuboid, cell_pos, r, group);
          r->cells.push_back(c);
          r->cells_by_pos[cellIndex] = r->cells.size()-1;
	  //MSG_DEBUG("ManagerCell::cellSubdivide::pc", "cellIndex = " << cellIndex);
        }else /* no cell was created. */
          r->cells_by_pos[cellIndex] = -1;
      }
  }

  for (int i=0; i<2; i++)
  {
    int_point_t pos;
    int lo=i*(r->n_cells.z-m_divby);
    int cellIndex = lo*xyPlane;
    cuboid_t cuboid; 
    for (pos.z = lo, cuboid.corner1.z = corner1.z + lo*width.z; pos.z < lo+m_divby; pos.z++, cuboid.corner1.z += width.z)
      for (pos.y = 0, cuboid.corner1.y = corner1.y; pos.y < r->n_cells.y; pos.y++, cuboid.corner1.y += width.y)
        for (pos.x = 0, cuboid.corner1.x = corner1.x; pos.x < r->n_cells.x; pos.x++, cellIndex++, cuboid.corner1.x += width.x)
          if (r->cells_by_pos[cellIndex] == -2){
            cuboid.corner2 = cuboid.corner1 + width;
            r->cells_by_pos[cellIndex] = -1;
            if (!m_phase->smartCells() || 
                m_phase->boundary()->isInside
                (cuboid, m_phase -> pairCreator() -> interactionCutoff()))
              addBoundaryCell(r, cuboid, pos, cellIndex, group);
	    //MSG_DEBUG("ManagerCell::cellSubdivide::z_planes", "cellIndex = " << cellIndex);
	    //MSG_DEBUG("ManagerCell::cellSubdivide::z_planes", "cuboid.corner1 = " << cuboid.corner1);
          }
  }
  int cellIndex, cellIndex2;
  cellIndex = m_divby*xyPlane;

  int_point_t pos;
  cuboid_t cuboid; 
  for (pos.z = m_divby, cuboid.corner1.z = corner1.z + m_divby*width.z; pos.z < r->n_cells.z-m_divby; pos.z++, cuboid.corner1.z += width.z)
  {
    for (int i=0; i<2; i++)
    {
      int lo=i*(r->n_cells.y-m_divby);
      cellIndex += (lo-m_divby*i)*r->n_cells.x;
      for (pos.y = lo, cuboid.corner1.y = corner1.y + lo*width.y; pos.y < lo+m_divby; pos.y++, cuboid.corner1.y += width.y)
        for (pos.x = 0, cuboid.corner1.x = corner1.x; pos.x < r->n_cells.x; pos.x++, cellIndex++, cuboid.corner1.x += width.x){
          if (r->cells_by_pos[cellIndex] == -2){
            cuboid.corner2 = cuboid.corner1 + width;
            r->cells_by_pos[cellIndex] = -1;
            if (!m_phase->smartCells() ||
                m_phase->boundary()->isInside
                (cuboid, m_phase -> pairCreator() -> interactionCutoff()))
              addBoundaryCell(r, cuboid, pos, cellIndex, group);
	    //MSG_DEBUG("ManagerCell::cellSubdivide::y_planes", "cellIndex = " << cellIndex);
	    //MSG_DEBUG("ManagerCell::cellSubdivide::y_planes", "cuboid.corner1 = " << cuboid.corner1);
          }
        }
    }
    cellIndex2 = pos.z*xyPlane + m_divby*r->n_cells.x;
    for (pos.y = m_divby, cuboid.corner1.y = corner1.y + m_divby * width.y; pos.y < r->n_cells.y-m_divby; pos.y++, cuboid.corner1.y += width.y)
      for (int i=0; i<2; i++)
      { 
        int lo=i*(r->n_cells.x-m_divby);
        cellIndex2 += i*(r->n_cells.x-2*m_divby);
        for (pos.x = lo, cuboid.corner1.x = corner1.x + lo*width.x; pos.x < lo+m_divby; pos.x++, cellIndex2++, cuboid.corner1.x += width.x){
          if (r->cells_by_pos[cellIndex2] == -2){
            cuboid.corner2 = cuboid.corner1 + width;
              r->cells_by_pos[cellIndex2] = -1;
            if (!m_phase->smartCells() ||
                m_phase->boundary()->isInside
                (cuboid, m_phase -> pairCreator() -> interactionCutoff()))
              addBoundaryCell(r, cuboid, pos, cellIndex2, group);
	    //MSG_DEBUG("ManagerCell::cellSubdivide::x_planes", "cellIndex2 = " << cellIndex2);
	    //MSG_DEBUG("ManagerCell::cellSubdivide::x_planes", "cuboid.corner1 = " << cuboid.corner1);
          }
        }
      }
  }

//  intervalSize /= 20;
  cellIndex = m_divby*(xyPlane + r->n_cells.x + 1);
  int_point_t cell_pos;
  cuboid.corner1 = corner1 + m_divby*width;
  for (cell_pos.z = m_divby, cuboid.corner1.z = corner1.z + m_divby*width.z; cell_pos.z < r->n_cells.z - m_divby; cell_pos.z++, cuboid.corner1[2] += width[2], cellIndex += 2*m_divby*r->n_cells.x)
    for (cell_pos.y = m_divby, cuboid.corner1.y = corner1.y + m_divby*width.y; cell_pos.y < r->n_cells.y - m_divby; cell_pos.y++, cuboid.corner1[1] += width[1], cellIndex += 2*m_divby)
      for (cell_pos.x = m_divby, cuboid.corner1.x = corner1.x + m_divby*width.x; cell_pos.x < r->n_cells.x - m_divby; cell_pos.x++, cuboid.corner1[0] += width[0],\
					     cuboid.corner2 = cuboid.corner1 + width, cellIndex++) {
	if((1+cell_pos.x)*(1+cell_pos.y)*(1+cell_pos.z) > intervalCounter*intervalSize)
	  {
	    cout << "done: " << intervalCounter*5 << "%" << endl;
	    ++intervalCounter;
	  }

	cuboid.corner2 = cuboid.corner1 + width;
        r->cells_by_pos[cellIndex] = -1;

	if (!m_phase->smartCells() ||
	    m_phase->boundary()->isInside
	    (cuboid, m_phase -> pairCreator() -> interactionCutoff())
// 	    (cuboid, ((Simulation*) (m_phase -> parent())) -> maxCutoff)
	    )
	{ 
          r->cells.push_back(new Cell(this, cuboid, cell_pos, group));
          r->cells_by_pos[cellIndex] += r->cells.size();
        }
      }

  
  AbstractCellLink *(Cell::*p_addNeighbor)(Cell*, int, bool);
  AbstractCellLink *(Cell::*p_establishLink)(Cell*, int, bool, bool, bool);
  if (!oneCellPeriodicDims){
    if (m_checkForFarNeighbors){
      p_addNeighbor = &Cell::addNeighbor<AddPairCheck_Regular, FarNeighborCheck_On>;
      p_establishLink = &Cell::establishLink<AddPairCheck_Regular, FarNeighborCheck_On>;
    }else{
      p_addNeighbor = &Cell::addNeighbor<AddPairCheck_Regular, FarNeighborCheck_Off>;
      p_establishLink = &Cell::establishLink<AddPairCheck_Regular, FarNeighborCheck_Off>;
    }
    for (vector<Cell*>::iterator i = r->cells.begin(); i != r->cells.end(); i++)
      (*i)->init<AddPairCheck_Regular>();
  }else{
    for (vector<Cell*>::iterator i = r->cells.begin(); i != r->cells.end(); i++)
      (*i)->init<AddPairCheck_OneCellDims>();
    if (m_checkForFarNeighbors){
      p_addNeighbor = &Cell::addNeighbor<AddPairCheck_OneCellDims, FarNeighborCheck_On>;
      p_establishLink = &Cell::establishLink<AddPairCheck_OneCellDims, FarNeighborCheck_On>;
    }else{
      p_addNeighbor = &Cell::addNeighbor<AddPairCheck_OneCellDims, FarNeighborCheck_Off>;
      p_establishLink = &Cell::establishLink<AddPairCheck_OneCellDims, FarNeighborCheck_Off>;
    }
  }
  // MSG_DEBUG("ManagerCell::cellSubdivide", "Looking for neighbors.");

  findRegionCellDirectNeighbors(r, periodic, p_addNeighbor);
  findRegionCellIndirectNeighbors(r, periodic, p_establishLink);
//  MSG_DEBUG("ManagerCell::cellSubdivide", r->cells.size() << " cells created.");

  /* Done... Insert into m_cells and m_regions. */
  m_cells.insert(m_cells.end(), r->cells.begin(), r->cells.end());
  m_regions.push_back(r);
//  MSG_DEBUG("ManagerCell::cellSubdivide", r->cells.size() << " AFTER. push back regions!");
//   MSG_DEBUG("ManagerCell::cellSubdivide", "m_regions.size() = " << m_regions.size());

//   LL_FOR_EACH__PARALLEL
//     (CellLink,
//      firtLink(),
//      activeLinks(),
//      NULL,
//
//
//      MSG_DEBUG("ManagerCell::cellSubdivide", i->actsOn() << "actsOnInfo")
//     );
//      MSG_DEBUG("ManagerCell::cellSubdivide", firstLink()->next->actsOn().first << "actsOnInfo");

  return r;
}

/* Fixme!!! a has to lie to the right of b. */ //migth be already fixed
void ManagerCell::connect(int adir, region_t *a, region_t *b, bool_point_t periodic, int t)
{
    int_point_t p, q, off, dir;

    if (adir >= 0) {
        dir[0] = adir;
        p[dir[0]] = a->n_cells[dir[0]]-1;
        q[dir[0]] = 0;
        off[dir[0]] = 1;
    } else {
        dir[0] = -adir-1;
        p[dir[0]] = 0;
        q[dir[0]] = b->n_cells[dir[0]]-1;
        off[dir[0]] = -1;
    }
    dir[1] = (dir[0]+1)%SPACE_DIMS;
    dir[2] = (dir[0]+2)%SPACE_DIMS;

    MSG_DEBUG("ManagerCell::connect", (t==OUTLET) << " dir " << dir[1] << ": N(a)=" << a->n_cells[dir[1]] << ", N(b)=" << b->n_cells[dir[1]]);
    MSG_DEBUG("ManagerCell::connect", (t==OUTLET) << " dir " << dir[2] << ": N(a)=" << a->n_cells[dir[2]] << ", N(b)=" << b->n_cells[dir[2]]);

    assert(a->n_cells[dir[1]] == b->n_cells[dir[1]]);
    assert(a->n_cells[dir[2]] == b->n_cells[dir[2]]);

/*    if (t == OUTLET){
        acts_on_first = false;
    }else{
        acts_on_first = true;
    }*/
    if (t == OUTLET)
      if (a->n_cells[dir[1]] != 1)
        if (m_checkForFarNeighbors)
          connectPlanePair(a, b, periodic, &Cell::addOutlet<AddPairCheck_Regular, FarNeighborCheck_On>, &ManagerCell::cellDist, 1, p, q, off, dir);
        else
          connectPlanePair(a, b, periodic, &Cell::addOutlet<AddPairCheck_Regular, FarNeighborCheck_Off>, &ManagerCell::cellDist, 1, p, q, off, dir);
      else
        if (m_checkForFarNeighbors)
          connectPlanePair(a, b, periodic, &Cell::addOutlet<AddPairCheck_OneCellDims, FarNeighborCheck_On>, &ManagerCell::cellDist, 1, p, q, off, dir);
        else
          connectPlanePair(a, b, periodic, &Cell::addOutlet<AddPairCheck_OneCellDims, FarNeighborCheck_Off>, &ManagerCell::cellDist, 1, p, q, off, dir); 
    else {
      AbstractCellLink *(Cell::*est)(Cell *, int, bool);
      if (a->n_cells[dir[1]] != 1){
        if (m_checkForFarNeighbors){
          connectPlanePair(a, b, periodic, &Cell::addNeighbor<AddPairCheck_Regular, FarNeighborCheck_On>, &ManagerCell::cellDist, m_divby, p, q, off, dir);
          est = &Cell::establishLink<AddPairCheck_Regular, FarNeighborCheck_On>;
        }else{
          connectPlanePair(a, b, periodic, &Cell::addNeighbor<AddPairCheck_Regular, FarNeighborCheck_On>, &ManagerCell::cellDist, m_divby, p, q, off, dir);
          est = &Cell::establishLink<AddPairCheck_Regular, FarNeighborCheck_On>;
        }
      }else{
        if (m_checkForFarNeighbors){
          connectPlanePair(a, b, periodic, &Cell::addNeighbor<AddPairCheck_OneCellDims, FarNeighborCheck_On>, &ManagerCell::cellDist, m_divby, p, q, off, dir);
          est = &Cell::establishLink<AddPairCheck_OneCellDims, FarNeighborCheck_On>;
        }else{
          connectPlanePair(a, b, periodic, &Cell::addNeighbor<AddPairCheck_OneCellDims, FarNeighborCheck_Off>, &ManagerCell::cellDist, m_divby, p, q, off, dir);
          est = &Cell::establishLink<AddPairCheck_OneCellDims, FarNeighborCheck_Off>;
        }
      }
      if (m_divby > 1) {
        if (adir >= 0) {
          off[dir[0]]++;
          connectPlanePair(a, b, periodic, est, &ManagerCell::cellDist2x, m_divby, p, q, off, dir, 1); //m_divby
          p[dir[0]]--;
        } else {
          off[dir[0]]--;
          connectPlanePair(a, b, periodic, est, &ManagerCell::cellDist2x, m_divby, p, q, off, dir, 1); //m_divby
          q[dir[0]]--;
        }
        connectPlanePair(a, b, periodic, est, &ManagerCell::cellDist2x, m_divby, p, q, off, dir, 0); 
      }
    }
}
/*
virtual void ManagerCell2x::ConnectTwoCells(Cell *c, region_t *b, int_point_t q, int_point_t off, bool acts_on_first)
{                         
  int where = offset2neighbor(off);

  if (Cell *d = b->cellByPos(q)){
    ManagerCell2x::cellDist(c, d, where, 1);
    c->establishLink(d, where, acts_on_first, true);
  }
}*/

void ManagerCell::connectPlanePair(region_t *a, region_t *b, bool_point_t periodic, AbstractCellLink *(Cell::*add)(Cell *, int, bool), void (ManagerCell::*dist)(Cell *, Cell*, int, int, int), int klbound, int_point_t p, int_point_t q, int_point_t off, int_point_t dir, int orderInter)
{
    int dir1 = dir[1]; int dir2 = dir[2];
    for (p[dir1] = 0; p[dir1] < a->n_cells[dir1]; p[dir1]++) {
        for (p[dir2] = 0; p[dir2] < a->n_cells[dir2]; p[dir2]++) {
            if (Cell *c = a->cellByPos(p)) {
                for (int k = -klbound; k <= klbound; k++) {
                    for (int l = -klbound; l <= klbound; l++) {
                        bool neighbourFound = true;
                        q[dir1] = p[dir1] + k;
                        q[dir2] = p[dir2] + l;

                        MSG_DEBUG("ManagerCell::connectPlanePair", q << " " << b->n_cells);
                        if (periodic[dir1])
                            q[dir1] = (q[dir1] + a->n_cells[dir1]) % a->n_cells[dir1];
                        if (periodic[dir2])
                            q[dir2] = (q[dir2] + a->n_cells[dir2]) % a->n_cells[dir2];
                        MSG_DEBUG("ManagerCell::connectPlanePair", q << " " << b->n_cells);
//deal with indirect neighbour, should call only establishlink, two meanings of word indirectoutlet, outlet over periodic end of region, another is outlet due to half cutoff cell width.
                        // if (neighbor_index >= 0 && neighbor_index < b->cells.size()) {
                        if (q.x >= 0 && q.x < b->n_cells.x &&
                            q.y >= 0 && q.y < b->n_cells.y &&
                            q.z >= 0 && q.z < b->n_cells.z){
                            off[dir1] = k; off[dir2] = l;
                            int where = offset2neighbor(off);
                            MSG_DEBUG("ManagerCell::connectPlanePair", off << " " << where);
                            if(a->n_cells[dir1] == 1 && c_offsets(where,dir1) != 0 || a->n_cells[dir2] == 1 && c_offsets(where, dir2) != 0) neighbourFound = false;

                            if (Cell *d = b->cellByPos(q)){
                              (this ->* dist)(c, d, where, 1, orderInter);
                              if (neighbourFound)
                                  m_links.push_back((c ->* add)(d, where, true));
                              else if (abs(c_offsets(where,dir1)) <= 1 && abs(c_offsets(where,dir2)) <= 1)
                                c->addPeriodic(d, where, true);
                            }
                       }
                   }
                }
            }
        }
    }
}



/*--- VTK export ---*/

void ManagerCell::toVTK(string filename)
{
    ofstream s;

    s.open(filename.c_str());
    toVTK(s);
    s.close();
}



void ManagerCell::toVTK(ostream &s)
{
    VertexList vertices;
    map<region_t*, int> first_vertex;

    MSG_DEBUG("ManagerCell::toVTK", "Computing vertex list...");

    for (list<region_t*>::iterator r = m_regions.begin(); r != m_regions.end(); r++) {
        point_t d = (*r)->corner2-(*r)->corner1;
        point_t s;

        first_vertex[*r] = vertices.vertices().size();

        for (int i = 0; i < SPACE_DIMS; i++)
            s[i] = d[i]/(*r)->n_cells[i];

        for (int k = 0; k < (*r)->n_cells.z+1; k++)
            for (int j = 0; j < (*r)->n_cells.y+1; j++)
                for (int i = 0; i < (*r)->n_cells.x+1; i++) {
                    point_t p;

                    p.x = i*s.x;
                    p.y = j*s.y;
                    p.z = k*s.z;

                    p += (*r)->corner1;

                    vertices.addVertex(p);
                }
    }

    /*
    for (vector<Cell*>::iterator i = m_cells.begin(); i != m_cells.end(); i++)  {
        point_t corners[8];

        corners[0] = (*i)->corner1;
        corners[1] = (*i)->corner1;
        corners[2] = (*i)->corner1;
        corners[4] = (*i)->corner1;

        corners[3] = (*i)->corner2;
        corners[5] = (*i)->corner2;
        corners[6] = (*i)->corner2;
        corners[7] = (*i)->corner2;

        corners[1].x = (*i)->corner2.x;
        corners[2].y = (*i)->corner2.y;
        corners[4].z = (*i)->corner2.z;

        corners[6].x = (*i)->corner1.x;
        corners[5].y = (*i)->corner1.y;
        corners[3].z = (*i)->corner1.z;

        for (int j = 0; j < 8; j++)
            vertices.needVertex(corners[j]);
    }
    */

    s << "# vtk DataFile Version 2.0" << endl;
    s << "SYMPLER simulation cell subdivision information " << endl;
    s << "ASCII" << endl;
    s << "DATASET UNSTRUCTURED_GRID" << endl;

    s << "POINTS " << vertices.vertices().size() << " double" << endl;

    MSG_DEBUG("ManagerCell::toVTK", "Writing vertex list...");
    for (vector<point_t>::iterator i = vertices.vertices().begin(); i != vertices.vertices().end(); i++) {
        s << i->x << " " << i->y << " " << i->z << endl;
    }
    s << endl;

    int n_cells = 0;
    for (vector<Cell*>::iterator i = m_cells.begin(); i != m_cells.end(); i++)  {
        n_cells ++;
    }

    MSG_DEBUG("ManagerCell::toVTK", "Writing cells...");
    s << "CELLS " << n_cells << " " << 9*n_cells << endl;
    //    for (vector<Cell*>::iterator i = m_cells.begin(); i != m_cells.end(); i++) {
    for (list<region_t*>::iterator r = m_regions.begin(); r != m_regions.end(); r++) {
        int_point_t n_vertices = (*r)->n_cells;
        int fv = first_vertex[*r];

        for (int i = 0; i < SPACE_DIMS; i++)
            n_vertices[i]++;

        for (vector<Cell*>::iterator i = (*r)->cells.begin(); i != (*r)->cells.end(); i++) {
            /*
            point_t corners[8];

            corners[0] = (*i)->corner1;
            corners[1] = (*i)->corner1;
            corners[2] = (*i)->corner1;
            corners[4] = (*i)->corner1;

            corners[3] = (*i)->corner2;
            corners[5] = (*i)->corner2;
            corners[6] = (*i)->corner2;
            corners[7] = (*i)->corner2;

            corners[1].x = (*i)->corner2.x;
            corners[2].y = (*i)->corner2.y;
            corners[4].z = (*i)->corner2.z;

            corners[6].x = (*i)->corner1.x;
            corners[5].y = (*i)->corner1.y;
            corners[3].z = (*i)->corner1.z;
            */

            int v[8];
            int_point_t p = (*i)->tag;

            /*  --- USING VOXELS ---  */
            v[0] = TOCELLINDEX(p, n_vertices)+fv;
            p.x++;
            v[1] = TOCELLINDEX(p, n_vertices)+fv;
            p.x--; p.y++;
            v[2] = TOCELLINDEX(p, n_vertices)+fv;
            p.x++;
            v[3] = TOCELLINDEX(p, n_vertices)+fv;
            p.x--; p.y--; p.z++;
            v[4] = TOCELLINDEX(p, n_vertices)+fv;
            p.x++;
            v[5] = TOCELLINDEX(p, n_vertices)+fv;
            p.x--; p.y++;
            v[6] = TOCELLINDEX(p, n_vertices)+fv;
            p.x++;
            v[7] = TOCELLINDEX(p, n_vertices)+fv;


            /* --- USING HEXAHEDRONS ---
            v[0] = TOCELLINDEX(p, n_vertices)+fv;
            p.x++;
            v[1] = TOCELLINDEX(p, n_vertices)+fv;
            p.y++;
            v[2] = TOCELLINDEX(p, n_vertices)+fv;
            p.x--;
            v[3] = TOCELLINDEX(p, n_vertices)+fv;
            p.y--; p.z++;
            v[4] = TOCELLINDEX(p, n_vertices)+fv;
            p.x++;
            v[5] = TOCELLINDEX(p, n_vertices)+fv;
            p.y++;
            v[6] = TOCELLINDEX(p, n_vertices)+fv;
            p.x--;
            v[7] = TOCELLINDEX(p, n_vertices)+fv;
            */

            s << "8";
            for (int j = 0; j < 8; j++) {
                /*
                int v = vertices.findVertex(corners[j]);
                assert(v != -1);
                s << " " << v;
                */

                s << " " << v[j];
            }
            s << endl;
        }
    }
    s << endl;

    s << "CELL_TYPES " << n_cells << endl;
    for (int i = 0; i < n_cells; i++) {
        s << "11" << endl;
    }
    s << endl;
}


/* Colour management. 
FIXME: why is this here and not, e.g., in the Phase ?!? */

size_t ManagerCell::addColour(string species)
{

  size_t c, offset;

  assert(m_cells.size() == 0);

  c = m_species.size();
  m_species.push_back(species);
  Particle::s_cached_properties.resize(m_species.size());
  for(size_t i = 0; i < Particle::s_cached_properties.size(); ++i)
    Particle::s_cached_properties[i].resize(/*PCA_MAX_STAGE+1*/3);
  Particle::s_cached_properties_0.resize(m_species.size());
  for(size_t i = 0; i < Particle::s_cached_properties_0.size(); ++i)
    Particle::s_cached_properties_0[i].resize(/*PCA_MAX_STAGE+1*/3);

  MSG_DEBUG("ManagerCell::addColour", "species = " << species << ", c = " << c);

  m_colourPairs.resize((c+1)*(c+2)/2, NULL);

  offset = c*(c+1)/2;
  for (size_t i = 0; i < c+1; i++) {
    ColourPair *cp = new ColourPair(this);

//     MSG_DEBUG
//       ("ManagerCell::addColour",
//        "offset = " << offset << ", " <<
//        "i = " << i << ", c = " << c << ", " <<
//        "m_species[i] = " << m_species[i] << ", m_species[c] = " << m_species[c]);

    cp->m_1stColour = i;
    cp->m_2ndColour = c;

    m_colourPairs[offset+i] = cp;
    cp->posInArr() = offset+i;
  }

  Particle::s_tag_format.resize(c+1);
  Particle::s_tag_format[c].setFormatAndAlloc(new DataFormat());

  m_phase->addColour(c);

  assert(m_cells.empty());

/*
  FOR_EACH
    (vector<Cell*>,
     m_cells,
     (*i)->addColour(c);
    );
*/
//   MSG_DEBUG(
//     "ManagerCell::addColour", "m_colourPairs.size() = " << m_colourPairs.size() << ", Particle::s_tag_format.size() == " << Particle::s_tag_format.size());

  return c;
}


ColourPair *ManagerCell::cp(string firstS, string secondS)
{
  int c1 = -1, c2 = -1;

//   MSG_DEBUG("ManagerCell::cp", "firstS = " << firstS << ", secondS = " << secondS);

  for (size_t i = 0; i < m_species.size(); i++) {
    if (m_species[i] == firstS) {
      c1 = i;
      assert(c1 > -1);
    }
    if (m_species[i] == secondS) {
      c2 = i;
      assert(c2 > -1);
    }
  }

  if (c1 == -1) {
    MSG_DEBUG("ManagerCell::cp", "Creating colour '" << firstS << "'.");

    c1 = addColour(firstS);

    if (firstS == secondS)
      c2 = c1;
  }

  if (c2 == -1) {
    MSG_DEBUG("ManagerCell::cp", "Creating colour '" << secondS << "'.");

    c2 = addColour(secondS);
  }

  MSG_DEBUG("ManagerCell::cp", "c1 = " << c1 << ", c2 = " << c2);

  return cp(c1, c2);
}



/* --- Cell activation, deactivation --- */

void ManagerCell::activateCell(Cell *c)
{
#ifdef ENABLE_PTHREADS
  pthread_mutex_lock(&m_cells__mutex);
#endif

  assert(c->next == NULL && c->prev == NULL);

  c->next = m_first_cell;
  if (c->next)
    c->next->prev = c;
  m_first_cell = c;
  c->prev = NULL;

  assert(m_n_active_cells < m_cells.size());
  ++m_n_active_cells;

#ifdef ENABLE_PTHREADS
  pthread_mutex_unlock(&m_cells__mutex);
#endif
}

void ManagerCell::deactivateCell(Cell *c)
{
#ifdef ENABLE_PTHREADS
  pthread_mutex_lock(&m_cells__mutex);
#endif

// MSG_DEBUG("ManagerCell::deactivateCell", "active cells before: " << m_n_active_cells);

// next line was commented out, because it could happen that only one cell was active
// and is now deactivated. The activation of cells happens later in the timestep
// in Cell::commitInjections(), so, temporarily, we have to allow for zero active cells
// see also ManagerCell::deactivateCellLink

//   assert(c->next != NULL || c->prev != NULL);

  if (c->prev) {
    c->prev->next = c->next;
  } else {
    m_first_cell = c->next;
  }
  if (c->next) {
    c->next->prev = c->prev;
  }

  c->next = c->prev = NULL;

  assert(m_n_active_cells > 0);
  --m_n_active_cells;

#ifdef ENABLE_PTHREADS
  pthread_mutex_unlock(&m_cells__mutex);
#endif
}


void ManagerCell::activateCellLink(AbstractCellLink *c)
{
#ifdef ENABLE_PTHREADS
  pthread_mutex_lock(&m_links__mutex);
#endif


#ifdef _OPENMP
  c->mThread() = thread_counter;
// MSG_DEBUG("ManagerCell::activateCellLink", "mthread = " << c->mThread());
  assert(c->next == NULL && c->prev == NULL);

  c->next = m_first_link[thread_counter];
  if (c->next)
    c->next->prev = c;
  m_first_link[thread_counter] = c;
  c->prev = NULL;
//MSG_DEBUG("ManagerCell::activateCellLink", "next = " << c->next << " size of links = " << m_links.size());
//This might make a problem in the parallel version
//  assert(m_n_active_links < m_links[thread_counter].size());
  ++m_n_active_links[thread_counter];
//  MSG_DEBUG("ManagerCell::activateCellLink", "active links for thread = " << m_n_active_links[thread_counter] << " thread counter = " << thread_counter);
 
  ++thread_counter;
  if (thread_counter == global::n_threads)
    thread_counter = 0;
    
#else    
  assert(c->next == NULL && c->prev == NULL);

  c->next = m_first_link;
//  MSG_DEBUG("ManagerCell::activateCellLink", "cell = " << c << "next = " << m_first_link);
  if (c->next)
    c->next->prev = c;
  m_first_link = c;
  c->prev = NULL;
//  MSG_DEBUG("ManagerCell::activateCellLink", "size of links = " << m_links.size());
  assert(m_n_active_links < m_links.size());
  ++m_n_active_links;
//  MSG_DEBUG("ManagerCell::activateCellLink", "active links = " << m_n_active_links);
#endif    


#ifdef ENABLE_PTHREADS
  pthread_mutex_unlock(&m_links__mutex);
#endif
}


void ManagerCell::deactivateCellLink(AbstractCellLink *c)
{
#ifdef ENABLE_PTHREADS
  pthread_mutex_lock(&m_links__mutex);
#endif

#ifdef _OPENMP
  int thread_no = c->mThread();

// next line was commented out, because it could happen that only one cellLink was active
// and is now deactivated. The activation of cellLinks happens later in the timestep
// in Cell::commitInjections(), so, temporarily, we have to allow for zero active cellLinks
// see also ManagerCell::deactivateCell

// assert(c->next != NULL || c->prev != NULL);

  if (c->prev) {
    c->prev->next = c->next;
    MSG_DEBUG("ManagerCell::deactivateCellLink", "CELL DEAKTIVATE case 1");
  } else {
    m_first_link[thread_no] = c->next;
        MSG_DEBUG("ManagerCell::deactivateCellLink", "CELL DEAKTIVATE case 2");
  }
  
  if (c->next) {
    c->next->prev = c->prev;
        MSG_DEBUG("ManagerCell::deactivateCellLink", "CELL DEAKTIVATE case 3");
  }

  c->next = c->prev = NULL; 
  
  assert(m_n_active_links[thread_no] > 0);
  --m_n_active_links[thread_no];
#else
    if (c->prev) {
    c->prev->next = c->next;
  } else {
    m_first_link = c->next;
  }
  
  if (c->next) {
    c->next->prev = c->prev;
  }

  c->next = c->prev = NULL;  

  assert(m_n_active_links > 0);
  --m_n_active_links;
#endif   

#ifdef ENABLE_PTHREADS
  pthread_mutex_unlock(&m_links__mutex);
#endif
}
