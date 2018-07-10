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

// necessary for arguments that are pointers to Phase
#include "phase.h"

using namespace std;

size_t ManagerCell::thread_counter = 0;

//---- Constructors/Destructor ----

#ifdef _OPENMP
ManagerCell::ManagerCell(Phase* p)
  : m_first_cell(NULL), m_n_active_cells(0),
    /*m_n_active_links(0), m_first_link(NULL),*/
    m_phase(p)/* m_distances_valid(false), m_pairsValid(false) */
{
  m_first_link.resize(global::n_threads);
  m_n_active_links.resize(global::n_threads);
  
  for (size_t t = 0; t < global::n_threads; ++t) {
    m_first_link[t]= NULL;
    m_n_active_links[t] = 0;
  }

}

#else
ManagerCell::ManagerCell(Phase* p)
  : m_first_cell(NULL), m_n_active_cells(0),
    m_n_active_links(0), m_first_link(NULL),
    m_phase(p)
{

}
#endif


ManagerCell::~ManagerCell()
{

  for (vector<CellLink*>::iterator i = m_links.begin(); i != m_links.end(); i++)
    delete *i;

  for (vector<Cell*>::iterator i = m_cells.begin(); i != m_cells.end(); i++)
    delete *i;

  for (list<region_t*>::iterator i = m_regions.begin(); i != m_regions.end(); i++)
    delete *i;

  for (vector<ColourPair*>::iterator i = m_colourPairs.begin(); i != m_colourPairs.end(); i++)
    delete *i;
}


void ManagerCell::invalidatePositions(size_t colour)
{
  
  /* Tell cells that positions of particles have changed. */
  LL_FOR_EACH__PARALLEL
    (Cell,
     m_first_cell,
     m_n_active_cells,
     NULL,
     i->updatePositions(colour);
    );

  /* All particles that left cell in the above step are in an
     additional commit list in the new cells to avoid updating
     them twice. Insert them into the real particle lists now. */
  FOR_EACH__PARALLEL
    (vector<Cell*>,
     m_cells,
     //NULL,
     (*__iFE)->commitInjections();
    );

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
     i->updatePositions((IntegratorPosition*) data);
    );

  /*Tell the ParticleCreators to create additional particles if they plan to*/
  // FIXME: This smells like a very dangerous call! What if this
  // function is called at a different place than where it is used to
  // be called right now? Maybe we do not wish new particles to be
  // created there? Seems as if you should disentangle that!
  phase()->boundary()->createMoreParticles();

  /* All particles that left cell in the above step are in an
     additional commit list in the new cells to avoid updating
     them twice. Insert them into the real particle lists now. */
  FOR_EACH__PARALLEL
    (vector<Cell*>,
     m_cells,
     //NULL,
     (*__iFE)->commitInjections();
    );

}


void ManagerCell::cellSubdivisionFinished()
{
  MSG_DEBUG("ManagerCell::cellSubdivisionFinished", m_links.size() << " cell links established.");
  MSG_DEBUG("ManagerCell::cellSubdivisionFinished", "APPROX. MEM USED. cells: " << m_cells.size() * (sizeof(Cell) + NUM_NEIGHBORS * 20) << ", links: " << m_links.size() * sizeof(CellLink));
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
          cout << "pos = " << pos << endl;
          cout << "c1 = " << cell->corner1 << ", c2 = " << cell->corner2 << endl;
          throw gError("ManagerCell::findCell", "FATAL: Point is not inside.");
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
     // NULL,
     (*__iFE)->clearTags();
    );
}


/*!
 * Assign given container to every cell. The cell is then looking for
 *    the walls crossing them and keeps track of these walls.
 *
 * @param container \a Container object to be assigned
 */
void ManagerCell::assignContainer(WallContainer *container)
{
  MSG_DEBUG("ManagerCell::assignContainer", "Sorting walls to cells.");

  FOR_EACH__PARALLEL
    (vector<Cell*>,
     m_cells,
     // container,
     // (*__iFE)->assignContainer((WallContainer*) data);
     (*__iFE)->assignContainer(container);
  );

//   MSG_DEBUG("ManagerCell::assignContainer", "setupWalls");

  /* setupWalls will determine which walls cross the *neighboring*
     cells and keep track of them as well. Needed for the test if
     a particle that leaves a cell actually did hit a wall.
  */
  FOR_EACH__PARALLEL
    (vector<Cell*>,
     m_cells,
     // NULL,
     (*__iFE)->setupWalls();
  );
}


/* cell subdivision tools */

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
      r->n_cells[i] = (int) (d[i] / cutoff);
    else 
      r->n_cells[i] = 2;
    // 2010-05-04: if the attempt above works, the following won't happen anymore
    if (r->n_cells[i] < 2) {

      // it follows some old code...

      /* If we have less than three cells we have to make sure
      pairs are not counted twice. */
      //      if (periodic[i])
	//        throw gError("ManagerCell::cellSubdivide", string(1, 'X'+i) + "-direction: Currently, box dimensions for periodic directions are not allowed to decrease 3*cutoff = " + ObjToString(3*cutoff));
	// 	r->n_cells[i] = 2;
      //      else {
      //        if (r->n_cells[i] <= 1)
	  //          throw gError("ManagerCell::cellSubdivide", string(1, 'X'+i) + "-direction: Currently, box dimensions for non-periodic directions are not allowed to decrease 2*cutoff = " + ObjToString(2*cutoff));
      //      }
      //      r->n_cells[i] = 2;

      // end of old code

      throw gError
	("ManagerCell::cellSubdivide",
  "Box length too small! No room, for at least two cells! Cells are always >= cutoff. Creating less than two cells is currently not allowed. Current cutoff is " + ObjToString(cutoff) + ". Does this make sense?");
    }

    r->inv_width[i] = r->n_cells[i] / d[i];
    width[i] = d[i] / r->n_cells[i];
  }

  /* Currently only works for 3 space dimensions. */

  int intervalSize = r->n_cells.z * r->n_cells.y * r->n_cells.x / 20;
  int intervalCounter = 1;

  for (int z = 0; z < r->n_cells.z; z++)
    for (int y = 0; y < r->n_cells.y; y++)
      for (int x = 0; x < r->n_cells.x; x++) {
	if((1+x)*(1+y)*(1+z) > intervalCounter*intervalSize)
	  {
	    cout << "done: " << intervalCounter*5 << "%" << endl;
	    ++intervalCounter;
	  }
	Cell *c = NULL;
	cuboid_t cuboid;
	int_point_t cell_pos;

	cell_pos.x = x;
	cell_pos.y = y;
	cell_pos.z = z;

	for (int i = 0; i < SPACE_DIMS; i++)
	  cuboid.corner1[i] = corner1[i] + cell_pos[i] * width[i];

	cuboid.corner2 = cuboid.corner1 + width;

	if (!m_phase->smartCells() ||
	    m_phase->boundary()->isInside
	    (cuboid, m_phase -> pairCreator() -> interactionCutoff())
// 	    (cuboid, ((Simulation*) (m_phase -> parent())) -> maxCutoff)
	    )
	  {

	    /* If a particle creator is given, determine which kind of
	       cell to create.

	       Right at the wall -> Inlet or Outlet
	    */

	    if (pc) {
	      if (dir == 1 && cell_pos[axis] == r->n_cells[axis]-1) {
		if (action == P_CREATE)
		  c = new InletCell_Create(this, pc, axis, dir, group);
		else
		  c = new InletCell_Delete(this, pc, axis, dir, group);
	      } else if (dir == -1 && cell_pos[axis] == 0) {
		if (action == P_CREATE)
		  c = new InletCell_Create(this, pc, axis, dir, group);
		else
		  c = new InletCell_Delete(this, pc, axis, dir, group);
	      } else
		c = new Cell(this, group);
	    } else
	      c = new Cell(this, group);

	    c->tag = cell_pos;
	    c->corner1 = cuboid.corner1;
	    c->corner2 = cuboid.corner2;
	  }

	/* If entry is NULL, cell does not exist. */
	if (c) {
	  r->cells.push_back(c);
	  r->cells_by_pos.push_back(r->cells.size()-1);
	} else
	  r->cells_by_pos.push_back(-1);
      }

  MSG_DEBUG("ManagerCell::cellSubdivide", "Initializing cells.");

  for (vector<Cell*>::iterator i = r->cells.begin(); i != r->cells.end(); i++)  {
    (*i)->init();
  }

//   MSG_DEBUG("ManagerCell::cellSubdivide", "Looking for neighbors.");

  /* Look for neighbors (which means periodicities) and connect
     cells appropriately
  */
  for (vector<Cell*>::iterator i = r->cells.begin(); i != r->cells.end(); i++)  {
    for (int n = 0; n < NUM_NEIGHBORS; n++) {
      if ((*i)->neighbors(n).empty()) {
	/* No neighbor on this side? Find one! */

	int_point_t neighbor;

	Cell *c;
//MSG_DEBUG("ManagerCell::cellSubdivide", "in neighbour loop");
	bool noNeighbour = false;

	for (int dim = 0; dim < SPACE_DIMS; dim++) {
	  neighbor[dim] = (*i)->tag[dim] + Cell::c_offsets[n][dim];

	  if (periodic[dim])
	    neighbor[dim] = (neighbor[dim] + r->n_cells[dim]) % r->n_cells[dim];
	  if(r->n_cells[dim] == 1 && Cell::c_offsets[n][dim] != 0) noNeighbour = true;
	}


	if (neighbor.x >= 0 && neighbor.x < r->n_cells.x &&
	    neighbor.y >= 0 && neighbor.y < r->n_cells.y &&
	    neighbor.z >= 0 && neighbor.z < r->n_cells.z) {

	  c = r->cellByPos(neighbor);
//MSG_DEBUG("ManagerCell::cellSubdivide", "if periodic");
	  if (c) {
	    if (noNeighbour) {
	      (*i)->addPeriodic(c, n);
	    } else {
	      (*i)->addNeighbor(c, n);
	    }
	  }
	}
      }
    }
  }

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


/* Fixme!!! a has to lie to the right of b. */
void ManagerCell::connect(int dir, region_t *a, region_t *b, bool_point_t periodic, int t)
{
    int_point_t p, q, off;
    int dir2, dir3;

    if (dir >= 0) {
        p[dir] = a->n_cells[dir]-1;
        q[dir] = 0;
        off[dir] = 1;
    } else {
        dir = -dir-1;
        p[dir] = 0;
        q[dir] = b->n_cells[dir]-1;
        off[dir] = -1;
    }

    dir2 = (dir+1)%SPACE_DIMS;
    dir3 = (dir+2)%SPACE_DIMS;

    MSG_DEBUG("ManagerCell::connect", "dir " << dir2 << ": N(a)=" << a->n_cells[dir2] << ", N(b)=" << b->n_cells[dir2]);
    MSG_DEBUG("ManagerCell::connect", "dir " << dir3 << ": N(a)=" << a->n_cells[dir3] << ", N(b)=" << b->n_cells[dir3]);

    assert(a->n_cells[dir2] == b->n_cells[dir2]);
    assert(a->n_cells[dir3] == b->n_cells[dir3]);

    for (p[dir2] = 0; p[dir2] < a->n_cells[dir2]; p[dir2]++) {
        for (p[dir3] = 0; p[dir3] < a->n_cells[dir3]; p[dir3]++) {
            Cell *c = a->cellByPos(p);

            if (c) {
                for (int k = -1; k <= 1; k++) {
                    for (int l = -1; l <= 1; l++) {
                        int where;

                        off[dir2] = k;
                        off[dir3] = l;
                        OFFSET2NEIGHBOR(off, where);

                        q[dir2] = p[dir2] + k;
                        q[dir3] = p[dir3] + l;

                        if (periodic[dir2])
                            q[dir2] = (q[dir2] + a->n_cells[dir2]) % a->n_cells[dir2];
                        if (periodic[dir3])
                            q[dir3] = (q[dir3] + a->n_cells[dir3]) % a->n_cells[dir3];


                        // if (neighbor_index >= 0 && neighbor_index < b->cells.size()) {
                        if (q.x >= 0 && q.x < b->n_cells.x &&
                            q.y >= 0 && q.y < b->n_cells.y &&
                            q.z >= 0 && q.z < b->n_cells.z) {

                            Cell *d;

                            d = b->cellByPos(q);

                            if (d) {
                                if (t == OUTLET)
                                    c->addOutlet(d, where);
                                else
                                    c->addNeighbor(d, where);
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

  assert(c->next == NULL && c->prev == NULL);

  c->next = m_first_cell;
  if (c->next)
    c->next->prev = c;
  m_first_cell = c;
  c->prev = NULL;

  assert(m_n_active_cells < m_cells.size());
  ++m_n_active_cells;

}

void ManagerCell::deactivateCell(Cell *c)
{

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

}


void ManagerCell::activateCellLink(CellLink *c)
{

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

}


void ManagerCell::deactivateCellLink(CellLink *c)
{

#ifdef _OPENMP
  int thread_no = c->mThread();

// next line was commented out, because it could happen that only one cellLink was active
// and is now deactivated. The activation of cellLinks happens later in the timestep
// in Cell::commitInjections(), so, temporarily, we have to allow for zero active cellLinks
// see also ManagerCell::deactivateCell

// assert(c->next != NULL || c->prev != NULL);

  if (c->prev) {
    c->prev->next = c->next;
    // MSG_DEBUG("ManagerCell::deactivateCellLink", "CELL DEACTIVATE case 1");
  } else {
    m_first_link[thread_no] = c->next;
        // MSG_DEBUG("ManagerCell::deactivateCellLink", "CELL DEACTIVATE case 2");
  }
  
  if (c->next) {
    c->next->prev = c->prev;
        // MSG_DEBUG("ManagerCell::deactivateCellLink", "CELL DEACTIVATE case 3");
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

}
