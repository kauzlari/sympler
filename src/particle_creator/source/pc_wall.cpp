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



#include "pc_wall.h"

#include "phase.h"
#include "boundary.h"
#include "simulation.h"
#include "manager_cell.h"

#define M_BOUNDARY ((Boundary*) m_parent)
#define M_PHASE ((Phase*) M_BOUNDARY->parent())
#define M_MANAGER M_PHASE->manager()
#define M_SIMULATION ((Simulation*) M_PHASE->parent())

/* Register this ParticleCreator with the factory. */
const ParticleCreator_Register<ParticleCreatorWall> particle_creator_wall("ParticleCreatorWall");


/* ---- WallCalculator ---- */

class WallCalculator: public calculator
{
protected:
point_t *m_position, *m_offset, *m_size;

void AnalyzeId(const string &id) {
static map<char, int> space_coords;
	
if (id == "x" || id == "y" || id == "z") {
curr_tok.key = DOUBLE_PTR;
curr_tok.x_ptr = &((*m_position)[id[0]-'x']);
} else if (id == "offX" || id == "offY" || id == "offZ") {
curr_tok.key = DOUBLE_PTR;
curr_tok.x_ptr = &(*m_offset)[id[3]-'X'];
} else if (id == "sizeX" || id == "sizeY" || id == "sizeZ") {
curr_tok.key = DOUBLE_PTR;
curr_tok.x_ptr = &(*m_size)[id[4]-'X'];
} else
throw gError("WallCalculator::AnalyzeId", "Id undefined. " \
		"Possible ids are: x, y, z, velX, velY, velZ");
}
	
public:
WallCalculator(): calculator(), m_position(NULL), m_offset(NULL), m_size(NULL) {
}
WallCalculator(point_t *position, point_t *offset, point_t *size)
: calculator(), m_position(position), m_offset(offset), m_size(size) {
}
WallCalculator(const WallCalculator& old)
: calculator(), m_position(old.m_position), m_offset(old.m_offset), m_size(old.m_size) {
FromString(old.text);
}
WallCalculator(point_t *position, point_t *offset, point_t *size, const string& old)
		: calculator(), m_position(position), m_offset(offset), m_size(size) {
FromString(old);
}
};

/* --- */


bool ParticleCreatorWall::already = false;

//---- Constructor/Destructor ----

ParticleCreatorWall::ParticleCreatorWall(Boundary *boundary): ParticleCreator(boundary)
{
if(already) throw gError("ParticleCreatorWall::ParticleCreatorWall: There can be only one.");
		already = true;
		init();
}


ParticleCreatorWall::~ParticleCreatorWall()
{
already = false;
}

//---- Methods ----

// This PC does not change point_t &size
void ParticleCreatorWall::adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend)
{
	bool dd_defined = false;
	if (m_distance != HUGE_VAL)
	{
		if (m_distance <= 0)
		throw gError("ParticleCreatorWall::adjustBoxSize", "Distance must be bigger than zero.");
	
		MSG_DEBUG("ParticleCreatorWall::adjustBoxSize", "distance = " << m_distance);
	
		dd_defined = true;
	}

	if (m_density != HUGE_VAL)
	{
		if (dd_defined)
		{
		//next is the relative difference between given distance and the one computed 
		//from the density
			double diff = fabs(1 - m_distance/pow(1/m_density, (1./3)));
			MSG_DEBUG("ParticleCreatorWall::adjustBoxSize", "diff = " << diff);
			//we only complain if diff > 0.1%
			if(diff > 0.001)
				throw gError("ParticleCreatorWall::adjustBoxSize", "You have specified both " \
				"distance and density and the distance computed from the density differs more " \
				"than 0.1% from the defined one. Aborting.");
			else MSG_DEBUG("ParticleCreatorWall::adjustBoxSize", "You have specified both " \
			"distance and density and they seem to be consistent. Taking the density and computing the distance from it.");
		}
	
		MSG_DEBUG("ParticleCreatorWall::adjustBoxSize", "density = " << m_density);
	
		m_distance = pow(1/m_density, (1./3));
		MSG_DEBUG("ParticleCreatorWall::adjustBoxSize", "m_distance = " << m_distance);
		dd_defined = true;
	}
	if (!dd_defined)
		throw gError("ParticleCreatorWall::adjustBoxSize", "Please specify one of: " \
		"distance, density");

// 	Phase *phase = (Phase*) ((Boundary*) m_parent)->parent();
// 	Simulation* sim = (Simulation*) phase -> parent();	
	point_t addFrame = {{{0, 0, 0}}};
	bool_point_t periodicityFront = ((Boundary*) m_parent)->periodicityFront();
	MSG_DEBUG("ParticleCreatorWall::adjustBoxSize", "periodicityFront = " << periodicityFront);
	bool_point_t periodicityBack = ((Boundary*) m_parent)->periodicityBack();
	MSG_DEBUG("ParticleCreatorWall::adjustBoxSize", "periodicityBack = " 
	<< periodicityBack);
	
	for (int i = 0; i < SPACE_DIMS; i++) 
	{
		frameRCfront[i] = frameRCfront[i] || !periodicityFront[i] || m_forceMin[i];
		frameRCend[i] = frameRCend[i] || !periodicityBack[i] || m_forceMax[i];
		if(frameRCfront[i]) addFrame[i] += myCutoff;
		if(frameRCend[i]) addFrame[i] += myCutoff;			
	}
	MSG_DEBUG("ParticleCreatorWall::adjustBoxSize","frameRCfront = " << frameRCfront << ", frameRCend = " << frameRCend);
	// for the # of lattice points, we have to consider the whole frame
	for (int i = 0; i < SPACE_DIMS; i++) 
	{
		int n = (int) ((size[i] + addFrame[i]) / m_distance + 0.5);
		if (n < 1)
			n = 1;
		m_nlattice_points[i] = n;
	}
}

void ParticleCreatorWall::createParticles()
{
 	MSG_DEBUG("ParticleCreatorWall::createParticles()", "start");    
	
	Phase *phase = (Phase*) ((Boundary*) m_parent)->parent();
	point_t size = ((Boundary*) m_parent)->boundingBox().size();
// MSG_DEBUG("ParticleCreatorWall::createParticles()", " size = " << size);
	point_t vel = { { { 0, 0, 0 } } };
	size_t colour = M_MANAGER->getColour/*AndAdd*/(m_species);
	
	MSG_DEBUG("ParticleCreatorWall::createParticles()",
	"Color = " << colour << ", species = " << m_species << ", "
	"maxInteractionCutoff = " << myCutoff);
	
	point_t boxOffset = phase->boundary()->boundingBox().corner1;
	point_t pos;
	size_t n = 0;
	//   __m_particle.resize(M_MANAGER->nColours());
	//   m_internal_dofs.resize(M_MANAGER->nColours());

	MSG_DEBUG("ParticleCreatorWall::createParticles()", "boxOffset = " << boxOffset);
	MSG_DEBUG("ParticleCreatorWall::createParticles()", "m_nlattice_points[0] = " << m_nlattice_points[0]);
	MSG_DEBUG("ParticleCreatorWall::createParticles()", "m_nlattice_points[1] = " << m_nlattice_points[1]);
	MSG_DEBUG("ParticleCreatorWall::createParticles()", "m_nlattice_points[2] = " << m_nlattice_points[2]);
	
  // build a list of functions corresponding to the unknown properties from m_properties.unknown()
  list<FunctionFixed> fList;
  for (map<string, string>::const_iterator i = m_properties.unknown().begin(); i != m_properties.unknown().end(); ++i) 
  {
    // the functions will be deleted automatically at the end of createParticles 
    // when the list fList is deleted
    FunctionFixed func;
    fList.push_back(func);
    FunctionFixed& lastF = fList.back(); 
    lastF.addVariable("x");
    lastF.addVariable("y");
    lastF.addVariable("z");
    lastF.addVariable("u");
    lastF.addVariable("v");
    lastF.addVariable("w");
    lastF.setOnlyExpr(i->second);
    lastF.compile();
  }
  for (int i = 0; i < m_nlattice_points[0]; i++) {
    for (int j = 0; j < m_nlattice_points[1]; j++) {
      for (int k = 0; k < m_nlattice_points[2]; k++) {
	pos.x = i*m_distance+m_distance/2;
	pos.y = j*m_distance+m_distance/2;
	pos.z = k*m_distance+m_distance/2;
	pos += boxOffset;
	// this may be used to partially create wall particles
	pos.x = m_posX(pos.x, pos.y, pos.z);
	pos.y = m_posY(pos.x, pos.y, pos.z);
	pos.z = m_posZ(pos.x, pos.y, pos.z);
	
	if (!(((Boundary*) m_parent)->isInside(pos))&&
	    ((Boundary*) m_parent)->isInWallRange( myCutoff, pos) )
	  { 
	    Particle *p;
	    if(m_freeze)
	      p = phase->addFrozenParticle(pos, vel, 0, colour);
	    
	    else 
	      p = phase->addParticle(pos, vel, 0, colour);
	    
	    p->v.x = m_velX(p->r.x, p->r.y, p->r.z, p->v.x, p->v.y, p->v.z);
	    p->v.y = m_velY(p->r.x, p->r.y, p->r.z, p->v.x, p->v.y, p->v.z);
	    p->v.z = m_velZ(p->r.x, p->r.y, p->r.z, p->v.x, p->v.y, p->v.z);
	    
	    list<FunctionFixed>::iterator fListIt = fList.begin();
	    for (map<string, string>::const_iterator i = m_properties.unknown().begin(); i != m_properties.unknown().end(); ++i) 
	      {
		assert(fListIt != fList.end());
		// MSG_DEBUG("ParticleCreatorWall::createParticles", "checking for unknown property " << i->first);
		// new style treatment of unknown properties
		double* value;
		DataFormat::datatype_t type;
		size_t offset;
		string name = i->first;        
		if(Particle::s_tag_format[colour].attrExists(name)) 
		  {
		    
		    offset = Particle::s_tag_format[colour].attrByName(name).offset;
		    type = Particle::s_tag_format[colour].attrByName(name).datatype;
		    
		    if(type == DataFormat::DOUBLE)
		      {
			// value is a reference
			value = &(p->tag.doubleByOffset(offset));
		      }
		    else 
		      throw gError("ParticleCreatorWall::createParticles", "For species '" + m_species + "': Scalar-type definition used for non-scalar property '" + name + "'. For vectors, you must add \"__x\", \"__y\" or \"__z\" and for tensors you must add \"__xx\", \"__xy\", \"__xz\", \"__yx\", \"__yy\", \"__yz\", \"__zx\", \"__zy\" or \"__zz\" to the property name.");
		  }
		else
		  { 
		    if(name[name.size()-2] == '_' && name[name.size()-3] == '_')
		      {
			string symbol = name;
			if(Particle::s_tag_format[colour].attrExists(symbol.erase(symbol.size()-3, symbol.size()-1)))
			  type = Particle::s_tag_format[colour].attrByName(symbol).datatype;
			else 
			  throw gError("ParticleCreatorWall::createParticles", "For species '" + m_species + ": Found a vector extension for attribute '" + name + "' but did not find symbol '" + symbol + "'.");
			if(type != DataFormat::POINT)
			  if(type != DataFormat::TENSOR)
			    throw gError("ParticleCreatorWall::createParticles", "For species '" + m_species + ": You use a vector-extension in attribute '" + name + "', but '" + symbol + "' is not a vector.");
			offset = Particle::s_tag_format[colour].attrByName(symbol).offset;
			
			switch(name[name.size()-1])
			  {
			  case 'x':
			    value = &(p->tag.pointByOffset(offset).x); 
			    break;
			  case 'y':
			    value = &(p->tag.pointByOffset(offset).y); 
			    break;
			  case 'z':
			    value = &(p->tag.pointByOffset(offset).z); 
			    break;
			  default:
			    throw gError("ParticleCreatorWall::createParticles", "For species '" + m_species + ": Unknown vector extension '__" + name[name.size()-1] + "' in attribute '" + name + "'. Possibilities are \"__x\", \"__y\" or \"__z\".");
			  }
		      } // end if(vector-case)
		    else
		      {
			if(name[name.size()-3] == '_' && name[name.size()-4] == '_')
			  {
			    string symbol = name;
			    if(Particle::s_tag_format[colour].attrExists(symbol.erase(symbol.size()-4, symbol.size()-1)))
			      type = Particle::s_tag_format[colour].attrByName(symbol).datatype;
			    else 
			      throw gError("ParticleCreatorWall::createParticles", "For species '" + m_species + ": Found a tensor extension for attribute '" + name + "' but did not find symbol '" + symbol + "'.");
			    if(type != DataFormat::TENSOR)
			      throw gError("ParticleCreatorWall::createParticles", "For species '" + m_species + ": You use a tensor-extension in attribute '" + name + "', but '" + symbol + "' is not a tensor.");
			    offset = Particle::s_tag_format[colour].attrByName(symbol).offset;
			    
			    bool error = false;
			    switch(name[name.size()-2])
			      {
			      case 'x':
				switch(name[name.size()-1])
				  {
				  case 'x':
				    value = &(p->tag.tensorByOffset(offset).xx);
				    break;
				  case 'y':
				    value = &(p->tag.tensorByOffset(offset).xy);
				    break;
				  case 'z':
				    value = &(p->tag.tensorByOffset(offset).xz);
				    break;
				  default:
				    error = true; 
				  }
				break;
			      case 'y':
				switch(name[name.size()-1])
				  {
				  case 'x':
				    value = &(p->tag.tensorByOffset(offset).yx);
				    break;
				  case 'y':
				    value = &(p->tag.tensorByOffset(offset).yy);
				    break;
				  case 'z':
				    value = &(p->tag.tensorByOffset(offset).yz);
				    break;
				  default:
				    error = true; 
				  }
				break;
			      case 'z':
				switch(name[name.size()-1])
				  {
				  case 'x':
				    value = &(p->tag.tensorByOffset(offset).zx);
				    break;
				  case 'y':
				    value = &(p->tag.tensorByOffset(offset).zy);
				    break;
				  case 'z':
				    value = &(p->tag.tensorByOffset(offset).zz);
				    break;
				  default:
				    error = true; 
				  }
				break;
			      }
			    if(error)
			      throw gError("ParticleCreatorWall::createParticles", "For species '" + m_species + ": Unknown tensor extension '__" + name[name.size()-2] + name[name.size()-1] + "' in attribute '" + name + "'. Possibilities are \"__xx\", \"__xy\", \"__xz\", \"__yx\", \"__yy\", \"__yz\", \"__zx\", \"__zy\" or \"__zz\".");
			    
			  } // end if(tensor-case)
			else
			  {
			    MSG_DEBUG("ParticleCreatorWall::createParticles", "NOT FOUND " << i->first);
			    throw gError("ParticleCreatorWall::createParticles", 
					 "No internal degree of freedom named '" + i->first +
					 "' for species '" + m_species + "'. Possibilities are: " + Particle::s_tag_format[colour].toStr/*ing*/());
			  }
		      } // end else of if(vector-case)
		  } // end else of if(attrExists)
		
		*value = (*fListIt)(p->r.x, p->r.y, p->r.z, p->v.x, p->v.y, p->v.z);
		++fListIt;
	      }
	    //           MSG_DEBUG("ParticleCreatorWall::createParticles", n << " particles created.");
	    ++n;
	    
	  }
      }			
    }
  }
  MSG_DEBUG("ParticleCreatorWall::createParticles", n << " particles created.");
}

void ParticleCreatorWall::init()
{
  m_properties.setClassName("ParticleCreatorWall");
  m_properties.setDescription(
  "Generates particles behind walls, i.e., outside the domain,"
      " a freely moving particle can reach.\n"
      "The initial conditions of user-defined symbols can be set by taking the symbol as an attribute and defining a mathematical expression for it. For the expression, the same variables are allowed as, e.g., for the attribute 'u'. For non-scalars you have to add the following to the attribute name:\n"
      "\"_x\", \"_y\" or \"_z\" for the respective components of a vector\n"
      "\"_xx\", \"_xy\", \"_xz\", \"_yx\", \"_yy\", \"_yz\", \"_zx\", \"_zy\" and \"_zz\" for the respective components of a tensor."
  );
  /* Allow unknown properties. Those ones have to be identified later.
  They are used to set the particles degrees of freedom initially. */
  m_properties.allowUnknown();
  
  BOOLPC(forceMaxX, m_forceMax[0],
         "Generate wall particles at the 'upper end' of the x-direction in any case");
  BOOLPC(forceMaxY, m_forceMax[1],
    "Generate wall particles at the 'upper end' of the y-direction in any case");
  BOOLPC(forceMaxZ, m_forceMax[2],
    "Generate wall particles at the 'upper end' of the z-direction in any case");
  BOOLPC(forceMinX, m_forceMin[0],
    "Generate wall particles at the 'lower end' of the x-direction in any case");
  BOOLPC(forceMinY, m_forceMin[1],
    "Generate wall particles at the 'lower end' of the y-direction in any case");
  BOOLPC(forceMinZ, m_forceMin[2],
    "Generate wall particles at the 'lower end' of the z-direction in any case");
  
  BOOLPC(freeze, m_freeze,
    "Should the particles be frozen particles, i.e., all quantities fixed to "
    "the initial value and unchangeable?");
  
  FUNCTIONFIXEDPC(x, m_posX, 
                  "This sets the x-component of the position to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z'.");
  
  FUNCTIONFIXEDPC(y, m_posY, 
                  "This sets the y-component of the position to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z'.");
  
  FUNCTIONFIXEDPC(z, m_posZ, 
                  "This sets the z-component of the position to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z'.");
  
  FUNCTIONFIXEDPC(u, m_velX, 
                  "This sets the x-component of the velocity to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z', u', 'v', 'w', where 'x', 'y', 'z' represent the new modified positions (if done so).");
  
  FUNCTIONFIXEDPC(v, m_velY, 
                  "This sets the y-component of the velocity to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z', u', 'v', 'w', where 'x', 'y', 'z' represent the new modified positions (if done so).");
  
  FUNCTIONFIXEDPC(w, m_velZ, 
                  "This sets the z-component of the velocity to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z', u', 'v', 'w', where 'x', 'y', 'z' represent the new modified positions (if done so).");
  
      
  DOUBLEPC
  (distance, m_distance, 0,
  "Spacing of the lattice points.");
  DOUBLEPC
  (density, m_density, 0,
  "Number density of the particles.");
    
  m_distance = m_density = HUGE_VAL;
	
	for(int i = 0; i < SPACE_DIMS; ++i)
	{
		m_forceMax[i] = false;
		m_forceMin[i] = false;
	}

  m_freeze = "true";

  m_posX.addVariable("x");
  m_posX.addVariable("y");
  m_posX.addVariable("z");
  
  m_posY.addVariable("x");
  m_posY.addVariable("y");
  m_posY.addVariable("z");
  
  m_posZ.addVariable("x");
  m_posZ.addVariable("y");
  m_posZ.addVariable("z");
  
  m_posX.setExpression("x");
  m_posY.setExpression("y");
  m_posZ.setExpression("z");

  m_velX.addVariable("x");
  m_velX.addVariable("y");
  m_velX.addVariable("z");
  m_velX.addVariable("u");
  m_velX.addVariable("v");
  m_velX.addVariable("w");
  
  m_velY.addVariable("x");
  m_velY.addVariable("y");
  m_velY.addVariable("z");
  m_velY.addVariable("u");
  m_velY.addVariable("v");
  m_velY.addVariable("w");
  
  m_velZ.addVariable("x");
  m_velZ.addVariable("y");
  m_velZ.addVariable("z");
  m_velZ.addVariable("u");
  m_velZ.addVariable("v");
  m_velZ.addVariable("w");
  
  m_velX.setExpression("u");
  m_velY.setExpression("v");
  m_velZ.setExpression("w");
}


void ParticleCreatorWall::setup()
{
  ParticleCreator::setup();
  ParticleCreator::s_createFrozenParts = true;

  myCutoff = M_PHASE->pairCreator()->interactionCutoff();
//   myCutoff = M_SIMULATION->maxCutoff;
}


// anyhow,ParticleCreatorWall inherits Node::write(), but the following definition reminds us
// that this PC writes ITSELF in the output file and NOT a PCFile 
ostream &ParticleCreatorWall::write(ostream &s, int shift)
{
	return Node::write(s, shift);
}

