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


#include "pc_free_f.h"

#include "simulation.h"
#include "manager_cell.h"

#define M_BOUNDARY ((Boundary*) m_parent)
#define M_PHASE ((Phase*) M_BOUNDARY->parent())
#define M_SIMULATION ((Simulation*) M_PHASE->parent())
#define M_MANAGER M_PHASE->manager()

ParticleCreatorFreeF::ParticleCreatorFreeF()
  : ParticleCreatorFree()
{
  throw gError("ParticleCreatorFreeF::ParticleCreatorFreeF()", "Do not use.");
}

ParticleCreatorFreeF::ParticleCreatorFreeF(Boundary *boundary)
  : ParticleCreatorFree(boundary)
{
  init();
}


ParticleCreatorFreeF::~ParticleCreatorFreeF()
{
}



//--- Methods ---

void ParticleCreatorFreeF::initTransform()
{
  // Do I still need this?
  
  /* Let's hope the size of the box doesn't change anymore after this point. */
  m_box_size = ((Boundary *) m_parent)->boundingBox().size();

}

void ParticleCreatorFreeF::init()
{
  m_properties.setClassName("ParticleCreatorFreeF");

  /* Allow unknown properties. Those ones have to be identified later.
  They are used to set the particles degrees of freedom initially. */
  m_properties.allowUnknown();
		
  m_properties.setDescription(
      "The initial conditions of user-defined symbols can be set by taking the symbol as an attribute and defining a mathematical expression for it. For the expression, the same variables are allowed as, e.g., for the attibute 'u'. For non-scalars you have to add the following to the attribute name (using two underscores):\n"
      "\"__x\", \"__y\" or \"__z\" for the respective components of a vector\n"
      "\"__xx\", \"__xy\", \"__xz\", \"__yx\", \"__yy\", \"__yz\", \"__zx\", \"__zy\" and \"__zz\" for the respective components of a tensor."
                             );
  
  FUNCTIONFIXEDPC(x, m_posX, 
                  "This sets the x-component of the postion to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z', u', 'v', 'w'. Notice that u', 'v', 'w' contain the ParticleCreator-specific initial values (e.g. zero or based on a kinetic temperature).");
  
  FUNCTIONFIXEDPC(y, m_posY, 
                  "This sets the y-component of the postion to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z', u', 'v', 'w'. Notice that u', 'v', 'w' contain the ParticleCreator-specific initial values (e.g. zero or based on a kinetic temperature).");
  
  FUNCTIONFIXEDPC(z, m_posZ, 
                  "This sets the z-component of the postion to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z', u', 'v', 'w'. Notice that u', 'v', 'w' contain the ParticleCreator-specific initial values (e.g. zero or based on a kinetic temperature).");
  
  
  FUNCTIONFIXEDPC(u, m_velX, 
                  "This sets the x-component of the velocity to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z', u', 'v', 'w'. Notice that 'x', 'y', 'z' could have been previously modified by the corresponding attributes.");
  
  FUNCTIONFIXEDPC(v, m_velY, 
                  "This sets the y-component of the velocity to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z', u', 'v', 'w'. Notice that 'x', 'y', 'z' could have been previously modified by the corresponding attributes.");
  
  FUNCTIONFIXEDPC(w, m_velZ, 
                  "This sets the z-component of the velocity to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z', u', 'v', 'w'. Notice that 'x', 'y', 'z' could have been previously modified by the corresponding attributes.");
  
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

  m_posX.addVariable("x");
  m_posX.addVariable("y");
  m_posX.addVariable("z");
  m_posX.addVariable("u");
  m_posX.addVariable("v");
  m_posX.addVariable("w");
  
  m_posY.addVariable("x");
  m_posY.addVariable("y");
  m_posY.addVariable("z");
  m_posY.addVariable("u");
  m_posY.addVariable("v");
  m_posY.addVariable("w");
  
  m_posZ.addVariable("x");
  m_posZ.addVariable("y");
  m_posZ.addVariable("z");
  m_posZ.addVariable("u");
  m_posZ.addVariable("v");
  m_posZ.addVariable("w");
  
  m_posX.setExpression("x");
  m_posY.setExpression("y");
  m_posZ.setExpression("z");
  
}

void ParticleCreatorFreeF::setup()
{
  ParticleCreatorFree::setup();
}

void ParticleCreatorFreeF::fList(list<FunctionFixed>& fl)
{
  // build a list of functions corresponding to the unknown properties from m_properties.unknown()
  list<FunctionFixed> fList;
  for (map<string, string>::const_iterator i = m_properties.unknown().begin(); i != m_properties.unknown().end(); ++i) 
  {
    FunctionFixed func;
    fl.push_back(func);
    FunctionFixed& lastF = fl.back(); 
    lastF.addVariable("x");
    lastF.addVariable("y");
    lastF.addVariable("z");
    lastF.addVariable("u");
    lastF.addVariable("v");
    lastF.addVariable("w");
    lastF.setOnlyExpr(i->second);
    lastF.compile();
  }
}

void ParticleCreatorFreeF::computeUnknown(list<FunctionFixed>::iterator ffi, Particle& p)
{
  for (map<string, string>::const_iterator i = m_properties.unknown().begin(); i != m_properties.unknown().end(); ++i) 
  {
            // MSG_DEBUG("ParticleCreatorFreeF::createParticles", "checking for unknown property " << i->first);
            // new style treatment of unknown properties
    double* value;
    DataFormat::datatype_t type;
    size_t offset;
    size_t colour = p.c;
    string name = i->first;        
    if(Particle::s_tag_format[colour].attrExists(name)) 
    {
							
      offset = Particle::s_tag_format[colour].attrByName(name).offset;
      type = Particle::s_tag_format[colour].attrByName(name).datatype;
              
      if(type == DataFormat::DOUBLE)
      {
                // value is a reference
        value = &(p.tag.doubleByOffset(offset));
      }
      else 
        throw gError("ParticleCreatorFreeF::createParticles", className() + ": For species '" + m_species + "': Scalar-type definition used for non-scalar property '" + name + "'. For vectors, you must add \"__x\", \"__y\" or \"__z\" and for tensors you must add \"__xx\", \"__xy\", \"__xz\", \"__yx\", \"__yy\", \"__yz\", \"__zx\", \"__zy\" or \"__zz\" to the property name.");
    }
    else
    { 
      if(name[name.size()-2] == '_' && name[name.size()-3] == '_')
      {
        string symbol = name;
        if(Particle::s_tag_format[colour].attrExists(symbol.erase(symbol.size()-3, symbol.size()-1)))
          type = Particle::s_tag_format[colour].attrByName(symbol).datatype;
        else 
          throw gError("ParticleCreatorFreeF::createParticles", "For species '" + m_species + ": Found a vector extension for attribute '" + name + "' but did not find symbol '" + symbol + "'.");
        if(type != DataFormat::POINT)
          if(type != DataFormat::TENSOR)
            throw gError("ParticleCreatorFreeF::createParticles", "For species '" + m_species + ": You use a vector-extension in attribute '" + name + "', but '" + symbol + "' is not a vector.");
        offset = Particle::s_tag_format[colour].attrByName(symbol).offset;
  
        switch(name[name.size()-1])
        {
          case 'x':
            value = &(p.tag.pointByOffset(offset).x); 
            break;
          case 'y':
            value = &(p.tag.pointByOffset(offset).y); 
            break;
          case 'z':
            value = &(p.tag.pointByOffset(offset).z); 
            break;
          default:
            throw gError("ParticleCreatorFreeF::createParticles", "For species '" + m_species + ": Unknown vector extension '__" + name[name.size()-1] + "' in attribute '" + name + "'. Possibilities are \"__x\", \"__y\" or \"__z\".");
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
            throw gError("ParticleCreatorFreeF::createParticles", "For species '" + m_species + ": Found a tensor extension for attribute '" + name + "' but did not find symbol '" + symbol + "'.");
          if(type != DataFormat::TENSOR)
            throw gError("ParticleCreatorFreeF::createParticles", "For species '" + m_species + ": You use a tensor-extension in attribute '" + name + "', but '" + symbol + "' is not a tensor.");
          offset = Particle::s_tag_format[colour].attrByName(symbol).offset;
                  
          bool error = false;
          switch(name[name.size()-2])
          {
            case 'x':
              switch(name[name.size()-1])
              {
                case 'x':
                  value = &(p.tag.tensorByOffset(offset).xx);
                  break;
                case 'y':
                  value = &(p.tag.tensorByOffset(offset).xy);
                  break;
                case 'z':
                  value = &(p.tag.tensorByOffset(offset).xz);
                  break;
                default:
                  error = true; 
              }
              break;
            case 'y':
              switch(name[name.size()-1])
              {
                case 'x':
                  value = &(p.tag.tensorByOffset(offset).yx);
                  break;
                case 'y':
                  value = &(p.tag.tensorByOffset(offset).yy);
                  break;
                case 'z':
                  value = &(p.tag.tensorByOffset(offset).yz);
                  break;
                default:
                  error = true; 
              }
              break;
            case 'z':
              switch(name[name.size()-1])
              {
                case 'x':
                  value = &(p.tag.tensorByOffset(offset).zx);
                  break;
                case 'y':
                  value = &(p.tag.tensorByOffset(offset).zy);
                  break;
                case 'z':
                  value = &(p.tag.tensorByOffset(offset).zz);
                  break;
                default:
                  error = true; 
              }
              break;
          }
          if(error)
            throw gError("ParticleCreatorFreeF::createParticles", "For species '" + m_species + ": Unknown tensor extension '__" + name[name.size()-2] + name[name.size()-1] + "' in attribute '" + name + "'. Possibilities are \"__xx\", \"__xy\", \"__xz\", \"__yx\", \"__yy\", \"__yz\", \"__zx\", \"__zy\" or \"__zz\".");

        } // end if(tensor-case)
        else
        {
          MSG_DEBUG("ParticleCreatorFreeF::createParticles", className() << ": NOT FOUND " << i->first);
          throw gError("ParticleCreatorFreeF::createParticles", 
                       className() + ": No internal degree of freedom named '" + i->first +
                           "' for species '" + m_species + "'. Possibilities are: " + Particle::s_tag_format[colour].toStr/*ing*/());
        }
      } // end else of if(vector-case)
    } // end else of if(attrExists)
            
    *value = (*ffi)(p.r.x, p.r.y, p.r.z, p.v.x, p.v.y, p.v.z);
    ++ffi;
  }
 
}
