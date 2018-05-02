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


#ifndef __PARTICLE_CREATOR_FREE_F_H
#define __PARTICLE_CREATOR_FREE_F_H

#include "pc_free.h"
#include "function_fixed.h"

//---- Forward declarations ----


//---- Classes ----

class ParticleCreatorFree;

/*!
 *\a ParticleCreatorFree using \a FunctionFixed objects for setting of initial values.
 */
class ParticleCreatorFreeF: public ParticleCreatorFree
{
  protected:
  
    /*!
   * Sets the x-component of the position to the specified algebraic expression. See init() for the known symbols.
     */
    FunctionFixed m_posX;
  
    /*!
     * Sets the y-component of the position to the specified algebraic expression. See init() for the known symbols.
     */
    FunctionFixed m_posY;
  
    /*!
     * Sets the z-component of the position to the specified algebraic expression. See init() for the known symbols.
     */
    FunctionFixed m_posZ;
  
    /*!
     * Sets the x-component of the velocity to the specified algebraic expression. See init() for the known symbols.
     */
    FunctionFixed m_velX;
  
    /*!
     * Sets the y-component of the velocity to the specified algebraic expression. See init() for the known symbols.
     */
    FunctionFixed m_velY;
  
    /*!
     * Sets the z-component of the velocity to the specified algebraic expression. See init() for the known symbols.
     */
    FunctionFixed m_velZ;
  
    /*!
    * Initialise the property list
    */
    void init();
    
    /*!
    * Probably not needed anymore
    */
    void initTransform();

    /*!
    * Applies the functions given by x...z 
    */
    void transformPos(Particle &p) {
        // so, if dependent on them, the positions are based on the OLD velocities
        p.r.x = m_posX(p.r.x, p.r.y, p.r.z, p.v.x, p.v.y, p.v.z);
        p.r.y = m_posY(p.r.x, p.r.y, p.r.z, p.v.x, p.v.y, p.v.z);
        p.r.z = m_posZ(p.r.x, p.r.y, p.r.z, p.v.x, p.v.y, p.v.z);
        
        // I think next should be done somewhere else
                
/*        for (vector<dof_info_t>::iterator i = m_internal_dofs[0].begin();
             i != m_internal_dofs[0].end(); ++i) {
               p.tag.doubleByOffset(i->offset) = i->pc->eval();
             }*/
    }
    
    /*!
     * Applies the functions given by u .. w 
     */
    void transformVel(Particle &p) {
        // so, if dependent on them, the velocities are based on the NEW positions
        p.v.x = m_velX(p.r.x, p.r.y, p.r.z, p.v.x, p.v.y, p.v.z);
        p.v.y = m_velY(p.r.x, p.r.y, p.r.z, p.v.x, p.v.y, p.v.z);
        p.v.z = m_velZ(p.r.x, p.r.y, p.r.z, p.v.x, p.v.y, p.v.z);
    }
    
    /*!
    * Builds a list of already compiled \a FunctionFixed from \a m_properties.unknown()
    */
    void fList(list<FunctionFixed>& fl);
    
    /*!
     * Compute all from \a m_properties.unknown()
    */
    void computeUnknown(list<FunctionFixed>::iterator ffi, Particle& p);
	
  public:
    /*!
    * Constructor, forbidden to be used
    */
    ParticleCreatorFreeF();
    
    /*!
    * Constructor for \a Node hierarchy
    */
    ParticleCreatorFreeF(Boundary *boundary);
    
    /*!
    * Destructor
    */
    virtual ~ParticleCreatorFreeF();

    /*!
    * Specific settings to be done
    */
    virtual void setup();
/*    virtual void flushParticles();
    virtual void flushParticles(Particle** first_p);
    virtual ostream &write(ostream &s, int shift = 0);*/
};

#endif
