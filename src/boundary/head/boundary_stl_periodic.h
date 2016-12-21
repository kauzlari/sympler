/* 
 * File:   boundary_stl_periodic.h
 * Author: janmees
 *
 * Created on December 5, 2016, 3:19 PM
 */

#ifndef BOUNDARY_STL_PERIODIC_H
#define	BOUNDARY_STL_PERIODIC_H
/*!
 * Description of a boundary conditions to be used with stl files when
 * periodicity is desired.  
 */
#include "boundary_arbitrary.h"
#include "cell.h"


class BoundaryStlPeriodic: public BoundaryArbitrary
{
    protected: 
   /*!
   * STL file name
   */
    string m_filename;

  /*!
   * Invert the normal vector of the STL. Usually needs to be set
   * to true, I never encountered an STL where the normal vectors
   * were pointing to the inside of the geometry, but I wouldn't count 
   * on it given the simplicity/crappyness of the STL format.
   */
    bool m_invert_normals;

  
  /*!
   * Scale factor for the boundary
   */
  point_t m_scale;
  
  /*!
   * Is the STL ascii or binary? Has to be specified by the user.
   */
  string m_format;
  
  /*!
   * Periodicity of the stl geometry in x, y and or z direction.
   */

  bool_point_t m_periodic;

   
  /*!
   * Initialize the property list
   */
    void init(); 
        
    public:     
   /*!
   * Constructor
   * @param phase Pointer to the parent \a Phase object this \a Boundary belongs to
   */
    BoundaryStlPeriodic(Phase *phase);

  /*!
     * Destructor
   */
    virtual ~BoundaryStlPeriodic();
	 
     
  /*!
  * Initialize cell subdivision
  */
   virtual void setup(Simulation* sim, ManagerCell *mgr);

  /*!
     * Read and scale the boundary
   */
    virtual void setup();
};

#endif	/* BOUNDARY_STL_PERIODIC_H */

