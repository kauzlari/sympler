
#include "boundary_stl_periodic.h"

#include "stl_ascii.h"
#include "stl_binary.h"
#include "simulation.h"
#include "manager_cell.h"
#include "wall_triangle.h"


/* Register this Boundary with the factory. */
const Boundary_Register<BoundaryStlPeriodic> boundary_stl_periodic("BoundaryStlPeriodic");

//---- Constructors/Destructor ----

BoundaryStlPeriodic::BoundaryStlPeriodic(Phase *phase): BoundaryArbitrary(phase)
{
  init();
}


BoundaryStlPeriodic::~BoundaryStlPeriodic()
{
}
#define STL_ASCII "ascii"
#define STL_BINARY "binary"

//---- Methods ----

void BoundaryStlPeriodic::init()
{
  MSG_DEBUG("BoundaryStlPeriodic", "Initializing boundary.");  
  m_properties.setClassName("BoundaryStlPeriodic");
  m_properties.setName("BoundaryStlPeriodic");
  m_properties.setDescription("A Boundary with periodic stl geometry. ");
  
  for (int i = 0; i < 3; i++) {
    /* Register all properties */
    m_properties.addProperty
      ("scale" + string(1, 'X'+i), PropertyList::DOUBLE, 
       &m_scale[i],
       new PLCDoubleGreater(0),
       string(1, 'x'+i) + "-scaling factor for the geometry.");
    
    m_properties.addProperty
      ("periodic" + string(1, 'X'+i), PropertyList::BOOLEAN, &m_periodic[i], NULL,
       "Simulation box is periodic in " + string(1, 'x'+i) + "-direction.");
    
    //Defaults
    m_periodic[i] = false;
    m_scale[i] = 1;
  }
  
  STRINGPCINF
    (name, m_filename,
     "Name of the .STL geometry file to be loaded. Geometry must have inflow and outflow defined by extremal points. Inflow and outflow must be at same height and have same shape.");

  BOOLPC
    (invertNormals, m_invert_normals,
     "Invert the direction of the surface normal vectors.");

  m_properties.addProperty
    ("format", PropertyList::STRING,
     &m_format,
     new PLCStringList2(STL_ASCII, STL_BINARY),
     "Type of the STL: 'ascii' or 'binary'.");
  
  
  
  //Defaults
  
  m_filename = "default.stl";
  m_invert_normals = true;
  m_format = "ascii";  
}

    /*!
     * Read and scale the boundary
    */
void BoundaryStlPeriodic::setup()
{
  
    BoundaryArbitrary::setup();

    m_container = new WallContainer(m_reflectors[0]); 
    
    if (m_format == STL_ASCII) {
    STLAsciiParser loader(m_container);

    loader.setInvertNormals(m_invert_normals);

    m_container = loader.read(m_filename);
    } 
    else {
    STLBinaryParser loader(m_container);

    loader.setInvertNormals(m_invert_normals);

    m_container = loader.read(m_filename);
    }
 
/*  bool_point_t periodicFrontandBack = { { { m_periodic[0], m_periodic[1], m_periodic[2] } } };
  m_container->setPeriodicityFront(periodicFrontandBack);
  m_container->setPeriodicityBack(periodicFrontandBack);

 */
    m_container->setPeriodicityFront(m_periodic);
    m_container->setPeriodicityBack(m_periodic);  
    m_container->stretchBy(m_scale);
    m_container->updateBoundingBox();
//MSG_DEBUG("BoundaryStlPeriodic::read", "box: " << m_container->boundingBox().corner1 << " - " << m_container->boundingBox().corner2);

  


/*ALGORITHM
 * 1.- Find inflow outflow surfaces. 
 *Loop over all walls and find minimum and maximum corners along the axis defined as periodic.
 *If a wall is on the inlet or outlet all corners will share the same distance in this direction. 
 * for example c1(x1,y1,z1) c2(x1,y2,z2) c3(x1,y3,z3)    : All corners same x1!
 * 2.-
 * If this condition is met the wall will be deleted and with it, its associated reflectors.  
 */
    for (int j=0 ; j<3; j++){
   
        if(m_periodic[j] == true){
      
            double min, max, temp;
            list<Wall*>::iterator d; 
            int countermax, countermin; 
               
            for (list<Wall*>::iterator i = m_container->walls().begin(); i != m_container->walls().end(); i++) { 
                for(int p =0; p<3;p++){  
                    temp= (*i)->findcorner(p)[j];
                    
                    if(temp > max)
                        max= temp;                        
                    
                    if(temp< min)
                        min= temp;                       
                     
                }   
            }
     
            for (list<Wall*>::iterator i = m_container->walls().begin(); i != m_container->walls().end(); ) { 
                countermin =0; countermax=0; 
                
                for(int p =0; p<3;p++){ 
                    temp= (*i)->findcorner(p)[j];
                                      
                    if(temp==min) 
                        ++countermin;                
                
                    if(temp==max)
                        ++countermax;    
                }

                if (countermin == 3 || countermax ==3 ){  
                    d=i; 
                    i++;
                    m_container->deleteWall(d);
              
                }
                else i++;
            }
       
  
        }  
    } 
     MSG_DEBUG("Boundary Stl Periodic", "Setup completed.");
}





//---- Cell subdivision ----


#define OUTLET 0
#define NEIGHBOR 1

void BoundaryStlPeriodic::setup(Simulation* sim, ManagerCell *mgr)
{
  BoundaryArbitrary::setup(sim, mgr);
   MSG_DEBUG
    ("BoundaryStlPeriodic::setup",
     "corner1 = " << m_container->boundingBox().corner1 << ", " <<
     "corner2 = " << m_container->boundingBox().corner2);

  MSG_DEBUG("BoundaryStlPeriodic::setup", "Initializing cells.");
  
  region_t *r;
  point_t c1, c2 ; 
  
  c1 = m_container->boundingBox().corner1;
  c2 = m_container->boundingBox().corner2; 
  
  r = mgr->cellSubdivide(sim->maxCutoff, c1, c2, m_periodic, 1, NULL); // 
    
  mgr->cellSubdivisionFinished();
  
  MSG_DEBUG("BoundaryStlPeriodic::setup", "Cell subdivision Finished.");
 
  mgr->assignContainer(m_container);
  
  m_container->toVTK(m_geometry_filename);
}