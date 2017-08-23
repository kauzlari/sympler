extern "C" {
  #include <freesteam/steam_ps.h>
  #include <freesteam/steam_pT.h>
  #include <freesteam/region4.h>
}
#include "density_calculation.h"
#include "simulation.h"
#include "manager_cell.h"
#include <iostream>
#include <string>

const SymbolRegister<DensityCalculation> DensityCalculation("DensityCalculation");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

DensityCalculation::DensityCalculation
  (size_t colour, size_t offset, string symbolName)
  : ParticleCache(colour, offset, symbolName)/*m_colour(colour), m_offset(offset), *//*m_wf(wf)*/ {
  m_stage = 0;
  m_datatype = DataFormat::DOUBLE;
}

DensityCalculation::DensityCalculation
    (/*Node*/Simulation* parent)
  : ParticleCache(parent) {
  m_stage = 0;
  m_datatype = DataFormat::DOUBLE;
  init();
}


DensityCalculation::~DensityCalculation() {
  for (int i = 0; i < m_arraysize_temperature ; ++i) {
    delete [] m_array_rho[i];
  }
  delete [] m_array_rho;
}

void DensityCalculation::setupLUT(double Tmin, double pmin, double Tmax, double pmax, int m_arraysize_pressure, int m_arraysize_temperature) {
  // Function intern auxiliary vaiables
  double pr = 0;
  double t = 0;
  double density;
  // Step sizes depending on the size of the Array and the ranges of the input values.
  m_calcstepP= (pmax-pmin)/(m_arraysize_pressure-1);
  m_calcstepT= (Tmax-Tmin)/(m_arraysize_temperature-1);
  // Initialization of the 2D Array, depending on the given temperature and pressure sizes.
  m_array_rho = new double*[m_arraysize_pressure+1];
  for (int i = 0; i <= m_arraysize_pressure; i++) {
    m_array_rho[i] = new double [m_arraysize_temperature+1];
  }
  // Loop calculates the densities in an adjusted temperature and pressure range and
  // stores them into the LUT.
  for(int j = 0; j < m_arraysize_pressure; j++) {
    for (int i = 0; i < m_arraysize_temperature; i++) {  
      SteamState S = freesteam_set_pT(pmin + pr,Tmin + t);
      density = freesteam_rho(S);
      m_array_rho[j][i] = density;
      t += m_calcstepT;
    }
    t = 0;
    pr += m_calcstepP;
  }
}

double DensityCalculation::calculateDensity(double inputT, double inputP, double Tmin, double pmin) {
  // Calculation of the surrounded sampling points    .
  int x_pressure_array_0 = (floor((inputP-pmin)/m_calcstepP));
  int y_temperature_array_0 =(floor((inputT-Tmin)/m_calcstepT));
  int x_pressure_array_1 = (floor((inputP-pmin)/m_calcstepP)+1);
  int y_temperature_array_1 =(floor((inputT-Tmin)/m_calcstepT)+1);
  // Normalization of the input values.
  double press_nominated =  ((inputP -(x_pressure_array_0*m_calcstepP+pmin)) / m_calcstepP);
  double temp_nominated =   ((inputT -(y_temperature_array_0*m_calcstepT+Tmin)) / m_calcstepT);
  // Bilinear interpolation of the normalized values.
  double x_0_y_0 = m_array_rho[x_pressure_array_0][y_temperature_array_0]*(1-press_nominated)*(1-temp_nominated);
  double x_1_y_0 = m_array_rho[x_pressure_array_1][y_temperature_array_0]*press_nominated*(1-temp_nominated);
  double x_0_y_1 = m_array_rho[x_pressure_array_0][y_temperature_array_1]*(1-press_nominated)*temp_nominated;
  double x_1_y_1 = m_array_rho[x_pressure_array_1][y_temperature_array_1]*temp_nominated*press_nominated;
  // Interpolated density value.
  m_density_interpolation = x_0_y_0 + x_1_y_0 + x_0_y_1 + x_1_y_1;
  return m_density_interpolation;
}



void DensityCalculation::init() {
  m_properties.setClassName("DensityCalculation");
  m_properties.setName("DensityCalculation");
  m_properties.setDescription("Contribution of the particle itself to its local density.");  
  STRINGPC
      (symbol, m_symbolName,
       "Name for the local density.");
  STRINGPC
      (pressure, m_pressureName,
       "Name for the local pressure.");
  STRINGPC
      (temperature, m_temperatureName,
       "Name for the local temperature.");
  DOUBLEPC
      (temperatureMax, m_Tmax, 0,
       "The local maximum temperature.");
  DOUBLEPC
      (temperatureMin, m_Tmin, 0,
       "The local maximum temperature.");
  DOUBLEPC
      (pressureMax, m_pmax, 0,
       "The local maximum pressure.");
  DOUBLEPC
      (pressureMin, m_pmin, 0,
       "The local minimum pressure.");
  INTPC
      (arraysize_pressure, m_arraysize_pressure, 0,
       "The size of the array for pressure values.");
  INTPC
      (arraysize_temperature, m_arraysize_temperature, 0,
       "The size of the array for temperature values.");

  m_temperatureName = "undefined";
  m_pressureName = "undefined";
  m_Tmin = 0;
  m_pmin = 0;
  m_Tmax = 0;
  m_pmax = 0;
  m_arraysize_pressure = 0;
  m_arraysize_temperature = 0;
}

void DensityCalculation::setup() {
  
  // Checks if all necessary input values are defined.
  if(m_species == "undefined")
    throw gError("DensityCalculation::setup", "Attribute 'species' has value \"undefined\""); 
  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("DensityCalculation::setup", "Attribute 'stage' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_temperatureName == "undefined")
    throw gError("DensityCalculation::setup", "Attribute 'temperature' has value \"undefined\""); 
  if(m_pressureName == "undefined")
    throw gError("DensityCalculation::setup", "Attribute 'pressure' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_Tmin == 0)
    throw gError("DensityCalculation::setup", "Attribute 'temperatureMin' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_Tmax == 0)
    throw gError("DensityCalculation::setup", "Attribute 'temperatureMax' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_pmin == 0)
    throw gError("DensityCalculation::setup", "Attribute 'pressureMin' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_pmax == 0)
    throw gError("DensityCalculation::setup", "Attribute 'pressureMax' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_pmax == 0)
    throw gError("DensityCalculation::setup", "Attribute 'arraysize_temperature' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_pmax == 0)
    throw gError("DensityCalculation::setup", "Attribute 'arraysize_pressure' has none of the allowed values \"0\", \"1\", \"2\".");

  
  pair<size_t, size_t> tempPair;
  
  // should we create a Cache for the other colours, too?
  if(m_species == "ALL")
  {
    // This calculator does not care if the symbol is already existing
    for (m_colour = 0; m_colour < M_MANAGER->nColours(); ++m_colour)
    {
      // at least, if the attribute already exists, we preserve the persistency and can hope that another module has set it correctly
      if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
      {
        if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_symbolName).datatype)
          throw gError("DensityCalculation::setup", "Symbol " + m_symbolName + " already exists as a non-scalar.");
        else m_offset = Particle::s_tag_format[m_colour]./*indexOf*/offsetByName(m_symbolName);

      }
      else
      // if it does not yet exist we set the persistency to false
        m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false/*persistency*/, m_symbolName).offset;
      if(Particle::s_tag_format[m_colour].attrExists(m_temperatureName) && Particle::s_tag_format[m_colour].attrExists(m_pressureName)) {
        if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_temperatureName).datatype)
          throw gError("DensityCalculation::setup", "Temperature " + m_temperatureName + " already exists as a non-scalar.");
        else m_temperatureOffset = Particle::s_tag_format[m_colour]./*indexOf*/offsetByName(m_temperatureName);
        if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_pressureName).datatype)
          throw gError("DensityCalculation::setup", "Pressure " + m_pressureName + " already exists as a non-scalar.");
        else m_pressureOffset = Particle::s_tag_format[m_colour]./*indexOf*/offsetByName(m_pressureName);
      } else 
          throw gError("DensityCalculation::setup", "Pressure " + m_pressureName + " Temperature " + m_temperatureName + " do not exist.");
     
      // is it the last cache to be created?
      if(m_colour == M_MANAGER->nColours()-1)
      {
        if(m_phaseUser == 0)
          Particle::registerCache_0(this);
        else if(m_phaseUser == 1) 
          Particle::registerCache(this);
        else // so it is 2
        {
          ParticleCache* pc = new DensityCalculation(*this);
          assert(pc->mySymbolName() == m_symbolName);
          assert(((DensityCalculation*) pc)->stage() == m_stage);
          assert(((DensityCalculation*) pc)->m_colour == m_colour);
          assert(((DensityCalculation*) pc)->m_offset == m_offset);
	  assert(((DensityCalculation*) pc)->m_temperatureOffset == m_temperatureOffset);
	  assert(((DensityCalculation*) pc)->m_pressureOffset == m_pressureOffset);
          Particle::registerCache(pc);
          Particle::registerCache_0(this);
        }
      }
      // No? Then make a copy
      else 
      {
        ParticleCache* pc = new DensityCalculation(*this);
        assert(pc->mySymbolName() == m_symbolName);
        assert(((DensityCalculation*) pc)->stage() == m_stage);
        assert(((DensityCalculation*) pc)->m_colour == m_colour);
        assert(((DensityCalculation*) pc)->m_offset == m_offset);
	assert(((DensityCalculation*) pc)->m_temperatureOffset == m_temperatureOffset);
	assert(((DensityCalculation*) pc)->m_pressureOffset == m_pressureOffset);

        if(m_phaseUser == 0)
          Particle::registerCache_0(pc);
        else if(m_phaseUser == 1)
          Particle::registerCache(pc);
        else // so it is 2
        {
          Particle::registerCache(pc);
          
          pc = new DensityCalculation(*this);
          assert(pc->mySymbolName() == m_symbolName);
          assert(((DensityCalculation*) pc)->stage() == m_stage);
          assert(((DensityCalculation*) pc)->m_colour == m_colour);
          assert(((DensityCalculation*) pc)->m_offset == m_offset);
	  assert(((DensityCalculation*) pc)->m_temperatureOffset == m_temperatureOffset);
	  assert(((DensityCalculation*) pc)->m_pressureOffset == m_pressureOffset);
          
          Particle::registerCache_0(pc);
        }
      
      }
    }
  }
  else /*it is a Symbol limited to one colour*/
  {
    m_colour = M_MANAGER->getColour/*AndAdd*/(m_species);
    
    // at least, if the attribute alrteady exists, we preserve the persistency and can hope that an other module has set it correctly
    if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
    {
      if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_symbolName).datatype)
        throw gError("ParticleCacheDensitySelfContribution::setup", "Symbol " + m_symbolName + " already exists as a non-scalar.");
      else m_offset = Particle::s_tag_format[m_colour]./*indexOf*/offsetByName(m_symbolName);
    }
    else
      // if it does not yet exist we set the persistency to false
      m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false/*persistency*/, m_symbolName).offset;

    if(Particle::s_tag_format[m_colour].attrExists(m_temperatureName) && Particle::s_tag_format[m_colour].attrExists(m_pressureName)) {
      if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_temperatureName).datatype)
        throw gError("PressureCalculation::setup", "Temperature " + m_temperatureName + " already exists as a non-scalar.");
      else m_temperatureOffset = Particle::s_tag_format[m_colour]./*indexOf*/offsetByName(m_temperatureName);
        if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_pressureName).datatype)
          throw gError("PressureCalculation::setup", "Pressure " + m_pressureName + " already exists as a non-scalar.");
        else m_pressureOffset = Particle::s_tag_format[m_colour]./*indexOf*/offsetByName(m_pressureName);
    } else 
          throw gError("PressureCalculation::setup", "Pressure " + m_pressureName + " Temperature " + m_temperatureName + " do not exist.");
      if(Particle::s_tag_format[m_colour].attrExists(m_temperatureName) && Particle::s_tag_format[m_colour].attrExists(m_pressureName)) {
        if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_temperatureName).datatype)
          throw gError("PressureCalculation::setup", "Temperature " + m_temperatureName + " already exists as a non-scalar.");
        else m_temperatureOffset = Particle::s_tag_format[m_colour]./*indexOf*/offsetByName(m_temperatureName);
        if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_pressureName).datatype)
          throw gError("PressureCalculation::setup", "Pressure " + m_pressureName + " already exists as a non-scalar.");
        else m_pressureOffset = Particle::s_tag_format[m_colour]./*indexOf*/offsetByName(m_pressureName);
      } else 
          throw gError("PressureCalculation::setup", "Pressure " + m_pressureName + " Temperature " + m_temperatureName + " do not exist.");


    if(m_phaseUser == 0)
      Particle::registerCache_0(this);
    else if(m_phaseUser == 1)
      Particle::registerCache(this);
    else
    {
      ParticleCache* pc = new DensityCalculation(*this);
      assert(pc->mySymbolName() == m_symbolName);
      assert(((DensityCalculation*) pc)->stage() == m_stage);
      assert(((DensityCalculation*) pc)->m_colour == m_colour);
      assert(((DensityCalculation*) pc)->m_offset == m_offset);
      assert(((DensityCalculation*) pc)->m_temperatureOffset == m_temperatureOffset);
      assert(((DensityCalculation*) pc)->m_pressureOffset == m_pressureOffset);

      Particle::registerCache(pc);
      Particle::registerCache_0(this);
    } 
  }
  /*  Particle* p;
    m_pmax = p->tag.doubleByOffset(m_pressureMaxOffset);
    m_Tmax = p->tag.doubleByOffset(m_temperatureMaxOffset);
    m_pmin = p->tag.doubleByOffset(m_pressureMinOffset);
    m_Tmin = p->tag.doubleByOffset(m_temperatureMinOffset);
    printf("%f", m_pmax); */
   setupLUT( m_Tmin,  m_pmin, m_Tmax, m_pmax, m_arraysize_pressure, m_arraysize_temperature);   
}


void DensityCalculation::registerWithParticle()
{
}
