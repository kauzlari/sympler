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



#include "function_pair_fixed.h"

#include "phase.h"
#include "fp_scalar.h"
#include "simulation.h"
#include "colour_pair.h"
#include "manager_cell.h"
#include "val_calculator_rho.h"
#include "val_calculator_kernel.h"
#include "weighting_function.h"

#include "pca_density_0oc.h"
#include "pca_volume_self_contribution.h"

// #include "val_calculator_shear.h"

#include "fp_vector.h"
#include "fp_tensor.h"


#define M_MANAGER m_cp->manager()
#define M_PHASE M_MANAGER->phase()
#define M_SIMULATION ((Simulation*) M_PHASE->parent())


/* FunctionPair */

FunctionPairFixed::FunctionPairFixed()
{
//   MSG_DEBUG("FunctionPair::FunctionPair()", "CALLED");
}


FunctionPairFixed::~FunctionPairFixed()
{
}



class FunctionPairFixedCallback: public UnknownSymbolCallback
{
  protected:
    ColourPair *m_cp;
//   WeightingFunction *m_wf;

  public:
    FunctionPairFixedCallback(ColourPair *cp/*, WeightingFunction *wf*/): m_cp(cp) /*, m_wf(wf)*/ 
    { 
/*    MSG_DEBUG("FunctionPairCallback::Constructor", "CALLED");
      MSG_DEBUG("FunctionPairCallback::Constructor", "m_wf=" << m_wf->name());
      MSG_DEBUG("FunctionPairCallback::Constructor", "SURVIVED");  */
    }

    virtual TypedValue *operator()(const string &e) {
      /* Now let's see... */

//     throw gError("FunctionPairCallback::operator()", "Unknown expression '" + e + "'.");
    
      return NULL;
    
#if 0
    
//     cout << "Request for symbol: '" << e << "'" << endl;

    // is it a strain rate with interpolated local density?
    if (e[0] == '{' && e[1] == 'D' && (e[2] == 'i' || e[2] == 'j') && e[3] == '}') {
  /* Request for a strain rate. */
  string tmpString[SPACE_DIMS][SPACE_DIMS];
      
  if (e.size() == 4 && m_wf != NULL) {
// MSG_DEBUG("FunctionPairCallback::operator()", "found special, size=2: " << e);
    pair<size_t, size_t> shear_offset;
      
        
        // we consider this to be a non-ColourPair-specific strain rate
        // first we register the value in the CPs != m_cp
    FOR_EACH_COLOUR_PAIR
        (
        M_MANAGER,
          // if this is m_cp then do nothing (will be done afterwards)
    if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
      cp->registerCalc(shear_offset, new ValCalculatorShear(m_wf, "none"), true);
        );
        // now consider m_cp and create and return the new variables
    m_cp->registerCalc(shear_offset, new ValCalculatorShear(m_wf, "none"), true);
        
    MSG_DEBUG("FunctionPairCallback::operator()", "registering D: offset=" << shear_offset.first << ", " << shear_offset.second);
        
    if (e[2] == 'i') {
      FunctionArbitrary::tensor2CExpression(tmpString, "first_tag", shear_offset.first);
      return new FPTensorVariable
          (e,
           NULL,
           tmpString);
    } else if (e[2] == 'j') {
      FunctionArbitrary::tensor2CExpression(tmpString, "second_tag", shear_offset.second);
      return new FPTensorVariable
          (e,
           NULL,
           tmpString);
    }
  } 
  else if (e.size() > 4 && e[4] == '[' && e[e.size()-1] == ']') {
    string wf_name = string(e, 5, e.size()-6);
    WeightingFunction *wf;
    pair<size_t, size_t> shear_offset;
      
    cout << "wf = " << wf_name << endl;
      
    wf = M_SIMULATION->findWeightingFunction(wf_name);
      
    if (!wf)
      throw gError
          ("FunctionPairCallback::operator()",
           "Unknown weighting function: '" + wf_name + "'.");
        
        // we consider this to be a non-ColourPair-specific strain rate
        // first we register the value in the CPs != m_cp
    FOR_EACH_COLOUR_PAIR
        (
        M_MANAGER,
          // if this is m_cp then do nothing (will be done afterwards)
    if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
      cp->registerCalc(shear_offset, new ValCalculatorShear(wf, "none"), true);
        );
        // now consider m_cp and create and return the new variables
    m_cp->registerCalc(shear_offset, new ValCalculatorShear(wf, "none"), true);
      
    if (e[2] == 'i') {
      FunctionArbitrary::tensor2CExpression(tmpString, "first_tag", shear_offset.first);
      return new FPTensorVariable
          (e,
           NULL,
           tmpString);
    } else if (e[2] == 'j') {
      FunctionArbitrary::tensor2CExpression(tmpString, "second_tag", shear_offset.first);
      return new FPTensorVariable
          (e,
           NULL,
           tmpString);
    }
  }
    }

    // is it a strain rate with local density taken from a scalar?
    // This one does not support the definition of a specific weighting function
    if (e[0] == '{' && e[1] == 'D' && e[2] == '_' && (e[e.size()-2] == 'i' || e[e.size()-2] == 'j') && e[e.size()-1] == '}') {
      
      
      // Cut the scalar name out of {Di_scalar}
      string scalar_name = string(e, 3, e.size()-5);
      
      /* Request for a strain rate. */
      string tmpString[SPACE_DIMS][SPACE_DIMS];
      
      pair<size_t, size_t> shear_offset;
    
      
      // we consider this to be a non-ColourPair-specific strain rate
      // first we register the value in the CPs != m_cp
      FOR_EACH_COLOUR_PAIR
          (
          M_MANAGER,
        // if this is m_cp then do nothing (will be done afterwards)
      if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
        cp->registerCalc(shear_offset, new ValCalculatorShear(m_wf, scalar_name), true);
          );
      // now consider m_cp and create and return the new variables
      m_cp->registerCalc(shear_offset, new ValCalculatorShear(m_wf, scalar_name), true);
      
      MSG_DEBUG("FunctionPairCallback::operator()", "registering D: offset=" << shear_offset.first << ", " << shear_offset.second);
      
      if (e[2] == 'i') {
        FunctionArbitrary::tensor2CExpression(tmpString, "first_tag", shear_offset.first);
        return new FPTensorVariable
            (e,
             NULL,
             tmpString);
      } else if (e[2] == 'j') {
        FunctionArbitrary::tensor2CExpression(tmpString, "second_tag", shear_offset.second);
        return new FPTensorVariable
            (e,
             NULL,
             tmpString);
      }
    }

    // is it an uncorrected density?
    // we consider this to be a non-ColourPair-specific local density
    if (e[0] == 'n' && (e[1] == 'i' || e[1] == 'j')) {

/*MSG_DEBUG("FunctionPairCallback::operator()", "found special: " << e);
      if(e.size() == 2) MSG_DEBUG("FunctionPairCallback::operator()", "2");
      if(m_wf == NULL) MSG_DEBUG("FunctionPairCallback::operator()", "NULL");*/
      
      /* Request for a density. */
      if (e.size() == 2 && m_wf != NULL) {
// MSG_DEBUG("FunctionPairCallback::operator()", "found special, size=2: " << e);
        pair<size_t, size_t> density_offset;
        
        // we consider this to be a non-ColourPair-specific local density
        // first we register the value in the CPs != m_cp
        FOR_EACH_COLOUR_PAIR
            (
            M_MANAGER,
          // if this is m_cp then do nothing (will be done afterwards)
        if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
          cp->registerCalc(density_offset, new ValCalculatorRho(m_wf, ""), true);
            );
        m_cp->registerCalc(density_offset, new ValCalculatorRho(m_wf, ""), true);
      
        if (e[1] == 'i') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("first_tag", density_offset.first));
        } else if (e[1] == 'j') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("second_tag", density_offset.second));
        }
      } 
      else if (e.size() > 2 && e[2] == '[' && e[e.size()-1] == ']') {
        string wf_name = string(e, 3, e.size()-4);
        WeightingFunction *wf;
        pair<size_t, size_t> density_offset;
      
        cout << "wf = " << wf_name << endl;
      
        wf = M_SIMULATION->findWeightingFunction(wf_name);
      
        if (!wf)
          throw gError
              ("FunctionPairCallback::operator()",
               "Unknown weighting function: '" + wf_name + "'.");
      
        // we consider this to be a non-ColourPair-specific local density
        // first we register the value in the CPs != m_cp
        FOR_EACH_COLOUR_PAIR
            (
            M_MANAGER,
          // if this is m_cp then do nothing (will be done afterwards)
        if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
          cp->registerCalc(density_offset, new ValCalculatorRho(wf, ""), true);
            );
        m_cp->registerCalc(density_offset, new ValCalculatorRho(wf, ""), true);
      
        if (e[1] == 'i') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("first_tag", density_offset.first));
        } else if (e[1] == 'j') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("second_tag", density_offset.second));
        }
      }
    }

    // is it a density with 0th order correction, which is used for all 
    // ColourPairs?
    if (e[0] == 'n' && e[1] == '0' && e[2] == 'c' && (e[3] == 'i' || e[3] == 'j')) {

/*MSG_DEBUG("FunctionPairCallback::operator()", "found special: " << e);
      if(e.size() == 2) MSG_DEBUG("FunctionPairCallback::operator()", "2");
      if(m_wf == NULL) MSG_DEBUG("FunctionPairCallback::operator()", "NULL");*/
      
      /* Request for a density. */
      if (e.size() == 4 && m_wf != NULL) {
// MSG_DEBUG("FunctionPairCallback::operator()", "found special, size=2: " << e);
        pair<size_t, size_t> density_offset;
      
        // first we register the value in the CPs != m_cp
        FOR_EACH_COLOUR_PAIR
            (
            M_MANAGER,
          // if this is m_cp then do nothing (will be done afterwards)
        if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
          cp->registerCalc(density_offset, new ValCalculatorRho(m_wf, "0c"), true);
            );
        m_cp->registerCalc(density_offset, new ValCalculatorRho(m_wf, "0c"), true);
        // next will also generate a ParticleCalculatorVolume
        Particle::registerCache(new ParticleCacheDensity0Oc(m_cp, m_cp->firstColour(), density_offset.first, m_wf, "0c", true));
        if(m_cp->firstColour() != m_cp->secondColour())
          Particle::registerCache(new ParticleCacheDensity0Oc(m_cp, m_cp->secondColour(), density_offset.second, m_wf, "0c", true));
        
        if (e[3] == 'i') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("first_tag", density_offset.first));
        } else if (e[3] == 'j') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("second_tag", density_offset.second));
        }
      } 
      else if (e.size() > 4 && e[4] == '[' && e[e.size()-1] == ']') {
        string wf_name = string(e, 5, e.size()-6);
        WeightingFunction *wf;
        pair<size_t, size_t> density_offset;
      
        MSG_DEBUG("FunctionCallback::operator()", "searching wf '" << wf_name << "'.");
      
        wf = M_SIMULATION->findWeightingFunction(wf_name);
      
        if (!wf)
          throw gError
              ("FunctionPairCallback::operator()",
               "Unknown weighting function: '" + wf_name + "'.");
      
        // first we register the value in the CPs != m_cp
        FOR_EACH_COLOUR_PAIR
            (
            M_MANAGER,
          // if this is m_cp then do nothing (will be done afterwards)
        if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
          cp->registerCalc(density_offset, new ValCalculatorRho(wf, "0c"), true);
            );
        m_cp->registerCalc(density_offset, new ValCalculatorRho(wf, "0c"), true);
        Particle::registerCache(new ParticleCacheDensity0Oc(m_cp, m_cp->firstColour(), density_offset.first, wf, "0c", true));
        if(m_cp->firstColour() != m_cp->secondColour())
          Particle::registerCache(new ParticleCacheDensity0Oc(m_cp, m_cp->secondColour(), density_offset.second, wf, "0c", true));
      
        if (e[3] == 'i') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("first_tag", density_offset.first));
        } else if (e[3] == 'j') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("second_tag", density_offset.second));
        }
      }
    }

    // is it a standard density, which is used for only this specific ColourPair?
    if (e[0] == 'n' && e[1] == 'C' && e[2] == 'P' && (e[3] == 'i' || e[3] == 'j')) {

/*MSG_DEBUG("FunctionPairCallback::operator()", "found special: " << e);
      if(e.size() == 2) MSG_DEBUG("FunctionPairCallback::operator()", "2");
      if(m_wf == NULL) MSG_DEBUG("FunctionPairCallback::operator()", "NULL");*/
      
      /* Request for a density. */
      if (e.size() == 4 && m_wf != NULL) {
// MSG_DEBUG("FunctionPairCallback::operator()", "found special, size=2: " << e);
        pair<size_t, size_t> density_offset;
      
        m_cp->registerCalc(density_offset, new ValCalculatorRho(m_wf, "CP"), false);
        
        if (e[3] == 'i') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("first_tag", density_offset.first));
        } else if (e[3] == 'j') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("second_tag", density_offset.second));
        }
      } 
      else if (e.size() > 4 && e[4] == '[' && e[e.size()-1] == ']') {
        string wf_name = string(e, 5, e.size()-6);
        WeightingFunction *wf;
        pair<size_t, size_t> density_offset;
      
        MSG_DEBUG("FunctionCallback::operator()", "searching wf '" << wf_name << "'.");
      
        wf = M_SIMULATION->findWeightingFunction(wf_name);
      
        if (!wf)
          throw gError
              ("FunctionPairCallback::operator()",
               "Unknown weighting function: '" + wf_name + "'.");
      
        m_cp->registerCalc(density_offset, new ValCalculatorRho(wf, "CP"), false);
        Particle::registerCache(new ParticleCacheDensity0Oc(m_cp, m_cp->firstColour(), density_offset.first, wf, "CP", false));
        if(m_cp->firstColour() != m_cp->secondColour())
          Particle::registerCache(new ParticleCacheDensity0Oc(m_cp, m_cp->secondColour(), density_offset.second, wf, "CP", false));
      
        if (e[3] == 'i') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("first_tag", density_offset.first));
        } else if (e[3] == 'j') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("second_tag", density_offset.second));
        }
      }
    }

    // is it a density with 0th order correction, which is used for only this specific ColourPair?
    if (e[0] == 'n' && e[1] == '0' && e[2] == 'c' && e[3] == 'C' && e[4] == 'P' && (e[5] == 'i' || e[5] == 'j')) {

/*MSG_DEBUG("FunctionPairCallback::operator()", "found special: " << e);
      if(e.size() == 2) MSG_DEBUG("FunctionPairCallback::operator()", "2");
      if(m_wf == NULL) MSG_DEBUG("FunctionPairCallback::operator()", "NULL");*/
      
      /* Request for a density. */
      if (e.size() == 6 && m_wf != NULL) {
// MSG_DEBUG("FunctionPairCallback::operator()", "found special, size=2: " << e);
        pair<size_t, size_t> density_offset;
      
        m_cp->registerCalc(density_offset, new ValCalculatorRho(m_wf, "0cCP"), false);
        // next will also generate a ParticleCalculatorVolume
        Particle::registerCache(new ParticleCacheDensity0Oc(m_cp, m_cp->firstColour(), density_offset.first, m_wf, "0cCP", false));
        if(m_cp->firstColour() != m_cp->secondColour())
          Particle::registerCache(new ParticleCacheDensity0Oc(m_cp, m_cp->secondColour(), density_offset.second, m_wf, "0cCP", false));
        
        if (e[5] == 'i') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("first_tag", density_offset.first));
        } else if (e[5] == 'j') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("second_tag", density_offset.second));
        }
      } 
      else if (e.size() > 6 && e[6] == '[' && e[e.size()-1] == ']') {
        string wf_name = string(e, 7, e.size()-8);
        WeightingFunction *wf;
        pair<size_t, size_t> density_offset;
      
        MSG_DEBUG("FunctionCallback::operator()", "searching wf '" << wf_name << "'.");
      
        wf = M_SIMULATION->findWeightingFunction(wf_name);
      
        if (!wf)
          throw gError
              ("FunctionPairCallback::operator()",
               "Unknown weighting function: '" + wf_name + "'.");
      
        m_cp->registerCalc(density_offset, new ValCalculatorRho(wf, "0cCP"), false);
        Particle::registerCache(new ParticleCacheDensity0Oc(m_cp, m_cp->firstColour(), density_offset.first, wf, "0cCP", false));
        if(m_cp->firstColour() != m_cp->secondColour())
          Particle::registerCache(new ParticleCacheDensity0Oc(m_cp, m_cp->secondColour(), density_offset.second, wf, "0cCP", false));
      
        if (e[5] == 'i') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("first_tag", density_offset.first));
        } else if (e[5] == 'j') {
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("second_tag", density_offset.second));
        }
      }
    }

    // test, whether PCAVolume will and PCADensity0OC are placed into the 
    // correct stages
    // first colour
    for(size_t i = 0; i < Particle::s_cached_properties[m_cp->firstColour()][1].size(); ++i)
    {
      if(typeid(ParticleCacheDensity0Oc) == typeid(Particle::s_cached_properties[m_cp->firstColour()][i]))
      {
        throw gError("FunctionCallback::operator()", "ParticleCacheDensity0Oc found in stage 1 !!!");
      }
    }
    for(size_t i = 0; i < Particle::s_cached_properties[m_cp->firstColour()][2].size(); ++i)
    {
      if(typeid(ParticleCacheVolumeSelfContribution) == typeid(Particle::s_cached_properties[m_cp->firstColour()][i]))
      {
        throw gError("FunctionCallback::operator()", "ParticleCacheVolumeSelfContribution found in stage 2 !!!");
      }
    }
        
        // second colour    
    for(size_t i = 0; i < Particle::s_cached_properties[m_cp->secondColour()][1].size(); ++i)
    {
      if(typeid(ParticleCacheDensity0Oc) == typeid(Particle::s_cached_properties[m_cp->secondColour()][i]))
      {
        throw gError("FunctionCallback::operator()", "ParticleCacheDensity0Oc found in stage 1 !!!");
      }
    }
    for(size_t i = 0; i < Particle::s_cached_properties[m_cp->secondColour()][2].size(); ++i)
    {
      if(typeid(ParticleCacheVolumeSelfContribution) == typeid(Particle::s_cached_properties[m_cp->secondColour()][i]))
      {
        throw gError("FunctionCallback::operator()", "ParticleCacheVolumeSelfContribution found in stage 2 !!!");
      }
    }

    
    // is it the interpolation value of the kernel?
    // we consider this to be a non-ColourPair-specific quantity
    if (e[0] == 'W' && (e[1] == 'i' && e[2] == 'j')) {

/*MSG_DEBUG("FunctionPairCallback::operator()", "found special: " << e);
      if(e.size() == 2) MSG_DEBUG("FunctionPairCallback::operator()", "2");
      if(m_wf == NULL) MSG_DEBUG("FunctionPairCallback::operator()", "NULL");*/
      
      /* Request for a density. */
      if (e.size() == 3 && m_wf != NULL) {
// MSG_DEBUG("FunctionPairCallback::operator()", "found special, size=2: " << e);
        size_t kernel_offset;
        
        // we consider this to be a non-ColourPair-specific local density
        // first we register the value in the CPs != m_cp
        FOR_EACH_COLOUR_PAIR
            (
            M_MANAGER,
          // if this is m_cp then do nothing (will be done afterwards)
        if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
          cp->registerCalc(kernel_offset, new ValCalculatorKernel(m_wf), true);
            );
        m_cp->registerCalc(kernel_offset, new ValCalculatorKernel(m_wf), true);
      
        return new FPScalarVariable
            (e,
             NULL,
             FunctionArbitrary::double2CExpression("pair_tag", kernel_offset));
      } 
      else if (e.size() > 3 && e[3] == '[' && e[e.size()-1] == ']') {
        string wf_name = string(e, 4, e.size()-5);
        WeightingFunction *wf;
        size_t kernel_offset;
      
        cout << "wf = " << wf_name << endl;
      
        wf = M_SIMULATION->findWeightingFunction(wf_name);
      
        if (!wf)
          throw gError
              ("FunctionPairCallback::operator()",
               "Unknown weighting function: '" + wf_name + "'.");
      
        // we consider this to be a non-ColourPair-specific local density
        // first we register the value in the CPs != m_cp
        FOR_EACH_COLOUR_PAIR
            (
            M_MANAGER,
          // if this is m_cp then do nothing (will be done afterwards)
        if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
          cp->registerCalc(kernel_offset, new ValCalculatorKernel(wf), true);
            );
        m_cp->registerCalc(kernel_offset, new ValCalculatorRho(wf), true);
      
        return new FPScalarVariable
            (e,
             NULL,
             FunctionArbitrary::double2CExpression("pair_tag", kernel_offset));
      }
    }

    /* Don't know either. Perhaps it's a number? */
    return NULL;
    
#endif
    
    }
};


void FunctionPairFixed::compile()
{
//   MSG_DEBUG("FunctionPairFixed::compile", "m_default_wf COMES1" << endl << "m_expression = " << m_expression);
//   MSG_DEBUG("FunctionPairFixed::compile", "m_default_wf1 = " << m_default_wf->name());
//   MSG_DEBUG("FunctionPairFixed::compile", "m_default_wf ENDS1");
  MSG_DEBUG("FunctionPairFixed::compile", "START: m_expression = " << m_expression << ", m_cp = " << m_cp);
  FunctionArbitrary::compile();

  assert(m_cp);

  /*
   * Particle-index
   */
  addInt("i", "first_particle", offsetof(Particle, mySlot));
  addInt("j", "second_particle", offsetof(Particle, mySlot));

  /*
  * Relative distance vector
  */
  addPoint(/*"[rij]"*/"rij", "dist", offsetof(dist_t, cartesian));
  /*
  * Absolute relative distance
  */
  addDouble("rij", "dist", offsetof(dist_t, abs));
  /*
  * Position
  */
  addPoint("ri"/*"[ri]"*/, "first_particle", offsetof(Particle, r));
  addPoint(/*"[rj]"*/"rj", "second_particle", offsetof(Particle, r));
  
  /*
  * Velocity
  */
  addPoint("vi"/*"[vi]"*/, "first_particle", offsetof(Particle, v));
  addPoint(/*"[vj]"*/"vj", "second_particle", offsetof(Particle, v));

  /*
  * Add all that's found in the pair tag
  */
  addAllFromDataFormat("pair_tag", &m_cp->tagFormat(), "ij");  

  /*
  * Add all that's found in the particle tag
  */
  addAllFromDataFormat
      ("first_tag", 
       Particle::s_tag_format[/*m_colour*/m_cp->firstColour()].format(), "i");
  addAllFromDataFormat
      ("second_tag", 
       Particle::s_tag_format[/*m_colour*/m_cp->secondColour()].format(), "j");
// MSG_DEBUG("FunctionPairFixed::compile", "m_default_wf COMES10" << endl << "m_expression = " << m_expression);
// MSG_DEBUG("FunctionPairFixed::compile", "m_default_wf10 = " << m_default_wf->name());
// MSG_DEBUG("FunctionPairFixed::compile", "m_default_wf ENDS10");
  m_parser.setCallback(new FunctionPairFixedCallback(m_cp/*, m_default_wf*/));

  // here starts the "FunctionFixed-part"
  string header = "";

  m_used_vars.clear();

  if (m_vars.empty())
//     m_vars.push_back("x");
    throw gError("FunctionPairFixed::compile", "for expresssion " + m_expression + ": My list of variables is empty because my client didn't add any. Should not happen. Contact the programmer.");
  
  // add the fixed variables defined by the client with addVartiable
  for (vector<string>::iterator i = m_vars.begin(); i != m_vars.end(); ++i) {
    if((*i)[0] == '[' && (*i)[i->size()-1] == ']')
    {
      string name(*i);
      // remove the first bracket
      name.erase(0, 1);
      // remove the last bracket; don't know why, but with these arguments it works
      name.erase(name.size()-1, name.size()-1);
      
      addPoint(name, name, 0);
/*      string str_vec[SPACE_DIMS];
      point2CExpression(str_vec, *i, 0);
      m_parser.addSymbol
      (new FPVectorVariable(*i, NULL, str_vec));*/
      header += ", double* " + name;
    }
    else if((*i)[0] == '{' && (*i)[i->size()-1] == '}')
    {
      string name(*i);
      // remove the first bracket
      name.erase(0, 1);
      // remove the last bracket; don't know why, but with these arguments it works
      name.erase(name.size()-1, name.size()-1);
      
      addTensor(name, name, 0);
/*      string str_tensor[SPACE_DIMS][SPACE_DIMS];

      tensor2CExpression(str_tensor, *i, 0);
      m_parser.addSymbol
          (new FPTensorVariable
          (*i, NULL, str_tensor));*/
      header += ", double* " + name;
    }
    else if((*i)[0] != '{' && (*i)[i->size()-1] != '}' && (*i)[0] != '[' && (*i)[i->size()-1] != ']')
    {
      addDouble(*i, *i, 0);
      header += ", double " + *i;
    }
    else
      throw gError("FunctionPairFixed::compile", "for expresssion " + m_expression + ": Invalid variable " + *i + ". Contact the programmer.");

    //     m_parser.addSymbol(new FPScalarVariable(*i, NULL));
  }

  
  // here ends the "FunctionFixed-part"
  
  
  m_parser.parse(m_expression);
  
  
  m_compiler.setHeader
      (string("void *result") + header + ", void *dist, void *pair_tag, "
      "void *first_particle, void *second_particle, void *first_tag, void *second_tag");
  m_compiler.setParserAndCompile(&m_parser);
}

int FunctionPairFixed::addVariable(string name)
{
//   MSG_DEBUG("FunctionFixed::addVariable", "m_expression = " << m_expression);
  m_vars.push_back(name);
//   MSG_DEBUG("FunctionFixed::addVariable", "added variable " << name << ", size now = " << m_vars.size());

  return m_vars.size()-1;
}

void FunctionPairFixed::deleteVariable(string name)
{
  vector<string>::iterator i;

  i = find(m_vars.begin(), m_vars.end(), name);
  if (i != m_vars.end())
    m_vars.erase(i);
  else
    throw gError("FuntionPairFixed::deleteVariable", "Variable '" + name + "' not found.");
}

