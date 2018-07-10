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



#include "function_particle.h"

#include "manager_cell.h"
#include "val_calculator_rho.h"
// #include "val_calculator_shear.h"
#include "weighting_function.h"

#include "fp_vector.h"
#include "fp_tensor.h"

// #define M_MANAGER m_cp->manager()
#define M_PHASE M_MANAGER->phase()
#define M_SIMULATION ((Simulation*) M_PHASE->parent())

/* FunctionParticle */

FunctionParticle::FunctionParticle()
{
}


FunctionParticle::~FunctionParticle()
{
}



class FunctionParticleCallback: public UnknownSymbolCallback
{
protected:
//   ColourPair *m_cp;
//   WeightingFunction *m_wf;
  size_t m_colour;

public:
  FunctionParticleCallback(/*ColourPair *cp,*/ /*WeightingFunction *wf,*/ size_t colour): /*m_cp(cp),*/ /*m_wf(wf),*/ m_colour(colour) { }

  virtual TypedValue *operator()(const string &e) {
    /* Now let's see... */

//     throw gError("FunctionParticleCallback::operator()", "Test-exception: this function shouldn't be called any more.");
    
    return NULL;
    
#if 0
    
//     MSG_DEBUG("FunctionParticleCallback::()","Request for symbol: " << e);

    // non-ColourPair-specific strain rate using interpolated density?
    if (e[0] == '{' && e[1] == 'D' && e[2] == '}') {
      /* Request for a strain rate. */

      string tmpString[SPACE_DIMS][SPACE_DIMS];
      
      if (e.size() == 3 && m_wf != NULL) {
        pair<size_t, size_t> shear_offset;
      
        // we consider this to be a non-ColourPair-specific strain rate
        // first we register the value in the CPs != m_cp, if m_oneProp = true
        FOR_EACH_COLOUR_PAIR
            (
            M_MANAGER,
      // if this is m_cp then do nothing (will be done afterwards)
        if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
        {
          cp->registerCalc(shear_offset, new ValCalculatorShear(m_wf, "none"), true);
        }
            );
            m_cp->registerCalc(shear_offset, new ValCalculatorShear(m_wf, "none"), true);
      
            if(m_colour == m_cp->firstColour())
            {  
              FunctionArbitrary::tensor2CExpression(tmpString, "particle_tag", shear_offset.first);
              return new FPTensorVariable
                  (e,
                   NULL,
//                FunctionArbitrary::tensor2CExpression("particle_tag", shear_offset.first)
                   tmpString
                  );
            }
            else
            {
              FunctionArbitrary::tensor2CExpression(tmpString, "particle_tag", shear_offset.second);
              return new FPTensorVariable
                  (e,
                   NULL,
//                FunctionArbitrary::tensor2CExpression("particle_tag", shear_offset.second)
                   tmpString
                  );
            }
      } 
      else if (e.size() > 3 && e[3] == '[' && e[e.size()-1] == ']') {
        string wf_name = string(e, 4, e.size()-5);
        WeightingFunction *wf;
        pair<size_t, size_t> shear_offset;
      
        cout << "wf = " << wf_name << endl;
      
        wf = M_SIMULATION->findWeightingFunction(wf_name);
      
        if (!wf)
          throw gError
              ("FunctionParticleCallback::operator()",
               "Unknown weighting function: '" + wf_name + "'.");
      
        // first we register the value in the CPs != m_cp, if m_oneProp = true
        FOR_EACH_COLOUR_PAIR
            (
            M_MANAGER,
      // if this is m_cp then do nothing (will be done afterwards)
        if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
        {
          cp->registerCalc(shear_offset, new ValCalculatorShear(wf, "none"), true);
        }
            );
            m_cp->registerCalc(shear_offset, new ValCalculatorShear(wf, "none"), true);
      
            if(m_colour == m_cp->firstColour())
            {
              FunctionArbitrary::tensor2CExpression(tmpString, "particle_tag", shear_offset.first);
              return new FPTensorVariable
                  (e,
                   NULL,
                   /*FunctionArbitrary::tensor2CExpression("particle_tag", shear_offset.first)*/
                   tmpString);
            }
            else
            {
              FunctionArbitrary::tensor2CExpression(tmpString, "particle_tag", shear_offset.second);
              return new FPTensorVariable
                  (e,
                   NULL,
                   tmpString
               /*FunctionArbitrary::tensor2CExpression("particle_tag", shear_offset.second)*/);
            }
      }
    }

    // non-ColourPair-specific strain rate using density from scalar DOF?
    if (e[0] == '{' && e[1] == 'D' && e[2] == '_' && e[e.size()-1] == '}') {
      /* Request for a strain rate. */

      string tmpString[SPACE_DIMS][SPACE_DIMS];
      
      // Cut the scalar name out of {D_scalar}
      string scalar_name = string(e, 3, e.size()-4);
      
        pair<size_t, size_t> shear_offset;
      
        // we consider this to be a non-ColourPair-specific strain rate
        // first we register the value in the CPs != m_cp, if m_oneProp = true
        FOR_EACH_COLOUR_PAIR
            (
            M_MANAGER,
      // if this is m_cp then do nothing (will be done afterwards)
        if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
        {
          cp->registerCalc(shear_offset, new ValCalculatorShear(m_wf, scalar_name), true);
        }
            );
            m_cp->registerCalc(shear_offset, new ValCalculatorShear(m_wf, scalar_name), true);
      
            if(m_colour == m_cp->firstColour())
            {  
              FunctionArbitrary::tensor2CExpression(tmpString, "particle_tag", shear_offset.first);
              return new FPTensorVariable
                  (e,
                  NULL,
//                FunctionArbitrary::tensor2CExpression("particle_tag", shear_offset.first)
                  tmpString
                  );
            }
            else
            {
              FunctionArbitrary::tensor2CExpression(tmpString, "particle_tag", shear_offset.second);
              return new FPTensorVariable
                  (e,
                  NULL,
//                FunctionArbitrary::tensor2CExpression("particle_tag", shear_offset.second)
                  tmpString
                  );
            }
    }

    
    if (e[0] == 'n') {
      /* Request for a density. */

      if (e.size() == 1 && m_wf != NULL) {
        pair<size_t, size_t> density_offset;
        // first we register the value in the CPs != m_cp, if m_oneProp = true
        FOR_EACH_COLOUR_PAIR
            (
          M_MANAGER,
      // if this is m_cp then do nothing (will be done afterwards)
          if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
          {
            cp->registerCalc(density_offset, new ValCalculatorRho(m_wf, ""), true);
          }
            );
      
        m_cp->registerCalc(density_offset, new ValCalculatorRho(m_wf, ""), true);
      
        if(m_colour == m_cp->firstColour())
          return new FPScalarVariable
            (e,
            NULL,
            FunctionArbitrary::double2CExpression("particle_tag", density_offset.first));
        else
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("particle_tag", density_offset.second));
      } 
      else if (e.size() > 1 && e[1] == '[' && e[e.size()-1] == ']') {
        string wf_name = string(e, 2, e.size()-3);
        WeightingFunction *wf;
        pair<size_t, size_t> density_offset;
      
        cout << "wf = " << wf_name << endl;
      
        wf = M_SIMULATION->findWeightingFunction(wf_name);
      
        if (!wf)
          throw gError
            ("FunctionParticleCallback::operator()",
            "Unknown weighting function: '" + wf_name + "'.");
      
        // first we register the value in the CPs != m_cp, if m_oneProp = true
        FOR_EACH_COLOUR_PAIR
        (
          M_MANAGER,
      // if this is m_cp then do nothing (will be done afterwards)
          if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
          {
            cp->registerCalc(density_offset, new ValCalculatorRho(wf, ""), true);
          }
        );
        m_cp->registerCalc(density_offset, new ValCalculatorRho(wf, ""), true);
      
        if(m_colour == m_cp->firstColour())
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("particle_tag", density_offset.first));
        else
          return new FPScalarVariable
              (e,
               NULL,
               FunctionArbitrary::double2CExpression("particle_tag", density_offset.second));
      }
    }

    /* Don't know either. Perhaps it's a number? */
    return NULL;
  
#endif
  
  }
};



void FunctionParticle::compile()
{
  FunctionArbitrary::compile();

//   assert(m_cp);

//   string species(m_cp->manager()->species(m_colour));
// 
//   assert(m_cp->firstColour() == m_cp->secondColour());

//   assert(m_cp->firstColour() == m_colour || m_cp->secondColour() == m_colour);
  
  /*
   * Particle-index
   */
  addInt("i", "particle", offsetof(Particle, mySlot));
  
  /*
   * Position
   */
  addPoint(/*"[r]"*/"r", "particle", offsetof(Particle, r));
  
  /*
   * Velocity
   */
  addPoint(/*"[v]"*/"v", "particle", offsetof(Particle, v));

  /*
   * Add all that's found in the tag
   */
  addAllFromDataFormat("particle_tag", Particle::s_tag_format[m_colour].format());

  m_parser.setCallback(new FunctionParticleCallback(/*m_cp,*//* m_default_wf,*/ m_colour));
  m_parser.parse(m_expression);

  m_compiler.setHeader("void *result, void *particle, void *particle_tag");
  m_compiler.setParserAndCompile(&m_parser);
}

