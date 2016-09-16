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


#include "threads.h"
#include "particle.h"
#include "simulation.h"
#include "triplet_calc_angular_dt2f.h"
#include "triplet.h"

#include "triplet_calculator.h"
#include "quintet_calculator.h"

#include "particle_cache.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
// #define PI 3.141592654
#define PI M_PI

const SymbolRegister<TripletCalcAngularDt2F> triplet_calc_angular_dt2f("TripletCalcAngularDt2F");

TripletCalcAngularDt2F::TripletCalcAngularDt2F(Simulation *simulation): TripletCalculator(simulation)
{
  // does not depend on other symbols
  m_stage = 0;

  m_datatype = DataFormat::POINT;
  
  init();
#ifdef _OPENMP
  m_particleCalculator = true;
#endif
}


void TripletCalcAngularDt2F::init()
{
  m_properties.setClassName("TripletCalcPart");
  m_properties.setName("TripletCalcAngularDt2F");

  m_properties.setDescription( 
			      "This TripletCalculator computes the first two time derivatives of a hookian angular force: F=k(cos(theta)-cos(theta_eq)) and stores them in the Symbols given by the attributes 'symbol' and 'symbol2', respectively.\n Note:\n This Calculator requires forces stored in a symbol as input. These are meant to represent the total forces on the respective particle to be computed before the call of this Calculator. Note further that these forces are also used to compute velocities v'=v~+dt*f/m/2, where v~ is the 'predicted' velocity, f the previously mentioned force, and m the mass of the respective particle. This assumes a call of this Calculator at stage = \"1\" and a previous predictor integration step v~=v+dt*f'/m/2, where v is the previous velocity and f' the previous force. Other use-cases are not guaranteed to work properly and most likely require a generalisation.");

  BOOLPC(twoD, m_2D, "Should we assume 2D? Influences only the computation of the expression 1-eijeij, where 1 is the unit-matrix and eijeij is the tensor-product of the unit-vector eij with itself.");

  m_2D = false;

  ////////// FORCES NOT DONE HERE, JUST THE DERIVATIVES!!! ///////////
//   STRINGPC
//       (dtFsymbol, m_DtFname,
//        "Name of the symbol for the calculated first force derivative.");

//   m_DtFname = "undefined";


  STRINGPC
      (forceSymbol, m_FinName,
       "Name of the symbol for the calculated second force derivative.");

  m_FinName = "undefined";

  STRINGPC
      (symbol2, m_Dt2Fname,
       "Name of the symbol for the calculated second force derivative.");

  m_Dt2Fname = "undefined";


  DOUBLEPC(K, m_k, -1,"Spring force constant.");
  DOUBLEPC(thetaEq, m_thetaEq, -1,"equilibrium angle (in degrees).");
  
  m_k = 0;
  m_thetaEq = 0;


}


void TripletCalcAngularDt2F::setup()
{
  TripletCalculator::setup();

  ////////// FORCES NOT DONE HERE, JUST THE DERIVATIVES!!! ///////////
//   if(m_DtFname == "undefined")
//     throw gError("TripletCalculator::setup", "Attribute 'dtFsymbol' undefined!");
  if(m_Dt2Fname == "undefined")
    throw gError("TripletCalculator::setup", "Attribute 'dt2Fsymbol' undefined!");
  
  if(m_FinName == "undefined")
    throw gError("TripletCalculator::setup", "Attribute 'forceSymbol' undefined!");

  // the force-input symbol must already exist
  for(size_t p = 0; p < 3; ++p) {
    try
      {
	m_slotsFin[p] = Particle::s_tag_format[M_MANAGER->getColour(m_species[p])].indexOf(m_FinName, m_datatype);
	m_slotsFin[p] = Particle::s_tag_format[M_MANAGER->getColour(m_species[p])].offsetByIndex(m_slotsFin[p]);
      }
    catch(gError& err)
      {
	throw gError("TripletCalculator::setup", ": search for symbol for species '" + m_species[p] + " failed. The message was " + err.message()); 
      }
  } // end for(size_t p = 0; p < 3; ++p)  

  if(m_overwrite)
    {
      
      for(size_t p = 0; p < 3; ++p) {
  ////////// FORCES NOT DONE HERE, JUST THE DERIVATIVES!!! ///////////
// 	try
// 	  {
// 	    m_slotsDtF[p] = Particle::s_tag_format[M_MANAGER->getColour(m_species[p])].indexOf(m_DtFname, m_datatype);
// 	    m_slotsDtF[p] = Particle::s_tag_format[M_MANAGER->getColour(m_species[p])].offsetByIndex(m_slotsDtF[p]);
// 	  }
// 	catch(gError& err)
// 	  {
// 	    throw gError("TripletCalculator::setup", ": search for symbol for species '" + m_species[p] + " failed. The message was " + err.message()); 
// 	  }
	
	try
	  {
	    m_slotsDt2F[p] = Particle::s_tag_format[M_MANAGER->getColour(m_species[p])].indexOf(m_Dt2Fname, m_datatype);
	    m_slotsDt2F[p] = Particle::s_tag_format[M_MANAGER->getColour(m_species[p])].offsetByIndex(m_slotsDt2F[p]);
	  }
	catch(gError& err)
	  {
	    throw gError("TripletCalculator::setup", ": search for symbol for species '" + m_species[p] + " failed. The message was " + err.message()); 
	  }
	
      } // end for(size_t p = 0
      
    } // end if(m_overwrite)
  else {
    
    // see CONVENTION5 for rule about persistencies
    ////////// FORCES NOT DONE HERE, JUST THE DERIVATIVES!!! ///////////
//     m_slotsDtF[0] = Particle::s_tag_format[M_MANAGER->getColour(m_species[0])].addAttribute(m_DtFname, m_datatype, false, m_DtFname).offset;
    m_slotsDt2F[0] = Particle::s_tag_format[M_MANAGER->getColour(m_species[0])].addAttribute(m_Dt2Fname, m_datatype, false, m_Dt2Fname).offset;
    
    if(m_species[0] != m_species[1])
      {
        // see CONVENTION5 for rule about persistencies
    ////////// FORCES NOT DONE HERE, JUST THE DERIVATIVES!!! ///////////
// 	m_slotsDtF[1] = Particle::s_tag_format[M_MANAGER->getColour(m_species[1])].addAttribute(m_DtFname, m_datatype, false, m_DtFname).offset;
	m_slotsDt2F[1] = Particle::s_tag_format[M_MANAGER->getColour(m_species[1])].addAttribute(m_Dt2Fname, m_datatype, false, m_Dt2Fname).offset;
      }
    else {
    ////////// FORCES NOT DONE HERE, JUST THE DERIVATIVES!!! ///////////
//       m_slotsDtF[1] = m_slotsDtF[0];
      m_slotsDt2F[1] = m_slotsDt2F[0];
    }
    if(m_species[2] == m_species[0]) {
    ////////// FORCES NOT DONE HERE, JUST THE DERIVATIVES!!! ///////////
//       m_slotsDtF[2] = m_slotsDtF[0];
      m_slotsDt2F[2] = m_slotsDt2F[0];
    }
    else if(m_species[2] == m_species[1]) {
    ////////// FORCES NOT DONE HERE, JUST THE DERIVATIVES!!! ///////////
//       m_slotsDtF[2] = m_slotsDtF[1];
      m_slotsDt2F[2] = m_slotsDt2F[1];
    }
    else {
      // see CONVENTION5 for rule about persistencies
    ////////// FORCES NOT DONE HERE, JUST THE DERIVATIVES!!! ///////////
//       m_slotsDtF[2] = Particle::s_tag_format[M_MANAGER->getColour(m_species[2])].addAttribute(m_DtFname, m_datatype, false, m_DtFname).offset;
      m_slotsDt2F[2] = Particle::s_tag_format[M_MANAGER->getColour(m_species[2])].addAttribute(m_Dt2Fname, m_datatype, false, m_Dt2Fname).offset;
    
    }
    
  } // end else of if(m_overwrite)
  
  
  //   m_cosEq = cos(PI-PI*m_thetaEq/180.);	
  m_cosEq = cos(PI*m_thetaEq/180.);	
}

void TripletCalcAngularDt2F::setupAfterParticleCreation()
{
  TripletCalculator::setupAfterParticleCreation();
}

#ifndef _OPENMP
void TripletCalcAngularDt2F::compute(triplet_t* tr) {
#else
  void TripletCalcAngularDt2F::compute(triplet_t* tr/*, size_t thread_no*/ /*FIXME: paralelise!*/) {
#endif
    // FIXME!!!: currently not parallelised
    point_t eij, ekj; // vector
    double rijSq = 0;
    double rijDrkj = 0;
    double rkjSq = 0; // scalar product
    double cosa; //the cosine of the angle
    // scalar component of force
    double F;
    // the forces and later their derivatives
    point_t forcei, forcej, forcek;

    Particle* pi = tr->a;
    Particle* pj = tr->b;
    Particle* pk = tr->c;

    Data& piTag = pi->tag;
    Data& pjTag = pj->tag;
    Data& pkTag = pk->tag;

    point_t temp_p;
    tensor_t temp_t;
    
    /*---------- START: COMPUTATIONS NECESSARY FOR FORCE ---------*/

    // b is the central atom
    for (int _i = 0; _i < SPACE_DIMS; _i++) {
      // here we compute rij, rkj vectors first
      eij[_i] = pi->r[_i] - pj-> r[_i]; // rij=ri-rj
      ekj[_i] = pk->r[_i] - pj-> r[_i]; // rkj=rk-rj
      // periodic BCs
      if(m_periodic) {
	if(eij[_i] > 0.5*m_boxSize[_i]) eij[_i] -= m_boxSize[_i]; 
	if(eij[_i] < -0.5*m_boxSize[_i]) eij[_i] += m_boxSize[_i]; 
	if(ekj[_i] > 0.5*m_boxSize[_i]) ekj[_i] -= m_boxSize[_i]; 
	if(ekj[_i] < -0.5*m_boxSize[_i]) ekj[_i] += m_boxSize[_i]; 
      }      
      
      // eij, ekj still contain rij, rkj
      rijSq += eij[_i] * eij[_i]; // rij.rij (dot-product)
      rijDrkj += eij[_i] * ekj[_i]; // rij.rkj
      rkjSq += ekj[_i] * ekj[_i]; // (rkj.rkj (dot-product)
    }

    // NOW we compute eij, ekj
    double invRijAbs = 1/sqrt(rijSq);
    double invRkjAbs = 1/sqrt(rkjSq);

    // short form
//     eij *= invRijAbs;
//     ekj *= invRkjAbs;
        
    // long (but faster?) form 
    g_pointTimesEqDbl(eij, invRijAbs);
    g_pointTimesEqDbl(ekj, invRkjAbs);

    // !!! this is the REAL bond angle !!! 
    cosa = rijDrkj*invRijAbs*invRkjAbs; // = rij.rkj/(|rij|*|rkj|)=eij.ekj
    
//            MSG_DEBUG("TripletCalcAngularDt2F::compute", "particle " << pj->mySlot << " m_thetaEq = " << m_thetaEq);
    //       MSG_DEBUG("TripletCalcAngularDt2F::compute", "particle " << pj->mySlot << " cosEQ = " << cos(PI*m_thetaEq/180));
//             MSG_DEBUG("TripletCalcAngularDt2F::compute", "particle " << pj->mySlot << " m_cosEq = " << m_cosEq);


// MSG_DEBUG("TripletCalcAngularDt2F::compute", "particle " << pj->mySlot << " cos(theta) = " << cosa);

    F = m_k*(cosa-m_cosEq);

    //         MSG_DEBUG("TripletCalcAngularDt2F::compute", "particle " << pi->mySlot << " F " << F);

    // FIRST DERIVATIVES OF COSA

    // short form
//     point_t dcadri = invRijAbs*(ekj-cosa*eij);
//     point_t dcadrk = invRkjAbs*(eij-cosa*ekj);

    // long (but faster?) form 
    point_t dcadri;
    point_t dcadrk;

    g_pointTimesDbl(temp_p, eij, -cosa);
    g_pointPlusEqPoint(temp_p, ekj);
    g_pointTimesDbl(dcadri, temp_p, invRijAbs);

    g_pointTimesDbl(temp_p, ekj, -cosa);
    g_pointPlusEqPoint(temp_p, eij);
    g_pointTimesDbl(dcadrk, temp_p, invRkjAbs);


    ////////// FORCES NOT DONE HERE, JUST THE DERIVATIVES!!! ///////////

//     // FORCES
    
//     // short form
// //     forcei = -F*dcadri; 
// //     forcek = -F*dcadrk;
// //     forcej = -forcei - forcek;

//     // long (but faster?) form 
//     g_pointTimesDbl(forcei, dcadri, -F);

//     g_pointTimesDbl(forcek, dcadrk, -F);

//     g_pointPlusPoint(forcej, forcei, forcek);
//     g_pointTimesEqInt(forcej, -1);

//     //             MSG_DEBUG("TripletCalcAngularDt2F::compute", "particle " << pi->mySlot << " forcei " << forcei);
//     //             MSG_DEBUG("TripletCalcAngularDt2F::compute", "particle " << pk->mySlot << " forcek " << forcek);
//     //             MSG_DEBUG("TripletCalcAngularDt2F::compute", "particle " << pj->mySlot << " forcej " << forcej);

//     // short form
// //     piTag.pointByOffset(m_slots[0]) += forcei;
// //     pjTag.pointByOffset(m_slots[1]) += forcej;
// //     pkTag.pointByOffset(m_slots[2]) += forcek;
    
// //     MSG_DEBUG("TripletCalcAngularDt2F::compute", "particle " << pi->mySlot << " slot-forcei BEFORE = " << piTag.pointByOffset(m_slots[0]));

//     // long (but faster?) form 
//     g_pointPlusEqPoint(piTag.pointByOffset(m_slots[0]), forcei);
//     g_pointPlusEqPoint(pjTag.pointByOffset(m_slots[1]), forcej);
//     g_pointPlusEqPoint(pkTag.pointByOffset(m_slots[2]), forcek);

// //     MSG_DEBUG("TripletCalcAngularDt2F::compute", "particle " << pi->mySlot << " slot-forcei AFTER = " << piTag.pointByOffset(m_slots[0]));

    
    /*---------- END: COMPUTATIONS NECESSARY FOR FORCE ---------*/

    
    /*---------- START: COMPUTATIONS NECESSARY FOR Dt1FORCE ---------*/

    // FIRST DERIVATIVES OF eij, ekj

    // short form (Id not defined for now)
//     tensor_t deijdri = invRijAbs*(Id-eij.outer(eij));
//     tensor_t dekjdrk = invRkjAbs*(Id-ekj.outer(ekj));

    // long (but faster?) form 
    tensor_t deijdri;
    tensor_t dekjdrk;

    g_outerPointPoint(deijdri, eij, eij);
    if(m_2D)
      g_idMinusEqTensor2D(deijdri);
    else
      g_idMinusEqTensor(deijdri);
    g_tensorTimesEqDbl(deijdri, invRijAbs);

    g_outerPointPoint(dekjdrk, ekj, ekj);
    if(m_2D)
      g_idMinusEqTensor2D(dekjdrk);
    else
      g_idMinusEqTensor(dekjdrk);
    g_tensorTimesEqDbl(dekjdrk, invRkjAbs);

    // SECOND DERIVATIVES OF COSA

    // short form
//     tensor_t d2cadridri = -invRijAbs*(eij.outer(dcadri)+dcadri.outer(eij)+cosa*deijdri);
//     tensor_t d2cadrkdrk = -invRkjAbs*(ekj.outer(dcadrk)+dcadrk.outer(ekj)+cosa*dekjdrk);
//     tensor_t d2cadrkdri = invRijAbs*(dekjdrk-eij.outer(dcadrk));
//     tensor_t d2cadridrk = invRkjAbs*(deijdri-ekj.outer(dcadri));

    // long (but faster?) form 
    tensor_t d2cadridri; 
    tensor_t d2cadrkdrk;
    tensor_t d2cadrkdri;
    tensor_t d2cadridrk;

    g_tensorTimesDbl(d2cadridri, deijdri, cosa);
    g_outerPointPoint(temp_t, eij, dcadri);
    g_transpose(d2cadridri, temp_t);
    g_tensorPlusEqTensor(d2cadridri, temp_t);
    g_tensorTimesEqDbl(d2cadridri, -invRijAbs);

    g_tensorTimesDbl(d2cadrkdrk, dekjdrk, cosa);
    g_outerPointPoint(temp_t, ekj, dcadrk);
    g_transpose(d2cadrkdrk, temp_t);
    g_tensorPlusEqTensor(d2cadrkdrk, temp_t);
    g_tensorTimesEqDbl(d2cadrkdrk, -invRkjAbs);

    g_outerPointPoint(temp_t, eij, dcadrk);
    g_tensorTimesEqInt(temp_t, -1);
    g_tensorPlusTensor(d2cadrkdri, dekjdrk, temp_t);
    g_tensorTimesEqDbl(d2cadrkdri, invRijAbs);

    g_outerPointPoint(temp_t, ekj, dcadri);
    g_tensorTimesEqInt(temp_t, -1);
    g_tensorPlusTensor(d2cadridrk, deijdri, temp_t);
    g_tensorTimesEqDbl(d2cadridrk, invRkjAbs);

    // DIRECTLY COMPUTED SECOND POTENTIAL DERIVATIVES

    // F already contains m_k!!!

    // short form
//     tensor_t d2Udridri;
//     tensor_t d2Udrkdri;
//     tensor_t d2Udridrk;
//     tensor_t d2Udrkdrk;

//     d2Udridri = m_k*dcadri.outer(dcadri)+F*d2cadridri;
//     d2Udrkdri = m_k*dcadri.outer(dcadrk)+F*d2cadrkdri;
// //     CAN WE BETTER USE A TRANSPOSE IN NEXT LINE? In the long version it is done
//     d2Udridrk = m_k*dcadrk.outer(dcadri)+F*d2cadridrk;
//     d2Udrkdrk = m_k*dcadrk.outer(dcadrk)+F*d2cadrkdrk;

    // long (but faster?) form 
    tensor_t d2Udridri;
    tensor_t d2Udrkdri;
    tensor_t d2Udridrk;
    tensor_t d2Udrkdrk;

    g_outerPointPoint(temp_t, dcadri, dcadri);
    g_tensorTimesEqDbl(temp_t, m_k);
    g_tensorTimesDbl(d2Udridri, d2cadridri, F);
    g_tensorPlusEqTensor(d2Udridri, temp_t);

    g_outerPointPoint(temp_t, dcadri, dcadrk);
    g_tensorTimesEqDbl(temp_t, m_k);
    g_tensorTimesDbl(d2Udrkdri, d2cadrkdri, F);
    g_tensorPlusEqTensor(d2Udrkdri, temp_t);

//     g_outerPointPoint(temp_t, dcadrk, dcadri);
//     g_tensorTimesEqDbl(temp_t, m_k);
    g_transpose(d2Udridrk, temp_t);
//     g_tensorTimesDbl(d2Udridrk, d2cadridrk, F);
    g_tensorTimesDbl(temp_t, d2cadridrk, F);
    g_tensorPlusEqTensor(d2Udridrk, temp_t);

    g_outerPointPoint(temp_t, dcadrk, dcadrk);
    g_tensorTimesEqDbl(temp_t, m_k);
    g_tensorTimesDbl(d2Udrkdrk, d2cadrkdrk, F);
    g_tensorPlusEqTensor(d2Udrkdrk, temp_t);

    // INDIRECTLY COMPUTED SECOND POTENTIAL DERIVATIVES

    // short form
//     tensor_t d2Udrjdri = -1*d2Udridri - d2Udrkdri;
//     tensor_t d2Udrjdrk = -1*d2Udridrk - d2Udrkdrk;

    // long (but faster?) form 
    tensor_t d2Udrjdri;
    tensor_t d2Udrjdrk;
    g_tensorPlusTensor(d2Udrjdri, d2Udridri, d2Udrkdri);
    g_tensorTimesEqInt(d2Udrjdri, -1);
    g_tensorPlusTensor(d2Udrjdrk, d2Udridrk, d2Udrkdrk);
    g_tensorTimesEqInt(d2Udrjdrk, -1);

    // FINAL FIRST DERIVATIVES

    double dt = M_CONTROLLER -> dt();

    point_t ai = piTag.pointByOffset(m_slotsFin[0]) / pi->m_mass;
    point_t aj = pjTag.pointByOffset(m_slotsFin[1]) / pj->m_mass;
    point_t ak = pkTag.pointByOffset(m_slotsFin[2]) / pk->m_mass;

    point_t vi = pi->v + dt*ai/2;
    point_t vj = pj->v + dt*aj/2;
    point_t vk = pk->v + dt*ak/2;

    // short form
//     forcei = -1*(d2Udridri*/*2.1*/vi+d2Udrjdri*/*2.1*/vj+d2Udrkdri*/*2.1*/vk);
//     forcek = -1*(d2Udridrk*/*2.1*/vi+d2Udrjdrk*/*2.1*/vj+d2Udrkdrk*/*2.1*/vk);
//     forcej = -forcei - forcek;

    // long (but faster?) form 
    g_tensorDotPoint(forcei, d2Udridri, vi);
    g_tensorDotPoint(temp_p, d2Udrjdri, vj);
    g_pointPlusEqPoint(forcei, temp_p);
    g_tensorDotPoint(temp_p, d2Udrkdri, vk);
    g_pointPlusEqPoint(forcei, temp_p);
//     g_pointTimesEqInt(forcei, -1);

    g_tensorDotPoint(forcek, d2Udridrk, vi);
    g_tensorDotPoint(temp_p, d2Udrjdrk, vj);
    g_pointPlusEqPoint(forcek, temp_p);
    g_tensorDotPoint(temp_p, d2Udrkdrk, vk);
    g_pointPlusEqPoint(forcek, temp_p);
//     g_pointTimesEqInt(forcek, -1);

    g_pointPlusPoint(forcej, forcei, forcek);
//     g_pointTimesEqInt(forcej, -1);

    g_pointTimesEqInt(forcei, -1);
    g_pointTimesEqInt(forcek, -1);


    // ADD TO TAG

    // short form
//     piTag.pointByOffset(m_slotsDtF[0]) += forcei;
//     pjTag.pointByOffset(m_slotsDtF[1]) += forcej;
//     pkTag.pointByOffset(m_slotsDtF[2]) += forcek;

    // long (but faster?) form 
    g_pointPlusEqPoint(piTag.pointByOffset(m_slots/*DtF*/[0]), forcei);
    g_pointPlusEqPoint(pjTag.pointByOffset(m_slots/*DtF*/[1]), forcej);
    g_pointPlusEqPoint(pkTag.pointByOffset(m_slots/*DtF*/[2]), forcek);
    
    /*---------- END: COMPUTATIONS NECESSARY FOR Dt1FORCE ---------*/
    
    
    /*---------- START: COMPUTATIONS NECESSARY FOR Dt2FORCE ---------*/

    // SECOND DERIVATIVES OF eij, ekj

    // short form 
    // (outer2 means the argument contributes the 2nd, (not the 3rd!) index)
//     tensor3_t d2eijdridri = 
//       -invAbsRij*(eij.outer(deijdri)+deijdri.outer(eij)+deijdri.outer2(eij));
//     tensor3_t d2eijdrkdrk = 
//       -invAbsRkj*(ekj.outer(dekjdrk)+dekjdrk.outer(ekj)+dekjdrk.outer2(ekj));

    // long (but faster?) form 
    // g_outer2TensorPoint means that second input(!) contributes to 2nd index instead of 3rd.
    tensor3_t temp_3t;
    tensor3_t d2eijdridri, d2ekjdrkdrk; 

//     if(m_listName != "angular") {
//       MSG_DEBUG("TripletCalcAngularDt2F::compute", "d2eijdridri: " << d2eijdridri);
//       MSG_DEBUG("TripletCalcAngularDt2F::compute", "eij: " << eij);
//       MSG_DEBUG("TripletCalcAngularDt2F::compute", "deijdri: " << deijdri);
//     }

    g_outerPointTensor(d2eijdridri, eij, deijdri);
    g_outerTensorPoint(temp_3t, deijdri, eij);
    g_tensor3PlusEqTensor3(d2eijdridri, temp_3t);
    g_outer2TensorPoint(temp_3t, deijdri, eij);
    g_tensor3PlusEqTensor3(d2eijdridri, temp_3t);
    g_tensor3TimesEqDbl(d2eijdridri, -invRijAbs);

//     if(m_listName != "angular")
//       MSG_DEBUG("TripletCalcAngularDt2F::compute", "d2eijdridri: " << d2eijdridri);

    g_outerPointTensor(d2ekjdrkdrk, ekj, dekjdrk);
    g_outerTensorPoint(temp_3t, dekjdrk, ekj);
    g_tensor3PlusEqTensor3(d2ekjdrkdrk, temp_3t);
    g_outer2TensorPoint(temp_3t, dekjdrk, ekj);
    g_tensor3PlusEqTensor3(d2ekjdrkdrk, temp_3t);
    g_tensor3TimesEqDbl(d2ekjdrkdrk, -invRkjAbs);


    // THIRD DERIVATIVES OF COSA

    // short form 
    // (outer2 means the argument contributes the 2nd, (not the 3rd!) index)
//     tensor3_t d3cadridridri 
//       = -invRijAbs*(d2cadridri.outer(eij)+deijdri.outer2(dcadri)
// 		    +eij.outer(d2cadridri)+d2cadridri.outer2(eij)
// 		    +dcadri.outer(deijdri)+deijdri.outer(dcadri)
// 		    +cosa*d2eijdridri
// 		    );
//     tensor3_t d3cadridridrk 
//       = invRkjAbs*(d2eijdridri-ekj.outer(d2cadridri));
//     tensor3_t d3cadridrkdri
//       = -invRijAbs*(d2cadrkdri.outer(eij)+deijdri.outer2(dcadrk)
// 		    +eij.outer(d2cadridrk)
// 		    );    
//     tensor3_t d3cadrkdridri
//       = -invRijAbs*(eij.outer(d2cadrkdri)+d2cadrkdri.outer2(eij)
// 		    +deijdri.outer(dcadrk)
// 		    ); 
//     // next four just first four with interchanged i-k
//     tensor3_t d3cadrkdrkdrk 
//       = -invRkjAbs*(d2cadrkdrk.outer(ekj)+dekjdrk.outer2(dcadrk)
// 		    +ekj.outer(d2cadrkdrk)+d2cadrkdrk.outer2(ekj)
// 		    +dcadrk.outer(dekjdrk)+dekjdrk.outer(dcadrk)
// 		    +cosa*d2ekjdrkdrk
// 		    );
//     tensor3_t d3cadrkdrkdri 
//       = invRijAbs*(d2ekjdrkdrk-eij.outer(d2cadrkdrk));
//     tensor3_t d3cadrkdridrk
//       = -invRkjAbs*(d2cadridrk.outer(ekj)+dekjdrk.outer2(dcadri)
// 		    +ekj.outer(d2cadrkdri)
// 		    );
//     tensor3_t d3cadridrkdrk
//       = -invRkjAbs*(ekj.outer(d2cadridrk)+d2cadridrk.outer2(ekj)
// 		    +dekjdrk.outer(dcadri)
// 		    ); 

    // long (but faster?) form 
    // g_outer2TensorPoint means that second input(!) contributes to 2nd index instead of 3rd.
    tensor3_t d3cadridridri, d3cadridridrk, d3cadridrkdri, d3cadrkdridri, d3cadrkdrkdrk, d3cadrkdrkdri, d3cadrkdridrk, d3cadridrkdrk;

    g_outerTensorPoint(d3cadridridri, d2cadridri, eij);
    g_outer2TensorPoint(temp_3t, deijdri, dcadri);
    g_tensor3PlusEqTensor3(d3cadridridri, temp_3t);
    g_outerPointTensor(temp_3t, eij, d2cadridri);
    g_tensor3PlusEqTensor3(d3cadridridri, temp_3t);
    g_outer2TensorPoint(temp_3t, d2cadridri, eij);
    g_tensor3PlusEqTensor3(d3cadridridri, temp_3t);
    g_outerPointTensor(temp_3t, dcadri, deijdri);
    g_tensor3PlusEqTensor3(d3cadridridri, temp_3t);
    g_outerTensorPoint(temp_3t, deijdri, dcadri);
    g_tensor3PlusEqTensor3(d3cadridridri, temp_3t);
    g_tensor3TimesDbl(temp_3t, d2eijdridri, cosa);
    g_tensor3PlusEqTensor3(d3cadridridri, temp_3t);
    g_tensor3TimesEqDbl(d3cadridridri, -invRijAbs);

    g_outerPointTensor(d3cadridridrk, ekj, d2cadridri);
    g_tensor3TimesEqInt(d3cadridridrk, -1);
    g_tensor3PlusEqTensor3(d3cadridridrk, d2eijdridri);
    g_tensor3TimesEqDbl(d3cadridridrk, invRkjAbs);

    g_outerTensorPoint(d3cadridrkdri, d2cadrkdri, eij);
    g_outer2TensorPoint(temp_3t, deijdri, dcadrk);
    g_tensor3PlusEqTensor3(d3cadridrkdri, temp_3t);
    g_outerPointTensor(temp_3t, eij, d2cadridrk);
    g_tensor3PlusEqTensor3(d3cadridrkdri, temp_3t);
    g_tensor3TimesEqDbl(d3cadridrkdri, -invRijAbs);

    g_outerPointTensor(d3cadrkdridri, eij, d2cadrkdri);
    g_outer2TensorPoint(temp_3t, d2cadrkdri, eij);
    g_tensor3PlusEqTensor3(d3cadrkdridri, temp_3t);
    g_outerTensorPoint(temp_3t, deijdri, dcadrk);
    g_tensor3PlusEqTensor3(d3cadrkdridri, temp_3t);
    g_tensor3TimesEqDbl(d3cadrkdridri, -invRijAbs);

    // next four just first four with interchanged i-k
    g_outerTensorPoint(d3cadrkdrkdrk, d2cadrkdrk, ekj);
    g_outer2TensorPoint(temp_3t, dekjdrk, dcadrk);
    g_tensor3PlusEqTensor3(d3cadrkdrkdrk, temp_3t);
    g_outerPointTensor(temp_3t, ekj, d2cadrkdrk);
    g_tensor3PlusEqTensor3(d3cadrkdrkdrk, temp_3t);
    g_outer2TensorPoint(temp_3t, d2cadrkdrk, ekj);
    g_tensor3PlusEqTensor3(d3cadrkdrkdrk, temp_3t);
    g_outerPointTensor(temp_3t, dcadrk, dekjdrk);
    g_tensor3PlusEqTensor3(d3cadrkdrkdrk, temp_3t);
    g_outerTensorPoint(temp_3t, dekjdrk, dcadrk);
    g_tensor3PlusEqTensor3(d3cadrkdrkdrk, temp_3t);
    g_tensor3TimesDbl(temp_3t, d2ekjdrkdrk, cosa);
    g_tensor3PlusEqTensor3(d3cadrkdrkdrk, temp_3t);
    g_tensor3TimesEqDbl(d3cadrkdrkdrk, -invRkjAbs);

    g_outerPointTensor(d3cadrkdrkdri, eij, d2cadrkdrk);
    g_tensor3TimesEqInt(d3cadrkdrkdri, -1);
    g_tensor3PlusEqTensor3(d3cadrkdrkdri, d2ekjdrkdrk);
    g_tensor3TimesEqDbl(d3cadrkdrkdri, invRijAbs);

    g_outerTensorPoint(d3cadrkdridrk, d2cadridrk, ekj);
    g_outer2TensorPoint(temp_3t, dekjdrk, dcadri);
    g_tensor3PlusEqTensor3(d3cadrkdridrk, temp_3t);
    g_outerPointTensor(temp_3t, ekj, d2cadrkdri);
    g_tensor3PlusEqTensor3(d3cadrkdridrk, temp_3t);
    g_tensor3TimesEqDbl(d3cadrkdridrk, -invRkjAbs);

    g_outerPointTensor(d3cadridrkdrk, ekj, d2cadridrk);
    g_outer2TensorPoint(temp_3t, d2cadridrk, ekj);
    g_tensor3PlusEqTensor3(d3cadridrkdrk, temp_3t);
    g_outerTensorPoint(temp_3t, dekjdrk, dcadri);
    g_tensor3PlusEqTensor3(d3cadridrkdrk, temp_3t);
    g_tensor3TimesEqDbl(d3cadridrkdrk, -invRkjAbs);

    // DIRECTLY COMPUTED THIRD POTENTIAL DERIVATIVES

    // F already contains m_k!!!

    // short form
//     tensor3_t d3Udridridri 
//       = m_k*(d2cadridri.outer(dcadri)+d2cadridri.outer2(dcadri)+dcadri.outer(d2cadridri))+F*d3cadridridri;
//     tensor3_t d3Udridrkdri 
//       = m_k*(d2cadridrk.outer(dcadri)+d2cadridri.outer2(dcadrk)+dcadri.outer(d2cadrkdri))+F*d3cadridrkdri;
//     tensor3_t d3Udrkdridri 
//       = m_k*(d2cadrkdri.outer(dcadri)+d2cadrkdri.outer2(dcadri)+dcadrk.outer(d2cadridri))+F*d3cadrkdridri;
//     tensor3_t d3Udrkdrkdri 
//       = m_k*(d2cadrkdrk.outer(dcadri)+d2cadrkdri.outer2(dcadrk)+dcadrk.outer(d2cadrkdri))+F*d3cadrkdrkdri;

//     tensor3_t d3Udrkdrkdrk 
//       = m_k*(d2cadrkdrk.outer(dcadrk)+d2cadrkdrk.outer2(dcadrk)+dcadrk.outer(d2cadrkdrk))+F*d3cadrkdrkdrk;
//     tensor3_t d3Udrkdridrk 
//       = m_k*(d2cadrkdri.outer(dcadrk)+d2cadrkdrk.outer2(dcadri)+dcadrk.outer(d2cadridrk))+F*d3cadrkdridrk;
//     tensor3_t d3Udridrkdrk 
//       = m_k*(d2cadridrk.outer(dcadrk)+d2cadridrk.outer2(dcadrk)+dcadri.outer(d2cadrkdrk))+F*d3cadridrkdrk;
//     tensor3_t d3Udridridrk 
//       = m_k*(d2cadridri.outer(dcadrk)+d2cadridrk.outer2(dcadri)+dcadri.outer(d2cadridrk))+F*d3cadridridrk;



    // long (but faster?) form 
    tensor3_t d3Udridridri, d3Udridrkdri, d3Udrkdridri, d3Udrkdrkdri, d3Udridridrk, d3Udridrkdrk, d3Udrkdridrk, d3Udrkdrkdrk;

    g_outerTensorPoint(d3Udridridri, d2cadridri, dcadri);
    g_outer2TensorPoint(temp_3t, d2cadridri, dcadri);
    g_tensor3PlusEqTensor3(d3Udridridri, temp_3t);
    g_outerPointTensor(temp_3t, dcadri, d2cadridri);
    g_tensor3PlusEqTensor3(d3Udridridri, temp_3t);
    g_tensor3TimesEqDbl(d3Udridridri, m_k);
    g_tensor3TimesDbl(temp_3t, d3cadridridri, F);
    g_tensor3PlusEqTensor3(d3Udridridri, temp_3t);

    g_outerTensorPoint(d3Udridrkdri, d2cadridrk, dcadri);
    g_outer2TensorPoint(temp_3t, d2cadridri, dcadrk);
    g_tensor3PlusEqTensor3(d3Udridrkdri, temp_3t);
    g_outerPointTensor(temp_3t, dcadri, d2cadrkdri);
    g_tensor3PlusEqTensor3(d3Udridrkdri, temp_3t);
    g_tensor3TimesEqDbl(d3Udridrkdri, m_k);
    g_tensor3TimesDbl(temp_3t, d3cadridrkdri, F);
    g_tensor3PlusEqTensor3(d3Udridrkdri, temp_3t);

    g_outerTensorPoint(d3Udrkdridri, d2cadrkdri, dcadri);
    g_outer2TensorPoint(temp_3t, d2cadrkdri, dcadri);
    g_tensor3PlusEqTensor3(d3Udrkdridri, temp_3t);
    g_outerPointTensor(temp_3t, dcadrk, d2cadridri);
    g_tensor3PlusEqTensor3(d3Udrkdridri, temp_3t);
    g_tensor3TimesEqDbl(d3Udrkdridri, m_k);
    g_tensor3TimesDbl(temp_3t, d3cadrkdridri, F);
    g_tensor3PlusEqTensor3(d3Udrkdridri, temp_3t);

    g_outerTensorPoint(d3Udrkdrkdri, d2cadrkdrk, dcadri);
    g_outer2TensorPoint(temp_3t, d2cadrkdri, dcadrk);
    g_tensor3PlusEqTensor3(d3Udrkdrkdri, temp_3t);
    g_outerPointTensor(temp_3t, dcadrk, d2cadrkdri);
    g_tensor3PlusEqTensor3(d3Udrkdrkdri, temp_3t);
    g_tensor3TimesEqDbl(d3Udrkdrkdri, m_k);
    g_tensor3TimesDbl(temp_3t, d3cadrkdrkdri, F);
    g_tensor3PlusEqTensor3(d3Udrkdrkdri, temp_3t);
    //
    g_outerTensorPoint(d3Udrkdrkdrk, d2cadrkdrk, dcadrk);
    g_outer2TensorPoint(temp_3t, d2cadrkdrk, dcadrk);
    g_tensor3PlusEqTensor3(d3Udrkdrkdrk, temp_3t);
    g_outerPointTensor(temp_3t, dcadrk, d2cadrkdrk);
    g_tensor3PlusEqTensor3(d3Udrkdrkdrk, temp_3t);
    g_tensor3TimesEqDbl(d3Udrkdrkdrk, m_k);
    g_tensor3TimesDbl(temp_3t, d3cadrkdrkdrk, F);
    g_tensor3PlusEqTensor3(d3Udrkdrkdrk, temp_3t);

    g_outerTensorPoint(d3Udrkdridrk, d2cadrkdri, dcadrk);
    g_outer2TensorPoint(temp_3t, d2cadrkdrk, dcadri);
    g_tensor3PlusEqTensor3(d3Udrkdridrk, temp_3t);
    g_outerPointTensor(temp_3t, dcadrk, d2cadridrk);
    g_tensor3PlusEqTensor3(d3Udrkdridrk, temp_3t);
    g_tensor3TimesEqDbl(d3Udrkdridrk, m_k);
    g_tensor3TimesDbl(temp_3t, d3cadrkdridrk, F);
    g_tensor3PlusEqTensor3(d3Udrkdridrk, temp_3t);

    g_outerTensorPoint(d3Udridrkdrk, d2cadridrk, dcadrk);
    g_outer2TensorPoint(temp_3t, d2cadridrk, dcadrk);
    g_tensor3PlusEqTensor3(d3Udridrkdrk, temp_3t);
    g_outerPointTensor(temp_3t, dcadri, d2cadrkdrk);
    g_tensor3PlusEqTensor3(d3Udridrkdrk, temp_3t);
    g_tensor3TimesEqDbl(d3Udridrkdrk, m_k);
    g_tensor3TimesDbl(temp_3t, d3cadridrkdrk, F);
    g_tensor3PlusEqTensor3(d3Udridrkdrk, temp_3t);

    g_outerTensorPoint(d3Udridridrk, d2cadridri, dcadrk);
    g_outer2TensorPoint(temp_3t, d2cadridrk, dcadri);
    g_tensor3PlusEqTensor3(d3Udridridrk, temp_3t);
    g_outerPointTensor(temp_3t, dcadri, d2cadridrk);
    g_tensor3PlusEqTensor3(d3Udridridrk, temp_3t);
    g_tensor3TimesEqDbl(d3Udridridrk, m_k);
    g_tensor3TimesDbl(temp_3t, d3cadridridrk, F);
    g_tensor3PlusEqTensor3(d3Udridridrk, temp_3t);


    // INDIRECTLY COMPUTED THIRD POTENTIAL DERIVATIVES

    // short form
//     tensor_t d3Udridrjdri = -1*d3Udridridri - d3Udridrkdri;
//     tensor_t d3Udrjdridri = -1*d3Udridridri - d3Udrkdridri;
//     tensor_t d3Udrkdrjdri = -1*d3Udrkdridri - d3Udrkdrkdri;
//     tensor_t d3Udrjdrkdri = -1*d3Udridrkdri - d3Udrkdrkdri;
//     tensor_t d3Udrjdrjdri = -1*d3Udrjdridri - d3Udrjdrkdri;

//     tensor_t d3Udridrjdrk = -1*d3Udridridrk - d3Udridrkdrk;
//     tensor_t d3Udrjdridrk = -1*d3Udridridrk - d3Udrkdridrk;
//     tensor_t d3Udrkdrjdrk = -1*d3Udrkdridrk - d3Udrkdrkdrk;
//     tensor_t d3Udrjdrkdrk = -1*d3Udridrkdrk - d3Udrkdrkdrk;
//     tensor_t d3Udrjdrjdrk = -1*d3Udrjdridrk - d3Udrjdrkdrk;

    // long (but faster?) form 
    tensor3_t d3Udridrjdri, d3Udrjdridri, d3Udrjdrjdri, d3Udrkdrjdri, d3Udrjdrkdri;
    tensor3_t d3Udridrjdrk, d3Udrjdridrk, d3Udrjdrjdrk, d3Udrkdrjdrk, d3Udrjdrkdrk;

    g_tensor3PlusTensor3(d3Udridrjdri, d3Udridridri, d3Udridrkdri);
    g_tensor3TimesEqInt(d3Udridrjdri, -1);
    g_tensor3PlusTensor3(d3Udrjdridri, d3Udridridri, d3Udrkdridri);
    g_tensor3TimesEqInt(d3Udrjdridri, -1);
    g_tensor3PlusTensor3(d3Udrkdrjdri, d3Udrkdridri, d3Udrkdrkdri);
    g_tensor3TimesEqInt(d3Udrkdrjdri, -1);
    g_tensor3PlusTensor3(d3Udrjdrkdri, d3Udridrkdri, d3Udrkdrkdri);
    g_tensor3TimesEqInt(d3Udrjdrkdri, -1);
    g_tensor3PlusTensor3(d3Udrjdrjdri, d3Udrjdridri, d3Udrjdrkdri);
    g_tensor3TimesEqInt(d3Udrjdrjdri, -1);

    g_tensor3PlusTensor3(d3Udridrjdrk, d3Udridridrk, d3Udridrkdrk);
    g_tensor3TimesEqInt(d3Udridrjdrk, -1);
    g_tensor3PlusTensor3(d3Udrjdridrk, d3Udridridrk, d3Udrkdridrk);
    g_tensor3TimesEqInt(d3Udrjdridrk, -1);
    g_tensor3PlusTensor3(d3Udrkdrjdrk, d3Udrkdridrk, d3Udrkdrkdrk);
    g_tensor3TimesEqInt(d3Udrkdrjdrk, -1);
    g_tensor3PlusTensor3(d3Udrjdrkdrk, d3Udridrkdrk, d3Udrkdrkdrk);
    g_tensor3TimesEqInt(d3Udrjdrkdrk, -1);
    g_tensor3PlusTensor3(d3Udrjdrjdrk, d3Udrjdridrk, d3Udrjdrkdrk);
    g_tensor3TimesEqInt(d3Udrjdrjdrk, -1);



    // FINAL SECOND DERIVATIVES

    ////////// DEBUGGING START ////////////////

//     g_tensor3DotPoint(temp_t, d3Udridridri, vi);

//     MSG_DEBUG("TripletCalcAngularDt2F::compute", "ai=" << ai << "\nd2Udridri=" << d2Udridri << "\nd3Udridridri.vj=" << temp_t);

    ////////// DEBUGGING END ////////////////


    // short form
//     forcei = 
//       -1*(
// 	  d2Udridri*/*2.1*/ai+d2Udrjdri*/*2.1*/aj+d2Udrkdri*/*2.1*/ak
// 	  +(
// 	    d3Udridridri*/*3.1*/vi+d3Udrjdridri*/*3.1*/vj+d3Udrkdridri*/*3.1*/vk
// 	    )*/*2.1*/vi
// 	  +(
// 	    d3Udridrjdri*/*3.1*/vi+d3Udrjdrjdri*/*3.1*/vj+d3Udrkdrjdri*/*3.1*/vk
// 	    )*/*2.1*/vj
// 	  +(
// 	    d3Udridrkdri*/*3.1*/vi+d3Udrjdrkdri*/*3.1*/vj+d3Udrkdrkdri*/*3.1*/vk
// 	    )*/*2.1*/vk	  
// 	  );
//     forcek = 
//       -1*(
// 	  d2Udridrk*/*2.1*/ai+d2Udrjdrk*/*2.1*/aj+d2Udrkdrk*/*2.1*/ak
// 	  +(
// 	    d3Udridridrk*/*3.1*/vi+d3Udrjdridrk*/*3.1*/vj+d3Udrkdridrk*/*3.1*/vk
// 	    )*/*2.1*/vi
// 	  +(
// 	    d3Udridrjdrk*/*3.1*/vi+d3Udrjdrjdrk*/*3.1*/vj+d3Udrkdrjdrk*/*3.1*/vk
// 	    )*/*2.1*/vj
// 	  +(
// 	    d3Udridrkdrk*/*3.1*/vi+d3Udrjdrkdrk*/*3.1*/vj+d3Udrkdrkdrk*/*3.1*/vk
// 	    )*/*2.1*/vk	  
// 	  );
//     forcej = -forcei - forcek;

    // long (but faster?) form 
    tensor_t temp_t_b;

    // from here forcei
    g_tensorDotPoint(forcei, d2Udridri, ai);
    g_tensorDotPoint(temp_p, d2Udrjdri, aj);
    g_pointPlusEqPoint(forcei, temp_p);
    g_tensorDotPoint(temp_p, d2Udrkdri, ak);
    g_pointPlusEqPoint(forcei, temp_p);
//     g_pointTimesEqInt(forcei, -1);

    g_tensor3DotPoint(temp_t, d3Udridridri, vi);
    g_tensor3DotPoint(temp_t_b, d3Udrjdridri, vj);
    g_tensorPlusEqTensor(temp_t, temp_t_b);
    g_tensor3DotPoint(temp_t_b, d3Udrkdridri, vk);
    g_tensorPlusEqTensor(temp_t, temp_t_b);
    g_tensorDotPoint(temp_p, temp_t, vi);
    g_pointPlusEqPoint(forcei, temp_p);

    g_tensor3DotPoint(temp_t, d3Udridrjdri, vi);
    g_tensor3DotPoint(temp_t_b, d3Udrjdrjdri, vj);
    g_tensorPlusEqTensor(temp_t, temp_t_b);
    g_tensor3DotPoint(temp_t_b, d3Udrkdrjdri, vk);
    g_tensorPlusEqTensor(temp_t, temp_t_b);
    g_tensorDotPoint(temp_p, temp_t, vj);
    g_pointPlusEqPoint(forcei, temp_p);

    g_tensor3DotPoint(temp_t, d3Udridrkdri, vi);
    g_tensor3DotPoint(temp_t_b, d3Udrjdrkdri, vj);
    g_tensorPlusEqTensor(temp_t, temp_t_b);
    g_tensor3DotPoint(temp_t_b, d3Udrkdrkdri, vk);
    g_tensorPlusEqTensor(temp_t, temp_t_b);
    g_tensorDotPoint(temp_p, temp_t, vk);
    g_pointPlusEqPoint(forcei, temp_p);

    // from here forcek
    g_tensorDotPoint(forcek, d2Udridrk, ai);
    g_tensorDotPoint(temp_p, d2Udrjdrk, aj);
    g_pointPlusEqPoint(forcek, temp_p);
    g_tensorDotPoint(temp_p, d2Udrkdrk, ak);
    g_pointPlusEqPoint(forcek, temp_p);
//     g_pointTimesEqInt(forcek, -1);

    g_tensor3DotPoint(temp_t, d3Udridridrk, vi);
    g_tensor3DotPoint(temp_t_b, d3Udrjdridrk, vj);
    g_tensorPlusEqTensor(temp_t, temp_t_b);
    g_tensor3DotPoint(temp_t_b, d3Udrkdridrk, vk);
    g_tensorPlusEqTensor(temp_t, temp_t_b);
    g_tensorDotPoint(temp_p, temp_t, vi);
    g_pointPlusEqPoint(forcek, temp_p);

    g_tensor3DotPoint(temp_t, d3Udridrjdrk, vi);
    g_tensor3DotPoint(temp_t_b, d3Udrjdrjdrk, vj);
    g_tensorPlusEqTensor(temp_t, temp_t_b);
    g_tensor3DotPoint(temp_t_b, d3Udrkdrjdrk, vk);
    g_tensorPlusEqTensor(temp_t, temp_t_b);
    g_tensorDotPoint(temp_p, temp_t, vj);
    g_pointPlusEqPoint(forcek, temp_p);

    g_tensor3DotPoint(temp_t, d3Udridrkdrk, vi);
    g_tensor3DotPoint(temp_t_b, d3Udrjdrkdrk, vj);
    g_tensorPlusEqTensor(temp_t, temp_t_b);
    g_tensor3DotPoint(temp_t_b, d3Udrkdrkdrk, vk);
    g_tensorPlusEqTensor(temp_t, temp_t_b);
    g_tensorDotPoint(temp_p, temp_t, vk);
    g_pointPlusEqPoint(forcek, temp_p);

    // from here forcej and negation of forcei/k
    g_pointPlusPoint(forcej, forcei, forcek);
//     g_pointTimesEqInt(forcej, -1);

    g_pointTimesEqInt(forcei, -1);
    g_pointTimesEqInt(forcek, -1);


    // ADD TO TAG


    ////////// DEBUGGING START ////////////////

//     if(m_listName != "angular")
//       MSG_DEBUG("TripletCalcAngularDt2F::compute", "BEFORE: list = " << m_listName << "\npforcei=" << piTag.pointByOffset(m_slotsDt2F[0]) << "\npforcej=" << pjTag.pointByOffset(m_slotsDt2F[1]) << "\npforcek=" << pkTag.pointByOffset(m_slotsDt2F[2]));

    ////////// DEBUGGING END ////////////////


    // short form
//     piTag.pointByOffset(m_slotsDt2F[0]) += forcei;
//     pjTag.pointByOffset(m_slotsDt2F[1]) += forcej;
//     pkTag.pointByOffset(m_slotsDt2F[2]) += forcek;

    // long (but faster?) form 
    g_pointPlusEqPoint(piTag.pointByOffset(m_slotsDt2F[0]), forcei);
    g_pointPlusEqPoint(pjTag.pointByOffset(m_slotsDt2F[1]), forcej);
    g_pointPlusEqPoint(pkTag.pointByOffset(m_slotsDt2F[2]), forcek);


    ////////// DEBUGGING START ////////////////

//     if(m_listName != "angular")
//       MSG_DEBUG("TripletCalcAngularDt2F::compute", "\nforcei=" << forcei << "\nforcej=" << forcej << "\nforcek=" << forcek << "\npforcei=" << piTag.pointByOffset(m_slotsDt2F[0]) << "\npforcej=" << pjTag.pointByOffset(m_slotsDt2F[1]) << "\npforcek=" << pkTag.pointByOffset(m_slotsDt2F[2]) << "\npiSlot=" << pi->mySlot << "\npjSlot=" << pj->mySlot << "\npkSlot=" << pk->mySlot << "\nri=" << pi->r << "\nrj=" << pj->r << "\nrk=" << pk->r);

    ////////// DEBUGGING END ////////////////

    
    /*---------- END: COMPUTATIONS NECESSARY FOR Dt2FORCE ---------*/
    
  }

#ifdef _OPENMP

int TripletCalcAngularDt2F::setNumOfDoubles() {
  // three point_t 
  return 9;
}

// FIXME: This module is not yet parallelised. If you parallelise, the following function could roughly do what is commented out now
void TripletCalcAngularDt2F::mergeCopies(size_t thread_no) {
//   size_t slot1 = m_slots[0];
//   size_t slot2 = m_slots[1];
//   size_t slot3 = m_slots[2];
//   size_t copySlot1 = m_copy_slots[thread_no][0];
//   size_t copySlot2 = m_copy_slots[thread_no][1];
//   size_t copySlot3 = m_copy_slots[thread_no][2];
//   size_t vecSlot1 = m_vector_slots[0];
//   size_t vecSlot2 = m_vector_slots[1];
//   size_t vecSlot3 = m_vector_slots[2];

//   size_t slotDtF1 = m_slotsDtF[0];
//   size_t slotDtF2 = m_slotsDtF[1];
//   size_t slotDtF3 = m_slotsDtF[2];
//   size_t copySlotDtF1 = m_copy_slotsDtF[thread_no][0];
//   size_t copySlotDtF2 = m_copy_slotsDtF[thread_no][1];
//   size_t copySlotDtF3 = m_copy_slotsDtF[thread_no][2];
//   size_t vecSlotDtF1 = m_vector_slotsDtF[0];
//   size_t vecSlotDtF2 = m_vector_slotsDtF[1];
//   size_t vecSlotDtF3 = m_vector_slotsDtF[2];

//   size_t slotDt2F1 = m_slotsDt2F[0];
//   size_t slotDt2F2 = m_slotsDt2F[1];
//   size_t slotDt2F3 = m_slotsDt2F[2];
//   size_t copySlotDt2F1 = m_copy_slotsDt2F[thread_no][0];
//   size_t copySlotDt2F2 = m_copy_slotsDt2F[thread_no][1];
//   size_t copySlotDt2F3 = m_copy_slotsDt2F[thread_no][2];
//   size_t vecSlotDt2F1 = m_vector_slotsDt2F[0];
//   size_t vecSlotDt2F2 = m_vector_slotsDt2F[1];
//   size_t vecSlotDt2F3 = m_vector_slotsDt2F[2];

//   FOR_EACH_PARTICLE_C 
//   (M_PHASE, m_firstColour,

//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.pointByOffset(slot1)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlot1))[vecSlot1 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlot1))[vecSlot1 + j] = 0;
//     }

//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.pointByOffset(slot2)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + j] = 0;
//     }

//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.pointByOffset(slot3)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlot3))[vecSlot3 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlot3))[vecSlot3 + j] = 0;
//     }


//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.pointByOffset(slotDtF1)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlotDtF1))[vecSlotDtF1 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlotDtF1))[vecSlotDtF1 + j] = 0;
//     }

//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.pointByOffset(slotDtF2)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlotDtF2))[vecSlotDtF2 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlotDtF2))[vecSlotDtF2 + j] = 0;
//     }

//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.pointByOffset(slotDtF3)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlotDtF3))[vecSlotDtF3 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlotDtF3))[vecSlotDtF3 + j] = 0;
//     }


//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.pointByOffset(slotDt2F1)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlotDt2F1))[vecSlotDt2F1 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlotDt2F1))[vecSlotDt2F1 + j] = 0;
//     }

//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.pointByOffset(slotDt2F2)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlotDt2F2))[vecSlotDt2F2 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlotDt2F2))[vecSlotDt2F2 + j] = 0;
//     }

//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.pointByOffset(slotDt2F3)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlotDt2F3))[vecSlotDt2F3 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlotDt2F3))[vecSlotDt2F3 + j] = 0;
//     }

//   );
}
#endif



bool TripletCalcAngularDt2F::findStage()
{
//   MSG_DEBUG("TripletCalcAngularDt2F::findStage", "START");
  if(m_stage == -1)
  {
  // in the following loops we check, whether there exists a
  // ValCalculator for the force-input symbol. If yes, we check the 
  // stage of it and try to set the stage of this ParticleCache 
  // consistently to it,
    
  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for setting the stage to '0' when there is no Calculator at all 
  // (but probably Integrators or s.th. like that)
    bool nothing = true;     
   
  // first, loop over ColourPairs
    FOR_EACH_COLOUR_PAIR
      (
       M_MANAGER, 
       assert(m_phaseUser == 2 || m_phaseUser == 1);
       vector<ValCalculator*>* vCs;

       // loop over non-bonded vCs
       vCs = &(cp->valCalculatorsFlat());
       // loop over ValCalculators
       for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	   (vCIt != vCs->end() && !tooEarly); ++vCIt)
	 {	   	   
	   list<string> symbols = (*vCIt)->mySymbolNames();
	   
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       if((*symIt) == m_FinName)
		 {
		   nothing = false;
		   int stage = (*vCIt)->stage();
		   if(stage == -1) 
		     {
		       MSG_DEBUG("TripletCalcAngularDt2F::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
		       tooEarly = true;
		       m_stage = -1;
		     }
		   else
		     {
		       if(stage >= m_stage) m_stage = stage+1;
		     }
		 }
	     }
	 } // end loop over non-bonded vCs


       // loop over bonded vCs
       vCs = &(cp->bondedValCalculatorsFlat());
       // loop over ValCalculators
       for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	   (vCIt != vCs->end() && !tooEarly); ++vCIt)
	 {
	   
	   
	   list<string> symbols = (*vCIt)->mySymbolNames();
	   //           MSG_DEBUG("TripletCalcAngularDt2F::findStage", "symbols of VC: ");
	   
	   /*          for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		       cout << *symIt << endl;*/
	   
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       if((*symIt) == m_FinName)
		 {
		   nothing = false;
		   int stage = (*vCIt)->stage();
		   if(stage == -1) 
		     {
		       MSG_DEBUG("TripletCalcAngularDt2F::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
		       tooEarly = true;
		       m_stage = -1;
		     }
		   else
		     {
		       if(stage >= m_stage) m_stage = stage+1;
		     }
		 }
	     }
	 } // end loop over bonded vCs

        // can we stop FOR_EACH_COLOUR_PAIR now?
       if(tooEarly)
	 {
	   __cp = __end;
	   // important because there still comes the ++__cp from the loop
	   --__cp;
	 }
       );
        
    if(!tooEarly)
        {
        // we have to search now in the ParticleCaches
        // loop over stages
          vector<ParticleCache*>* pCs;
          pCs = &(Particle::s_cached_flat_properties);
          FOR_EACH
	    (
	     vector<ParticleCache*>,
	     (*pCs),
	     list<string> symbols = (*__iFE)->mySymbolNames();
	     for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	       {
		 if ((*symIt) == m_FinName)
		   {
		     nothing = false;
		     int stage = (*__iFE)->stage();
		     if(stage == -1) 
		       {
			 MSG_DEBUG("TripletCalcAngularDt2F::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
			 tooEarly = true;
			 m_stage = -1;
		       }
		     else
		       {
			 if(stage >= m_stage) m_stage = stage+1;
		       }
		   }
	       }
	     // may we abort the loop over the ParticleCaches?
	     if(tooEarly)
	       {
		 __iFE = __end;
		 // important because there still comes the ++__iFE from the loop
		 --__iFE; 
	       }
	     );
        } // end of if(!tooEarly) (for loop over ParticleCaches)

    if(!tooEarly)
        {
        // we have to search now in the TripletCalculators
        // loop over stages
          vector<TripletCalculator*>* tCs;
          tCs = M_PHASE->bondedTripletCalculatorsFlat();
          FOR_EACH
	    (
	     vector<TripletCalculator*>,
	     (*tCs),
	     list<string> symbols = (*__iFE)->mySymbolNames();
	     for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	       {
		 
		 if ((*symIt) == m_FinName)
		   {
		     nothing = false;
		     int stage = (*__iFE)->stage();
		     if(stage == -1) 
		       {
			 MSG_DEBUG("TripletCalcAngularDt2F::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
			 tooEarly = true;
			 m_stage = -1;
		       }
		     else
		       {
			 if(stage >= m_stage) m_stage = stage+1;
		       }
		   }
	       }
	     // may we abort the loop over the ParticleCaches?
	     if(tooEarly)
	       {
		 __iFE = __end;
		 // important because there still comes the ++__iFE from the loop
		 --__iFE; 
	       }
	     );
        } // end of if(!tooEarly) (for loop over TripletCalculators)

    if(!tooEarly)
        {
        // we have to search now in the QuintetCalculators
        // loop over stages
          vector<QuintetCalculator*>* tCs;
          tCs = M_PHASE->bondedQuintetCalculatorsFlat();
          FOR_EACH
	    (
	     vector<QuintetCalculator*>,
	     (*tCs),
	     list<string> symbols = (*__iFE)->mySymbolNames();
	     for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	       {
		 
		 if ((*symIt) == m_FinName)
		   {
		     nothing = false;
		     int stage = (*__iFE)->stage();
		     if(stage == -1) 
		       {
			 MSG_DEBUG("TripletCalcAngularDt2F::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
			 tooEarly = true;
			 m_stage = -1;
		       }
		     else
		       {
			 if(stage >= m_stage) m_stage = stage+1;
		       }
		   }
	       }
	     // may we abort the loop over the ParticleCaches?
	     if(tooEarly)
	       {
		 __iFE = __end;
		 // important because there still comes the ++__iFE from the loop
		 --__iFE; 
	       }
	     );
        } // end of if(!tooEarly) (for loop over QuintetCalculators)


        if(tooEarly)
          return false;
        if(m_stage == -1)
        {
          if(nothing)
          {
            m_stage = 0;
            MSG_DEBUG("TripletCalcAngularDt2F::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
            return true; 
          }
          else return false;
        }
        else 
        {
          MSG_DEBUG("TripletCalcAngularDt2F::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
          return true;
        }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("TripletCalcAngularDt2F::findStage", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}

bool TripletCalcAngularDt2F::findStage_0()
{
//   MSG_DEBUG("TripletCalcAngularDt2F::findStage_0", "START");
  if(m_stage == -1)
  {
  // in the following loops we check, whether there exists a
  // ValCalculator for the tensor symbol. If yes, we check the 
  // stage of it and try to set the stage of this ParticleCache 
  // consistently to it,
    
  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for setting the stage to '0' when there is no Calculator at all 
  // (but probably Integrators or s.th. like that)
    bool nothing = true;     
   
  // first, loop over ColourPairs
    FOR_EACH_COLOUR_PAIR
      (
       M_MANAGER, 
       assert(m_phaseUser == 0 || m_phaseUser == 2);
       vector<ValCalculator*>* vCs;
       
       // loop over non-bonded vCs
       vCs = &(cp->valCalculatorsFlat_0());
       // loop over ValCalculators
       for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	   (vCIt != vCs->end() && !tooEarly); ++vCIt)
	 {
	   
	   
	   list<string> symbols = (*vCIt)->mySymbolNames();
	   //           MSG_DEBUG("TripletCalcAngularDt2F::findStage_0", "symbols of VC: ");
	   
	   /*          for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		       cout << *symIt << endl;*/
	   
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       if((*symIt) == m_FinName)
		 {
		   nothing = false;
		   int stage = (*vCIt)->stage();
		   if(stage == -1) 
		     {
		       MSG_DEBUG("TripletCalcAngularDt2F::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
		       tooEarly = true;
		       m_stage = -1;
		     }
		   else
		     {
		       if(stage >= m_stage) m_stage = stage+1;
		     }
		 }
	     }
	 } // end loop over non-bonded vCs


       // loop over bonded vCs
       vCs = &(cp->bondedValCalculatorsFlat_0());
       // loop over ValCalculators
       for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	   (vCIt != vCs->end() && !tooEarly); ++vCIt)
	 {
	   
	   
	   list<string> symbols = (*vCIt)->mySymbolNames();
	   //           MSG_DEBUG("TripletCalcAngularDt2F::findStage_0", "symbols of VC: ");
	   
	   /*          for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		       cout << *symIt << endl;*/
	   
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       if((*symIt) == m_FinName)
		 {
		   nothing = false;
		   int stage = (*vCIt)->stage();
		   if(stage == -1) 
		     {
		       MSG_DEBUG("TripletCalcAngularDt2F::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
		       tooEarly = true;
		       m_stage = -1;
		     }
		   else
		     {
		       if(stage >= m_stage) m_stage = stage+1;
		     }
		 }
	     }
	 } // end loop over bonded vCs



       // can we stop FOR_EACH_COLOUR_PAIR now?
    if(tooEarly)
    {
      __cp = __end;
          // important because there still comes the ++__cp from the loop
      --__cp;
    }
        );
        
        if(!tooEarly)
        {
        // we have to search now in the ParticleCaches
        // loop over stages
          vector<ParticleCache*>* pCs;
//           if(m_phaseUser == 1) pCs = &(Particle::s_cached_flat_properties);
          /*if(m_phaseUser == 0)*/ pCs = &(Particle::s_cached_flat_properties_0);
          FOR_EACH
	    (
	     vector<ParticleCache*>,
	     (*pCs),
	     list<string> symbols = (*__iFE)->mySymbolNames();
	     for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	       {
		 if ((*symIt) == m_FinName)
		   {
		     nothing = false;
		     int stage = (*__iFE)->stage();
		     if(stage == -1) 
		       {
			 MSG_DEBUG("TripletCalcAngularDt2F::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
			 tooEarly = true;
			 m_stage = -1;
		       }
		     else
		       {
			 if(stage >= m_stage) m_stage = stage+1;
		       }
		   }
	       }
	     // may we abort the loop over the ParticleCaches?
	     if(tooEarly)
	       {
		 __iFE = __end;
		 // important because there still comes the ++__iFE from the loop
		 --__iFE; 
	       }
	     );
        } // end of if(!tooEarly) (for loop over ParticleCaches)

	if(!tooEarly)
	  {
	    // we have to search now in the TripletCalculators
	    // loop over stages
          vector<TripletCalculator*>* tCs;
          tCs = M_PHASE->bondedTripletCalculatorsFlat_0();
          FOR_EACH
	    (
	     vector<TripletCalculator*>,
	     (*tCs),
	     list<string> symbols = (*__iFE)->mySymbolNames();
	     for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	       {
		 
		 if ((*symIt) == m_FinName)
		   {
		     nothing = false;
		     int stage = (*__iFE)->stage();
		     if(stage == -1) 
		       {
			 MSG_DEBUG("TripletCalcAngularDt2F::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
			 tooEarly = true;
			 m_stage = -1;
		       }
		     else
		       {
			 if(stage >= m_stage) m_stage = stage+1;
		       }
		   }
	       }
	     // may we abort the loop over the ParticleCaches?
	     if(tooEarly)
	       {
		 __iFE = __end;
		 // important because there still comes the ++__iFE from the loop
		 --__iFE; 
	       }
	     );
	  } // end of if(!tooEarly) (for loop over TripletCalculators)
	
	if(!tooEarly)
	  {
	    // we have to search now in the QuintetCalculators
	    // loop over stages
          vector<QuintetCalculator*>* tCs;
          tCs = M_PHASE->bondedQuintetCalculatorsFlat_0();
          FOR_EACH
	    (
	     vector<QuintetCalculator*>,
	     (*tCs),
	     list<string> symbols = (*__iFE)->mySymbolNames();
	     for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	       {
		 
		 if ((*symIt) == m_FinName)
		   {
		     nothing = false;
		     int stage = (*__iFE)->stage();
		     if(stage == -1) 
		       {
			 MSG_DEBUG("TripletCalcAngularDt2F::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
			 tooEarly = true;
			 m_stage = -1;
		       }
		     else
		       {
			 if(stage >= m_stage) m_stage = stage+1;
		       }
		   }
	       }
	     // may we abort the loop over the ParticleCaches?
	     if(tooEarly)
	       {
		 __iFE = __end;
		 // important because there still comes the ++__iFE from the loop
		 --__iFE; 
	       }
	     );
	  } // end of if(!tooEarly) (for loop over QuintetCalculators)
	
        if(tooEarly)
          return false;
        if(m_stage == -1)
        {
          if(nothing)
          {
            m_stage = 0;
            MSG_DEBUG("TripletCalcAngularDt2F::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
            return true; 
          }
          else return false;
        }
        else 
        {
          MSG_DEBUG("TripletCalcAngularDt2F::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
          return true;
        }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("TripletCalcAngularDt2F::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}

