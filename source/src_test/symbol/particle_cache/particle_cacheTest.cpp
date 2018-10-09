/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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


#include "particle_cacheTest.h"

#include "simulation.h"
#include "cppunit_helper.h"

CPPUNIT_TEST_SUITE_REGISTRATION (ParticleCacheTest);

void ParticleCacheTest :: setUp (void)
{
  // This is required since ParticleCache::setup() calls the
  // Controller of Simulation. This starts getting bulky, but the pain
  // is still below threshold...
  m_simulation = new Simulation();
  Controller* cont = new Controller(m_simulation);
  m_simulation -> setController(cont);
  setUpParticleCache();
}

void ParticleCacheTest :: setUpParticleCache (void)
{
	m_particleCache = new FakeParticleCache(m_simulation);
}

void ParticleCacheTest :: tearDown (void)
{
	delete m_simulation -> controller();
  delete m_simulation;
  delete m_particleCache;
}

void ParticleCacheTest :: initTest (void)
{

	testPropClassName();

	testSpeciesAttr();

	CPPUNIT_ASSERT_EQUAL
  	(string("undefined"), m_particleCache -> species());
}

void ParticleCacheTest :: testPropClassName (void)
{
	CPPUNIT_ASSERT_EQUAL
		(string("ParticleCache"), m_particleCache -> returnProperties().className());
}

void ParticleCacheTest :: testSpeciesAttr()
{

	CppunitHelper::testPropAttr(
		string("species"),
		string(
				"Name for the species of the particles, this Symbol is used for. "
				 "If set to \"ALL\", the Symbol will be used for all registered "
				 "species."),
		(Node*) m_particleCache);

}

void ParticleCacheTest :: setupTest (void)
{

  string exceptMsg =
  		"ERROR: Unspecified exception in ParticleCacheTest :: setupTest for "
  		"module " + m_particleCache->className() + ". Please check!";

  // We catch the first gError and do nothing else, since we do not yet
  // (2018-08-22) test xml-input.
  // FIXME: could we do that with some FakeClass?
  try {
    // better use a copy such that setup() is not run twice
  	ParticleCache* copy =
      (ParticleCache*) (m_particleCache -> returnCopy());
    copy -> setup();
  } catch(gError& err) {
    exceptMsg = err.message(); 
  }

  // Hence, currently (2018-08-22) we expect exceptions (thrown by
  // ParticleCache or parents) due to unset xml-input
  CPPUNIT_ASSERT_THROW_MESSAGE(exceptMsg, m_particleCache -> setup(), gError);
}

size_t colour = 0;

// FIXME: check which design flaw caused this long test function. To be tested
// function too complex?
void ParticleCacheTest :: checkOutputSymbolExistenceTest (void)
{
	FakeParticleCache* fakeParticleCache = (FakeParticleCache*) m_particleCache;

	Particle::s_tag_format.resize(2);

  string exceptMsg;

  string symbolName = m_particleCache -> mySymbolName();

  // Let's test for more than one colour
	for (size_t colour = 0; colour < 2; ++colour) {

		Particle::s_tag_format[colour].setFormatAndAlloc(new DataFormat());

		// Case 1. test function while attribute does NOT exist

		// Case 1.1 m_overwrite = true:

		fakeParticleCache -> setOverwrite(true);

		// set default exception message
	  exceptMsg =
	  		"ERROR: Unspecified exception in ParticleCacheTest :: "
	  		"checkOutputSymbolExistenceTest for module "
	  		+ m_particleCache->className() + ". Please check!";

	  // Get the exception message, if thrown
	  try {
	  	m_particleCache -> checkOutputSymbolExistence(colour);
	  } catch(gError& err) {
	    exceptMsg = err.message();
	  }

	  // If correct, exception should be thrown
	  CPPUNIT_ASSERT_THROW_MESSAGE(exceptMsg,
	  		m_particleCache -> checkOutputSymbolExistence(colour), gError);

	  string referenceString =
	  		//string("test");
	  		string("[ParticleCache::setup for module ")
	  	  		+ m_particleCache -> className()
						+ string("] ERROR: You have chosen 'overwrite = "
	  				"\"yes\"', but symbol '") + symbolName +
	  	  		 string("' does not yet exist. Aborting.");

	  // test exception message
	  CPPUNIT_ASSERT_EQUAL(referenceString, exceptMsg);

		// Case 1.2 m_overwrite = false:

		fakeParticleCache -> setOverwrite(false);

		// Let's add a symbol before the one that should be added next
		Particle::s_tag_format[colour].addAttribute
		    ("firstDummySymbol", DataFormat::POINT, false, "firstDummySymbol");

		// type of symbol that should be added in next step
		fakeParticleCache -> setDatatype(DataFormat::DOUBLE);

		// now this should add a symbol
		m_particleCache -> checkOutputSymbolExistence(colour);

		// check offset (should be location behind the already added data)
		CPPUNIT_ASSERT_EQUAL(
				sizeof(point_t), // since we added a POINT with offsets 0, 8, 16 before
				Particle::s_tag_format[colour].attrByName(m_particleCache
						-> mySymbolName()).offset
		);

		// Case 2. test function while attribute DOES exist
		// (Yes it does now, due to the previous test)

		// Case 2.1 m_overwrite = true

		fakeParticleCache -> setOverwrite(true);

		// Case 2.1.1 correct m_datatype: check offset

		m_particleCache -> checkOutputSymbolExistence(colour);

		// We reached this line. That's test enough for this case.
		CPPUNIT_ASSERT_EQUAL(true, true);

		// Case 2.1.2 wrong m_datatype: check that exception thrown

		// set wrong datatype
		fakeParticleCache -> setDatatype(DataFormat::POINT);

		// set default exception message
	  exceptMsg =
	  		"ERROR: Unspecified or no exception in ParticleCacheTest :: "
	  		"checkOutputSymbolExistenceTest for module "
	  		+ m_particleCache -> className() + ". Colour = " + ObjToString(colour)
				+ ". Please check!";

	  // Get the exception message, if thrown
	  try {
	  	m_particleCache -> checkOutputSymbolExistence(colour);
	  } catch(gError& err) {

	    exceptMsg = err.message();
	  }

	  // If correct, exception should be thrown
	  CPPUNIT_ASSERT_THROW_MESSAGE(exceptMsg,
	  		m_particleCache -> checkOutputSymbolExistence(colour), gError);

	  DataFormat::attribute_t tempAttr
	  	= Particle::s_tag_format[colour].attrByName(symbolName);

	  // test exception message
	  referenceString =
	  		string("[ParticleCache::setup for module ")
				+ m_particleCache -> className()
				+ string("] ERROR: Symbol '" + symbolName)
				+ "' already exists, but with different datatype '"
				+ tempAttr.datatypeAsString()
				+ string("' to be used due to ")
				+ "your choice of 'overwrite = \"yes\"', instead of your "
				+ string("desired datatype '")
				+ DataFormat::attribute_t::datatypeAsString
				(m_particleCache -> symbolType()) + "'. Aborting.";

	  CPPUNIT_ASSERT_EQUAL(referenceString, exceptMsg);

		// Case 2.2 m_overwrite = false: check that exception thrown

	  fakeParticleCache -> setOverwrite(false);

		// set default exception message
	  exceptMsg =
	  		"ERROR: Unspecified exception in ParticleCacheTest :: "
	  		"checkOutputSymbolExistenceTest for module "
	  		+ m_particleCache -> className() + ". Please check!";

	  // Get the exception message, if thrown
	  try {
	  	m_particleCache -> checkOutputSymbolExistence(colour);
	  } catch(gError& err) {
	    exceptMsg = err.message();
	  }

	  // If correct, exception should be thrown
	  CPPUNIT_ASSERT_THROW_MESSAGE(exceptMsg,
	  		m_particleCache -> checkOutputSymbolExistence(colour), gError);

	  // test exception message
	  referenceString =
	  		string("[ParticleCache::setup for module ")
	  		+ m_particleCache -> className() + "] ERROR: Symbol '"
				+ symbolName + "' was already created by other module "
				"for colour " + ObjToString(colour) + ", and you "
				"have chosen overwrite = 'no'.";

	  CPPUNIT_ASSERT_EQUAL(referenceString, exceptMsg);
	}
}

void ParticleCacheTest :: cleanSymbolTest (void)
{
	FakeParticleCache* fakeParticleCache = (FakeParticleCache*) m_particleCache;

	string name = string("{hello}");
	fakeParticleCache -> cleanSymbolPublic(name);
	CPPUNIT_ASSERT_EQUAL(string("hello"), name);

	name = string("[hello]");
	fakeParticleCache -> cleanSymbolPublic(name);
	CPPUNIT_ASSERT_EQUAL(string("hello"), name);

	name = string("hello");
	fakeParticleCache -> cleanSymbolPublic(name);
	CPPUNIT_ASSERT_EQUAL(string("hello"), name);

}

