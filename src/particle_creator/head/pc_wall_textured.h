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



#ifndef __PC_WALL_TEXTURED_H_
#define __PC_WALL_TEXTURED_H_

#include "pc_wall.h"
#include "boundary.h"
#include "particle_creator.h"
#include "function_fixed.h"

/*!
 * Class representing a raster which can be read from
 * a PGM file.
 */
class Raster {
public:
	size_t width;
	size_t height;
	u_int16_t* data;
	int defaultValue;
	Raster();
	Raster(string fname);
	Raster(string fname, int def);
	~Raster();
	void readFromPGMFile(string fname);
	unsigned int operator()(int row, int col);
protected:
	std::istream& skipComments(std::istream& in);
};


/* ---- ParticleCreatorWallTextured ---- */

/*!
 * \brief Creates particle texture behind walls from pgm file.
 * 
 * Generates particles behind walls, i.e., outside the domain
 * that a freely moving particle can reach. The particles are
 * made from the same species. The pattern is read from a file
 * in PGM format (http://netpbm.sourceforge.net/doc/pgm.html) 
 * and stored in an attribute; a mapping function can be used
 * to modify the pixel value before.
 * 
 * For each side of the cube, a different pattern can be used
 * 
 */
class ParticleCreatorWallTextured: public ParticleCreatorWall {
protected:
	
	/*
	 * File names with textures. Must be PGM files with
	 * no more than 256 colors (1 byte/pixel for binary
	 * format). 
	 */
	string m_filenameMin[SPACE_DIMS], m_filenameMax[SPACE_DIMS];

	/*
	 * The textures of the different sides of the cube.
	 */
	Raster *textureMin[SPACE_DIMS], *textureMax[SPACE_DIMS];
	/*!
	 * This expression maps the pixel values \a v (0..255) to a double which than can be used in expressions.
	 */
	FunctionFixed m_mappixel;

	/*!
	 * Name of the symbol containing the texture of the particle.
	 */
	string m_textureSymbol;
	
	/*!
	 * Initialisation of the \a PropertyList
	 */
	void init();

	/*!
	 * Setup for the \a ParticleCreatorWall
	 */
	virtual void setup();
	
	/*!
	 * Get the value corresponding to the current position and
	 * the nearest available texture
	 */
	inline double getNearestRaster(int i, int j, int k);


public:
	/*!
	 * Constructor
	 * @param boundary The \a Boundary, this \a ParticleCreator belongs to 
	 */
	ParticleCreatorWallTextured(Boundary *boundary);

	/*!
	 * Destructor
	 */
	virtual ~ParticleCreatorWallTextured();

	/*!
	 * The routine for creating wall \a Particles
	 */
	virtual void createParticles();

};


#endif /*__PC_WALL_TEXTURED_H_*/
