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



#include "pc_wall_textured.h"
#include "phase.h"
#include "simulation.h"
#include "manager_cell.h"
#include <iostream>
#include <sys/types.h>

/* Register this ParticleCreator with the factory. */
const ParticleCreator_Register<ParticleCreatorWallTextured>
		particle_creator_wall_textured("ParticleCreatorWallTextured");

/* --- */

std::istream& Raster::skipComments(std::istream& in) {
	int c = in.peek();
	while (c == '#'|| c == '\n') {
		if (c == '#')
			in.ignore(INT_MAX, '\n');
		else
			c = in.get();
		c = in.peek();
	}
	return in;
}
Raster::Raster() {
}
Raster::Raster(string fname) {
	data = NULL;
	readFromPGMFile(fname);
	defaultValue = 0;
}
Raster::Raster(string fname, int def) {
	data = NULL;
	readFromPGMFile(fname);
	defaultValue = def;
}
Raster::~Raster() {
	if (data != NULL)
		delete data;
}
void Raster::readFromPGMFile(string fname) {
	std::ifstream f;
	unsigned int maxgray;
	unsigned long wh;
	string s;

	f.exceptions(std::ifstream::eofbit| std::ifstream::failbit
			| std::ifstream::badbit);

	try {
		f.open(fname.c_str(), ios::in);

		// Check header
		skipComments(f >> skipws) >> s;
		if(s == "P2" || s == "P5") {
			skipComments(f >> skipws) >> width;
			skipComments(f >> skipws) >> height;
			skipComments(f >> skipws) >> maxgray;
			wh = width*height;
			data = new u_int16_t[height*width];
		}
		else throw gError("Raster::"+PRETTYFUNC_INFO, "Problem reading "
				"from file "+fname+": wrong file format "+s);
		if(s == "P2") { // plain, ASCII
			for(size_t row = 0; row < height; row++)
			for(size_t col = 0; col < width; col++)
			skipComments(f >> skipws) >> data[row*width+col];
		}
		else { // raw, binary
			if(maxgray >= 256)
			f.read(reinterpret_cast<char*>(data), wh * sizeof(u_int16_t)/sizeof(char));
			else {
				u_int8_t* cdata = new u_int8_t[height*width];
				f.read(reinterpret_cast<char*>(cdata), wh * sizeof(u_int8_t)/sizeof(char));
				for(size_t row = 0; row < height; row++)
				for(size_t col = 0; col < width; col++) {
					data[row*width+col] = cdata[row*width+col];
				}
				delete cdata;
			}
		}
	}
	catch (std::ifstream::failure e) {
		throw gError("Raster::"+PRETTYFUNC_INFO, "Problem reading "
				"from file "+fname+": " + e.what());
	}
}
unsigned int Raster::operator()(int row, int col) {
	if (row >= height || col >= width || row < 0|| col < 0)
	return defaultValue;
	else
	return data[row*width+col];
}


//---- Constructor/Destructor ----
ParticleCreatorWallTextured::ParticleCreatorWallTextured(Boundary *boundary) :
	ParticleCreatorWall(boundary) {
	// init() should be virtual, but isn't. So it is called already
	// by ParticleCreatorWall(boundary).
	init();
}

ParticleCreatorWallTextured::~ParticleCreatorWallTextured() {
	ParticleCreatorWall::already = false; //shared with superclass
}

//---- Methods ----
inline double ParticleCreatorWallTextured::getNearestRaster(int i, int j, int k) {
	int dmin= INT_MAX;
	double v = 0;
	int c = 1;
	if (textureMin[0]) {
			v = double((*textureMin[0])(j, k));
			dmin = i;
		}
	if (textureMin[1])
		if (j < dmin) {
			v = double((*textureMin[1])(i, k));
			c = 1;
			dmin = j;
		} else if (j == dmin) {
			v = (double((*textureMin[1])(i, k)) + v*c) / (c+1);
			++c;
		}
	if (textureMin[2])
		if (k < dmin) {
			v = double((*textureMin[2])(i, j));
			c = 1;
			dmin = k;
		} else if (k == dmin) {
			v = (double((*textureMin[2])(i, j)) + v*c) / (c+1);
			++c;
		}
	if (textureMax[0])
		if (m_nlattice_points[0]-1-i < dmin) {
			v = double((*textureMax[0])(j, k));
			c = 1;
			dmin = m_nlattice_points[0]-1-i;
		} else if (m_nlattice_points[0]-1-i == dmin) {
			v = (double((*textureMax[0])(j, k)) + v*c) / (c+1);
			++c;
		}
	if (textureMax[1])
		if (m_nlattice_points[1]-1-j < dmin) {
			v = (*textureMax[1])(i, k);
			c = 1;
			dmin = m_nlattice_points[1]-1-j;
		} else if (m_nlattice_points[1]-1-j == dmin) {
			v = (double((*textureMax[1])(i, k)) + v*c) / (c+1);
			++c;
		}
	if (textureMax[2])
		if (m_nlattice_points[2]-1-k < dmin) {
			v = (*textureMax[2])(i, j);
			c = 1;
			dmin = m_nlattice_points[2]-1-k;
		} else if (m_nlattice_points[2]-1-k == dmin) {
			v = (double((*textureMax[2])(i, j)) + v*c) / (c+1);
			++c;
		}
	return v;
}


void ParticleCreatorWallTextured::createParticles() {
	// 	MSG_DEBUG("ParticleCreatorWallTextured::createParticles()", "start");    


	/*!
	 *  calculate which of the cube's planes are nearest, then assign
	 *  the value of that plane's texture. When more than one plane
	 *  is nearest, assign the average.
	 */

	Phase *phase = (Phase*) ((Boundary*) m_parent)->parent();
	point_t size = ((Boundary*) m_parent)->boundingBox().size();
	point_t vel = { { { 0, 0, 0 } } };
	size_t colour = phase->manager()->getColour(m_species);

	point_t boxOffset = phase->boundary()->boundingBox().corner1;
	point_t pos;
	size_t n = 0;
	ptrdiff_t attOffs;

	MSG_DEBUG("ParticleCreatorWallTextured::"+PRETTYFUNC_INFO,
			"Color = " << colour << ", species = " << m_species << ", "
			"interactionCutoff = " << myCutoff << ", "
			"textures = " << m_filenameMin[0] << ", " << m_filenameMin[1] << ", "
			<< m_filenameMin[2] << ", " << m_filenameMax[0] << ", "
			<< m_filenameMax[1] << ", "<< m_filenameMax[2]);

	for(int i=0; i<SPACE_DIMS; i++)
		textureMin[i] = textureMax[i] = NULL;

	// real PGM files
	for (int i = 0; i < SPACE_DIMS; ++i) {
		if (m_filenameMin[i] != "")
			textureMin[i] = new Raster(m_filenameMin[i]);
		if (m_filenameMax[i] != "")
			textureMax[i] = new Raster(m_filenameMax[i]);
	}

	// add an attribute for the texture to the particle data tag format
	if (Particle::s_tag_format[colour].attrExists(m_textureSymbol))
		throw gError("ParticleCreatorWallTextured::"+PRETTYFUNC_INFO,
				"Symbol " + m_textureSymbol + " is already existing "
				"for species '" + phase->manager()->species(colour) +
				"'. Second definition is not allowed.");
	attOffs = Particle::s_tag_format[colour].addAttribute(m_textureSymbol,
			DataFormat::DOUBLE, /*persist.first*/false, m_textureSymbol).offset;

	MSG_DEBUG("ParticleCreatorWallTextured::createParticles()", "boxOffset = " << boxOffset);
	MSG_DEBUG("ParticleCreatorWallTextured::createParticles()", "m_nlattice_points[0] = " << m_nlattice_points[0]);
	MSG_DEBUG("ParticleCreatorWallTextured::createParticles()", "m_nlattice_points[1] = " << m_nlattice_points[1]);
	MSG_DEBUG("ParticleCreatorWallTextured::createParticles()", "m_nlattice_points[2] = " << m_nlattice_points[2]);

	// build a list of functions corresponding to the unknown properties from m_properties.unknown()
	list<FunctionFixed> fList;
	for (map<string, string>::const_iterator i = m_properties.unknown().begin(); i != m_properties.unknown().end(); ++i) {
		// the functions will be deleted automatically at the end of createParticles 
		// when the list fList is deleted
		FunctionFixed func;
		fList.push_back(func);
		FunctionFixed& lastF = fList.back();
		lastF.addVariable("x");
		lastF.addVariable("y");
		lastF.addVariable("z");
		lastF.addVariable("u");
		lastF.addVariable("v");
		lastF.addVariable("w");
		lastF.setOnlyExpr(i->second);
		lastF.compile();
	}

	m_mappixel.compile();

	// create the particles
	for (int i = 0; i < m_nlattice_points[0]; i++) {
		for (int j = 0; j < m_nlattice_points[1]; j++) {
			for (int k = 0; k < m_nlattice_points[2]; k++) {
				
				// calculate position
				pos.x = i*m_distance+m_distance/2;
				pos.y = j*m_distance+m_distance/2;
				pos.z = k*m_distance+m_distance/2;
				pos += boxOffset;

				if (!(((Boundary*) m_parent)->isInside(pos))&&((Boundary*) m_parent)->isInWallRange(
						myCutoff, pos) ) {
					Particle *p;
					if (m_freeze)
						p = phase->addFrozenParticle(pos, vel, 0, colour);

					else
						p = phase->addParticle(pos, vel, 0, colour);

					p->v.x = m_velX(p->r.x, p->r.y, p->r.z, p->v.x, p->v.y,
							p->v.z);
					p->v.y = m_velY(p->r.x, p->r.y, p->r.z, p->v.x, p->v.y,
							p->v.z);
					p->v.z = m_velZ(p->r.x, p->r.y, p->r.z, p->v.x, p->v.y,
							p->v.z);
					p->tag.doubleByOffset(attOffs)
							= m_mappixel(getNearestRaster(i, j, k));

					list<FunctionFixed>::iterator fListIt = fList.begin();
					for (map<string, string>::const_iterator
							i = m_properties.unknown().begin(); i != m_properties.unknown().end(); ++i) {
						assert(fListIt != fList.end());
						// MSG_DEBUG("ParticleCreatorWallTextured::createParticles", "checking for unknown property " << i->first);
						// new style treatment of unknown properties
						double* value;
						DataFormat::datatype_t type;
						size_t offset;
						string name = i->first;
						if (Particle::s_tag_format[colour].attrExists(name)) {

							offset = Particle::s_tag_format[colour].attrByName(name).offset;
							type = Particle::s_tag_format[colour].attrByName(name).datatype;

							if (type == DataFormat::DOUBLE) {
								// value is a reference
								value = &(p->tag.doubleByOffset(offset));
							} else
								throw gError("ParticleCreatorWallTextured::createParticles", "For species '" + m_species + "': Scalar-type definition used for non-scalar property '" + name + "'. For vectors, you must add \"__x\", \"__y\" or \"__z\" and for tensors you must add \"__xx\", \"__xy\", \"__xz\", \"__yx\", \"__yy\", \"__yz\", \"__zx\", \"__zy\" or \"__zz\" to the property name.");
						} else {
							if (name[name.size()-2] == '_'&& name[name.size()-3] == '_') {
								string symbol = name;
								if (Particle::s_tag_format[colour].attrExists(symbol.erase(
										symbol.size()-3, symbol.size()-1)))
									type
											= Particle::s_tag_format[colour].attrByName(symbol).datatype;
								else
									throw gError("ParticleCreatorWallTextured::createParticles", "For species '" + m_species + ": Found a vector extension for attribute '" + name + "' but did not find symbol '" + symbol + "'.");
								if (type != DataFormat::POINT)
									if (type != DataFormat::TENSOR)
										throw gError("ParticleCreatorWallTextured::createParticles", "For species '" + m_species + ": You use a vector-extension in attribute '" + name + "', but '" + symbol + "' is not a vector.");
								offset
										= Particle::s_tag_format[colour].attrByName(symbol).offset;

								switch (name[name.size()-1]) {
								case 'x':
									value = &(p->tag.pointByOffset(offset).x);
									break;
								case 'y':
									value = &(p->tag.pointByOffset(offset).y);
									break;
								case 'z':
									value = &(p->tag.pointByOffset(offset).z);
									break;
								default:
									throw gError("ParticleCreatorWallTextured::createParticles", "For species '" + m_species + ": Unknown vector extension '__" + name[name.size()-1] + "' in attribute '" + name + "'. Possibilities are \"__x\", \"__y\" or \"__z\".");
								}
							} // end if(vector-case)
							else {
								if (name[name.size()-3] == '_'&& name[name.size()-4] == '_') {
									string symbol = name;
									if (Particle::s_tag_format[colour].attrExists(symbol.erase(
											symbol.size()-4, symbol.size()-1)))
										type
												= Particle::s_tag_format[colour].attrByName(symbol).datatype;
									else
										throw gError("ParticleCreatorWallTextured::createParticles", "For species '" + m_species + ": Found a tensor extension for attribute '" + name + "' but did not find symbol '" + symbol + "'.");
									if (type != DataFormat::TENSOR)
										throw gError("ParticleCreatorWallTextured::createParticles", "For species '" + m_species + ": You use a tensor-extension in attribute '" + name + "', but '" + symbol + "' is not a tensor.");
									offset
											= Particle::s_tag_format[colour].attrByName(symbol).offset;

									bool error = false;
									switch (name[name.size()-2]) {
									case 'x':
										switch (name[name.size()-1]) {
										case 'x':
											value = &(p->tag.tensorByOffset(offset).xx);
											break;
										case 'y':
											value = &(p->tag.tensorByOffset(offset).xy);
											break;
										case 'z':
											value = &(p->tag.tensorByOffset(offset).xz);
											break;
										default:
											error = true;
										}
										break;
									case 'y':
										switch (name[name.size()-1]) {
										case 'x':
											value = &(p->tag.tensorByOffset(offset).yx);
											break;
										case 'y':
											value = &(p->tag.tensorByOffset(offset).yy);
											break;
										case 'z':
											value = &(p->tag.tensorByOffset(offset).yz);
											break;
										default:
											error = true;
										}
										break;
									case 'z':
										switch (name[name.size()-1]) {
										case 'x':
											value = &(p->tag.tensorByOffset(offset).zx);
											break;
										case 'y':
											value = &(p->tag.tensorByOffset(offset).zy);
											break;
										case 'z':
											value = &(p->tag.tensorByOffset(offset).zz);
											break;
										default:
											error = true;
										}
										break;
									}
									if (error)
										throw gError("ParticleCreatorWallTextured::createParticles", "For species '" + m_species + ": Unknown tensor extension '__" + name[name.size()-2] + name[name.size()-1] + "' in attribute '" + name + "'. Possibilities are \"__xx\", \"__xy\", \"__xz\", \"__yx\", \"__yy\", \"__yz\", \"__zx\", \"__zy\" or \"__zz\".");

								} // end if(tensor-case)
								else {
									MSG_DEBUG("ParticleCreatorWallTextured::createParticles", "NOT FOUND " << i->first);
									throw gError("ParticleCreatorWallTextured::createParticles",
											"No internal degree of freedom named '" + i->first +
											"' for species '" + m_species + "'. Possibilities are: " + Particle::s_tag_format[colour].toString());
								}
							} // end else of if(vector-case)
						} // end else of if(attrExists)

						*value = (*fListIt)(p->r.x, p->r.y, p->r.z, p->v.x,
								p->v.y, p->v.z);
						++fListIt;
					}
					//           MSG_DEBUG("ParticleCreatorWallTextured::createParticles", n << " particles created.");
					++n;

				}
			}
		}
	}
	MSG_DEBUG("ParticleCreatorWallTextured::createParticles", n << " particles created.");
}

void ParticleCreatorWallTextured::init() {
	// ParticleCreatorWall::init() is called before automatically by the constructor 
	m_properties.setClassName("ParticleCreatorWallTextured");
	m_properties.setDescription
		("Generates particles behind boundary walls, i.e., outside the domain "
		"that a freely moving particle can reach.\n"
		"The user can give a texture to the walls. The textures are "
		"defined in PGM files (see http://netpbm.sourceforge.net/doc/pgm.html)"
		"whose gray values are read in. These gray values can further be modified"
		"by a function given in mapPixel. For each side of the bounding cube, "
		"a different texture can be assigend. The particles inside the cube "
		"will get the value of their projection to the nearest side of the cube.\n"
		"The initial conditions of user-defined symbols can be set by "
		"taking the symbol as an attribute and defining a mathematical "
		"expression for it. For the expression, the same variables are "
		"allowed as, e.g., for the attribute 'u'. For non-scalars you "
		"have to add the following to the attribute name:\n"
		"\"__x\", \"__y\" or \"__z\" for the respective components of a vector\n"
		"\"__xx\", \"__xy\", \"__xz\", \"__yx\", \"__yy\", \"__yz\", \"__zx\", "
		"\"__zy\" and \"__zz\" for the respective components of a tensor.");

	FUNCTIONFIXEDPC(mapPixel, m_mappixel,
			"This expression maps the pixel values 'v' (0..255) to a double which than can be used in expressions. Default for v is 0.")
	;

	for (int i = 0; i < SPACE_DIMS; ++i) {
		m_properties.addProperty(string("nameMin")+char('X'+i)+string("InputFile"), PropertyList::STRING,
				&(m_filenameMin[i]), NULL,
				string("PGM File containing the texture"
					" for the 'lower end' of the ")+char('x'+i)+string("-direction."
					"Default for position outside is 0."));
		m_properties.addProperty(string("nameMax")+char('X'+i)+string("InputFile"), PropertyList::STRING,
				&(m_filenameMax[i]), NULL,
				string("PGM File containing the texture"
					" for the 'upper end' of the ")+char('x'+i)+string("-direction. "
					"Default for position outside is 0."));
	}
	STRINGPC(textureSymbol, m_textureSymbol,
			"Name of the symbol containing the texture of the particle.")
	;
	
	m_textureSymbol = "UNDEF";

	for (int i = 0; i < SPACE_DIMS; ++i) {
		m_filenameMin[i] = "";
		m_filenameMax[i] = "";
	}

	m_mappixel.addVariable("v");
	m_mappixel.setExpression("v");
}

void ParticleCreatorWallTextured::setup() {
	ParticleCreatorWall::setup();
}
