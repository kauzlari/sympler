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



#include "meter_live_geomview.h"
#include "phase.h"
#include "simulation.h"
// WARNING: This is not very portable!
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <arpa/inet.h>
#include <netinet/in.h> //for hton*
#include <string.h>

#include <assert.h>

// default name of executable
#define GEOMVIEWEXEC "geomview -c -"

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<MeterLiveGeomview>
		meter_live_geomview("MeterLiveGeomview");

MeterLiveGeomview::MeterLiveGeomview(Simulation *simulation) :
	Meter(simulation), m_geomexec(GEOMVIEWEXEC) {
	init();
}

MeterLiveGeomview::~MeterLiveGeomview() {
	if (m_pipe)
		pclose(m_pipe);
}

//---- Methods ----

void MeterLiveGeomview::init() {

	m_properties.setClassName("MeterLiveGeomview");

	m_properties.setDescription("Graphical view of the simulation at runtime using geomview. The colour of the particles represents their species.");

	STRINGPC(geomexec, m_geomexec, "Name of the geomview executable.")
	;
	m_pipe = 0;
	m_geomexec = "geomview";
}

void MeterLiveGeomview::setup() {

	Meter::setup();

	if ((m_pipe = popen(m_geomexec.c_str(), "w")) == NULL)
		throw gError("MeterLiveGeomview::setup",
				string("Could not open geomview: ") + strerror(errno));
	if ((fprintf(m_pipe, "(geometry sympler { : particles })\n")) < 0)
		throw gError("MeterLiveGeomview::setup",
				string("Could not write to geomview: ") + strerror(errno));
	fflush(m_pipe);
}

 /*!
  * Convert a float to big endian
  * @param Float to convert
  */

#ifdef IS_BIG_ENDIAN
float MeterLiveGeomview::htonf(float f) {
	return f;
}
#else
float MeterLiveGeomview::htonf(float f) {
	char f2[4];
	char c[4];
	memcpy(&f2, &f, sizeof(float));
	c[0] = f2[3];
	c[1] = f2[2];
	c[2] = f2[1];
	c[3] = f2[0];
	return *((float*)&c);
}
#endif

void MeterLiveGeomview::measureNow(const double &time) {
	Phase *phase = ((Simulation*) m_parent)->phase();
	point_t c1 = phase->boundary()->boundingBox()[0];
	point_t c2 = phase->boundary()->boundingBox()[1];
	size_t colors = phase->nColours();
	bool freecolors[8] = {false,false,false,false,false,false,false,false};
	bool frozencolors[8] = {false,false,false,false,false,false,false,false};
	size_t filecolors = 0;
	size_t npart = phase->returnNofPart();
	size_t nfrozenpart = phase->returnNofFrozenP();
	int32_t h[3];
	size_t p;
	point_t r;

	struct {
		int32_t hd[3]; // NPolylines  NVertices  NColors
		int16_t nv1[2]; // No. of vertices per polyline 
		int16_t nc1[2]; // No. of colors per polyline
		float v1[3]; // Vertices
		float v2[3]; // Vertices
		float c1[4]; // Colors
	} h1;
	float colorlist[] = { 
			htonf(0.f), htonf(0.f), htonf(1.f), htonf(1.f),
			htonf(0.f), htonf(1.f), htonf(0.f), htonf(1.f), 
			htonf(1.f), htonf(0.f), htonf(0.f), htonf(1.f),
			htonf(1.f), htonf(1.f), htonf(0.f), htonf(1.f),
			htonf(1.f), htonf(0.f), htonf(1.f),	htonf(1.f),
			htonf(0.f), htonf(1.f), htonf(1.f), htonf(1.f),
			htonf(1.f), htonf(1.f), htonf(1.f), htonf(1.f),
			htonf(0.f),	htonf(0.f), htonf(0.f), htonf(1.f) };
	float frozencolorlist[] = { 
			htonf(0.f), htonf(0.f), htonf(0.5f), htonf(1.f),
			htonf(0.f), htonf(0.5f), htonf(0.f), htonf(1.f), 
			htonf(0.5f), htonf(0.f), htonf(0.f), htonf(1.f),
			htonf(0.5f), htonf(0.5f), htonf(0.f), htonf(1.f),
			htonf(0.5f), htonf(0.f), htonf(0.5f),	htonf(1.f),
			htonf(0.f), htonf(0.5f), htonf(0.5f), htonf(1.f),
			htonf(0.5f), htonf(0.5f), htonf(0.5f), htonf(1.f),
			htonf(0.f),	htonf(0.f), htonf(0.f), htonf(1.f) };
	int16_t *i2s;
	int16_t *j2s;
	float *fs;

	// create two points for the bounding box
	h1.hd[0] = htonl(2);
	h1.hd[1] = htonl(2);
	h1.hd[2] = htonl(1);
	h1.nv1[0] = htons(1);
	h1.nv1[1] = htons(1);
	h1.nc1[0] = htons(1);
	h1.nc1[1] = htons(0);
	h1.v1[0] = htonf(c1[0]);
	h1.v1[1] = htonf(c1[1]);
	h1.v1[2] = htonf(c1[2]);
	h1.v2[0] = htonf(c2[0]);
	h1.v2[1] = htonf(c2[1]);
	h1.v2[2] = htonf(c2[2]);
	h1.c1[0] = htonf(0.0f);
	h1.c1[1] = htonf(0.0f);
	h1.c1[2] = htonf(0.0f);
	h1.c1[3] = htonf(0.0f);

#define METERLIVEGEOMVIEW_FWRITE(what,size) if(fwrite(what,1,size,m_pipe)<0) throw gError("MeterLiveGeomview::setup",\
				string("Could not write to geomview: ") + strerror(errno))	; while(0)
#define METERLIVEGEOMVIEW_FPRINTF(what, ...) if(fprintf(m_pipe,what,##__VA_ARGS__)<0) throw gError("MeterLiveGeomview::setup",\
				string("Could not write to geomview: ") + strerror(errno))	; while(0)

	METERLIVEGEOMVIEW_FPRINTF("(read geometry { define particles appearance { linewidth 2 } LIST { VECT BINARY\n")
	;
	METERLIVEGEOMVIEW_FWRITE(&h1,3*4+4*2+10*4)
	;
	METERLIVEGEOMVIEW_FPRINTF("\nVECT BINARY\n")
	;

	if ((i2s = (int16_t*)malloc((npart+nfrozenpart)*2)) == NULL)
		throw gError("MeterLiveGeomview::setup",
				string("Could not allocate memory: ") + strerror(errno));
	if ((j2s = (int16_t*)calloc((npart+nfrozenpart+1),2)) == NULL)
		throw gError("MeterLiveGeomview::setup",
				string("Could not allocate memory: ") + strerror(errno));
	if ((fs = (float*)malloc((npart+nfrozenpart)*3*4)) == NULL)
		throw gError("MeterLiveGeomview::setup",
				string("Could not allocate memory: ") + strerror(errno));

	p = 0;
	for (size_t c = 0; c < colors; c++) {
		if (c<8) {
			if (phase->returnNofPartC(c) > 0)
			{
				freecolors[c] = true;
				filecolors++;
				j2s[p] = htons(1);
			} else
				freecolors[c] = false;
		}
		SL_FOR_EACH(Particle, phase->particles(c),
				r = __iSLFE->r;
					fs[3*p] = htonf(r[0]); fs[3*p+1] = htonf(r[1]); fs[3*p+2] = htonf(r[2]);
					p++;
			)
			;
		if (c<8) {
			if (phase->returnNofFrozenPC(c) > 0)
			{
				frozencolors[c] = true;
				filecolors++;
				j2s[p] = htons(1);
			} else
				frozencolors[c] = false;
		}
		SL_FOR_EACH(Particle, phase->frozenParticles(c),
				r = __iSLFE->r;
				fs[3*p] = htonf(r[0]); fs[3*p+1] = htonf(r[1]); fs[3*p+2] = htonf(r[2]);
				p++;
		)
		;
	}
	assert(npart+nfrozenpart == p);

	h[0] = h[1] = htonl((int32_t)npart+nfrozenpart);
	h[2] = htonl((int32_t)(filecolors));

	METERLIVEGEOMVIEW_FWRITE(&h,3*4)
	;
	// NPolylines  NVertices  NColors

	for (int i = 0; i < npart+nfrozenpart; i++)
		i2s[i] = htons(1);
	METERLIVEGEOMVIEW_FWRITE(i2s,2 * (npart+nfrozenpart))
	;
	// No. of vertices per polyline 

	METERLIVEGEOMVIEW_FWRITE(j2s,2 * (npart+nfrozenpart))
	;
	// No. of colors per polyline
	free(i2s);
	free(j2s);

	METERLIVEGEOMVIEW_FWRITE(fs,3 * 4 * (npart+nfrozenpart))
	;
	// vectors
	free(fs);

	for (int i = 0; i < 8; i++)
	{
		if(freecolors[i])
			METERLIVEGEOMVIEW_FWRITE(&colorlist[4*i],4 * 4);
		if(frozencolors[i])
			METERLIVEGEOMVIEW_FWRITE(&frozencolorlist[4*i],4 * 4);
	}	// colors 	
	METERLIVEGEOMVIEW_FPRINTF("}})\n")
	;
	fflush(m_pipe);
	// -- non-binary --
	//	METERLIVEGEOMVIEW_FPRINTF("\nVECT\n");
	//	h[0] = h[1] = npart+nfrozenpart;
	//	h[2] = (colors>8?8:colors);
	//
	//	METERLIVEGEOMVIEW_FPRINTF("%i %i %i\n",h[0],h[1],h[2]); // NPolylines  NVertices  NColors
	//
	//	for(int i = 0; i < npart+nfrozenpart; i++)
	//		METERLIVEGEOMVIEW_FPRINTF("1 ");
	//	METERLIVEGEOMVIEW_FPRINTF("\n");
	//
	//	for(int i = 0; i < 2; i++)
	//		METERLIVEGEOMVIEW_FPRINTF("1 ");
	//	for(int i = 2; i < npart+nfrozenpart; i++)
	//		METERLIVEGEOMVIEW_FPRINTF("0 ");
	//
	//	p = 0;
	//	for (size_t c = 0; c < colors; c++) {
	//	    SL_FOR_EACH(Particle, phase->particles(c), 
	//	    	METERLIVEGEOMVIEW_FPRINTF("%f %f %f\n",__iSLFE->r[0],__iSLFE->r[1],__iSLFE->r[2]);
	//	    	p++;
	//	    ); 
	//	    SL_FOR_EACH(Particle, phase->frozenParticles(c), 
	//		    	METERLIVEGEOMVIEW_FPRINTF("%f %f %f\n",__iSLFE->r[0],__iSLFE->r[1],__iSLFE->r[2]);
	//		    	p++;
	//	    );
	//	}
	//	
	//	METERLIVEGEOMVIEW_FPRINTF("1 0 0 1 0 1 0 1"); // colors 	
	//	cout << "**** p" << p << " colors " << colors << " npart " << npart << " nfrozenpart " << nfrozenpart <<endl;
	//	METERLIVEGEOMVIEW_FPRINTF("}})\n");
	//	fflush(m_pipe);
}
