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



#include <fstream>

#include "sat_parser.h"

#include "sat_unknown.h"

#define END_OF_DATA_STR "End-of-ACIS-data"

/*--- SATParser ---*/

SATParser::SATParser()
{
}


SATParser::~SATParser()
{
    for (vector<SATEntity*>::iterator i = m_entities.begin(); i != m_entities.end(); i++) {
        delete *i;
    }
}


void SATParser::read(string filename)
{
    ifstream s;

    s.open(filename.c_str());
    read(&s);
    s.close();
}


void SATParser::read(istream *s)
{
    string name;

    m_stream = s;

    m_header.version = readInt();
    m_header.n_data_records = readInt();
    m_header.n_entities = readInt();
    m_header.has_history = readBool();
    finishLine();

    m_header.product_id = readString();
    m_header.acis_version = readString();
    m_header.date = readString();
    finishLine();

    m_header.unit = readDouble();
    m_header.resabs = readDouble();
    m_header.resnor = readDouble();
    MSG_DEBUG("SATParser::read", "resnor = " << m_header.resnor);
    finishLine();

    //    cout << readString() << endl;

    name = readWord();
    while (!m_stream->eof() && name != END_OF_DATA_STR) {
        SATEntity *e;

        if (SATEntity_Factory::exists(name)) {
            e = SATEntity_Factory::byName(name).instantiate();
        } else {
            MSG_DEBUG("SATParser::read", "Unknown entity: '" << name << "'");

            e = new SATUnknown();
        }

        e->read(this);

        m_entities.push_back(e);

        finishLine();

        name = readWord();
    }

    fixPointers();
}


void SATParser::fixPointers()
{
    for (vector<SATEntity*>::iterator i = m_entities.begin(); i != m_entities.end(); i++) {
        (*i)->fixPointers(this);
    }
}


/*--- Parser functions ---*/

#define MAX_STR_LEN 80

string SATParser::readString()
{
    int length;
    char s[MAX_STR_LEN+1];

    (*m_stream) >> skipws >> length >> skipws;

    MSG_DEBUG("SATParser::readString", "length = " << length);

    if (length > MAX_STR_LEN)
        throw SATParseError("SATParser::readString", "Maximum string length exceeded.");

    m_stream->get();
    (*m_stream).get(s, length+1);
    s[length] = 0;

    MSG_DEBUG("SATParser::readString", "--> |" << s << "|");

    return string(s);
}

string SATParser::readWord()
{
    string s;

    (*m_stream) >> skipws >> s;

    return s;
}

int SATParser::readInt()
{
    int i;

    (*m_stream) >> skipws >> i;

    return i;
}

SATPointer SATParser::readPointer()
{
    SATPointer i;

    (*m_stream) >> skipws;

    m_stream->get();
    if (m_stream->get() == '$') {
        (*m_stream) >> i.index;
    } else {
        throw SATParseError
            ("SATParser::readPointer", "This is no pointer.");
    }

    return i;
}

double SATParser::readDouble()
{
    double d;

    (*m_stream) >> skipws >> d;

    return d;
}

bool SATParser::readBool()
{
    int i;

    (*m_stream) >> skipws >> i;

    return i&1;
}

void SATParser::finishLine()
{
    while (m_stream->get() != '\n' && !m_stream->eof());
}
