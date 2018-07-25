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



#include "help.h"

/* --- HelpNode --- */

HelpNode::HelpNode(): m_parent(NULL), m_text("ROOT")
{
}


HelpNode::HelpNode(HelpNode *parent, string text): m_parent(parent), m_text(text)
{
    if (parent)
        parent->addChild(this);
}


HelpNode::~HelpNode()
{
    for (list<HelpNode*>::iterator i = m_children.begin(); i != m_children.end(); i++)
        delete *i;
}



void HelpNode::addChild(HelpNode *child)
{
    m_children.push_back(child);
}



/* --- Functions --- */

string block(string str, size_t indent, size_t width) {
    string s;
    size_t i, size = indent;

    for (i = str.find_first_of(" \n"); str.size() > 0; i = str.find_first_of(" \n")) {
        bool newline = (str[i] == '\n');

        if (i == string::npos)
            i = str.size()-1;

        string chunk(str, 0, i+1);
        str.erase(0, i+1);

        if (size+i+1 >= width) {
            s += "\n" + string(indent, ' ');
            size = indent;
        }

        s += chunk;
        size += i+1;

        if (newline) {
            s += string(indent, ' ');
            size = indent;
        }
    }

    // commented out on 2013-06-17 because it does not allow to use "_" in the help-text of nodes
//     while ((i = s.find('_')) != string::npos)  s[i] = ' ';

    return s;
}


void __HelpFormatScreen(ostream &s, const HelpNode *root, int width, int depth) {
    static string bullets[] = {"", ">>>", "", "--", "-", "", ""};
    static int indents[] = {0, 0, 0, 2, 6, 10, 14};
    list<HelpNode*>::const_iterator obj;

    //    cout << "depth = " << depth << ", root->text() = " << root->text() << endl;

    if (depth == 1 || depth == 2)
        s << endl << endl;

    s << string(indents[depth], ' ') + bullets[depth] + " "
        + block(root->text(), indents[depth]+bullets[depth].size()+1, width) + " ";

    if (depth == 3)
        s << string(width - indents[depth] - bullets[depth].size() - 6 - root->text().size(), '-');

    s << endl;

    for (obj = root->children().begin(); obj != root->children().end(); obj++) {
        __HelpFormatScreen(s, *obj, width, depth+1);
    }
}


void HelpFormatScreen(ostream &s, const HelpNode &root, int start_depth) {
    __HelpFormatScreen(s, &root, 80, start_depth);
}
