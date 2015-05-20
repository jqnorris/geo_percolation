// Collection of useful tools

#include "Tools.h"
#include "Abstract_Classes.h"
#include "Forward_Declarations.h"
#include "Site.cpp"
#include "Bond.cpp"
#include "Strength.cpp"
#include "Lattice.cpp"
#include "Algorithm.cpp"
#include <iostream>
#include <fstream>
#include <string>
#include <map>

namespace tools
{
    // Function to split strings
    typedef std::pair<std::string, std::string> string_pair;

    string_pair split_string(std::string to_split, std::string delimiter)
    {
        std::size_t d_index = to_split.find(delimiter);
        std::string left, right;
        if(to_split.compare(d_index-1, 1, " ") == 0 || to_split.compare(d_index-1, 1, "\t") == 0)
        {
            left = to_split.substr(0, d_index - 1);
        }
        else
        {
            left = to_split.substr(0, d_index);
        }
        if(to_split.compare(d_index+1, 1, " ") == 0 || to_split.compare(d_index+1, 1, "\t") == 0)
        {
            right = to_split.substr(d_index + 2);
        }
        else
        {
            right = to_split.substr(d_index + 1);
        }

        return std::make_pair(left, right);
    }
}

