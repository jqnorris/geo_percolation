#ifndef TOOLS_H
#define TOOLS_H

#include "Abstract_Classes.h"
#include <map>
#include <string>

namespace tools
{
    typedef std::map<std::string, std::string> params;
    // Function to load configuration parameters from file
    params load_cofiguration_from_file(void);

    // Function to split strings
    typedef std::pair<std::string, std::string> string_pair;
    string_pair split_string(std::string to_split, std::string delimiter);
}


#endif // TOOLS_H
