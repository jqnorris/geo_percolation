#ifndef TOOLS_H
#define TOOLS_H

#include "Abstract_Classes.h"
#include <string>

namespace tools
{
    // Function to load sim parameters from file
    void configure_from_file(Simulation & sim);

    // Function to split strings
    typedef std::pair<std::string, std::string> string_pair;
    string_pair split_string(std::string to_split, std::string delimiter);
}


#endif // TOOLS_H
