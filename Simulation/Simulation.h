#ifndef SIMULATION_H
#define SIMULATION_H

#include "Forward_Declarations.h"
#include "Tools.h"
#include <string>

// Container for storing poiters to abstract classes
class Simulation
{
public:
    Algorithm * Algorithm_ptr;
    Lattice * Lattice_ptr;
    Site * Site_ptr;
    Bond * Bond_ptr;
    Timing * Timing_ptr;
    Strength * Strength_ptr;

    long int size_of_run;
    long int number_of_runs;
    double p_c;
    std::string name;

    void configure_sim(tools::params params);
    void run_sim(void);

    ~Simulation();
};

#endif // SIMULATION_H
