#ifndef ANALYSIS_H
#define ANALYSIS_H
#include "Abstract_Classes.h"
#include "Tools.h"
#include <map>

class Statistics
{
public:
    Simulation * sim;
    tools::params params;
    Bond * this_bond;
    double p_c = 0.24885;
    std::map<long int, long int> burst_distribution;
    std::map<long int, long int> mass_r_distribution;
    std::map<long int, long int> mass_l_distribution;
    std::map<std::pair<int, int>, std::map<int, int> > branch_distribution;
    bool network_initialized = false;
    bool parameters_set = false;

    Statistics(Simulation * this_sim);
    Statistics(Simulation * this_sim, tools::params this_params);
    void calculate(void);
    void set_p_c(double value);
    void get_burst_distribution(void);
    void write_burst_distribution_to_file(void);
    void get_mass_r_distribution(void);
    void write_mass_r_distribution_to_file(void);
    void get_mass_l_distribution(void);
    void write_mass_l_distribution_to_file(void);
    void get_branching_statistics(void);
    void write_branching_statistics_to_file(void);

    long int steps_before_fault = 0;
    long int on_fault_count = 0;
    long int off_fault_count = 0;
    double fault_fraction;

    void get_fault_data(void);
    void write_fault_data_to_file(void);
};


#endif // ANALYSIS_H
