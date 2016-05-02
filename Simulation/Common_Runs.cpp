#include "Common_Runs.h"
#include <time.h>
#include "Analysis.h"
#include "Tools.h"
#include "Abstract_Classes.h"
#include "Forward_Declarations.h"
#include "Site.cpp"
#include "Bond.cpp"
#include "Strength.cpp"
#include "Lattice.cpp"
#include "Algorithm.cpp"

namespace Common_Runs
{

void normal_run(int size_of_run, double p_c)
{
    time_t start, end;

    time(&start);

    Simulation current_sim;

    ip_central_Algorithm algorithm;

    current_sim.Algorithm_ptr = & algorithm;
    current_sim.Algorithm_ptr->set_sim(&current_sim);

    // unbound_square_Lattice_2D_with_faults lattice;
    // lattice.add_fault(atof(argv[2]), atoi(argv[3]), 'v');

    // Lattice_1D lattice;

    unbound_square_Lattice_2D lattice;

    //unbound_cubic_Lattice_3D lattice;

    //unbound_cubic_Lattice_3D_faults_anisotropy lattice;
    lattice.setup_lattice();

    // unbound_hypercubic_Lattice_6D lattice;

    lattice.sim = & current_sim;

    current_sim.Lattice_ptr = & lattice;

    uniform_Strength strength;
    //testing_Strength strength;
    current_sim.Strength_ptr = & strength;

    simple_Bond bond;

    current_sim.Bond_ptr = & bond;

    simple_Site site;
    current_sim.Site_ptr = & site;


    current_sim.Algorithm_ptr->initialize_sim();

    for(long int i=0; i<size_of_run; i++)
    {
        current_sim.Algorithm_ptr->advance_sim();
    }

    current_sim.Algorithm_ptr->write_sim_to_file();
    current_sim.Algorithm_ptr->free_up_memory();

    Statistics * stats;

    stats = new Statistics(&current_sim);
    stats->set_p_c(p_c);

    stats->get_burst_distribution();

    stats->write_burst_distribution_to_file();

    stats->get_mass_r_distribution();

    stats->write_mass_r_distribution_to_file();

    stats->get_mass_l_distribution();

    stats->write_mass_l_distribution_to_file();

    //stats->get_branching_statistics();

    //stats->write_branching_statistics_to_file();

    delete stats;

    time(&end);

    std::cout << difftime(end,start) << std::endl;
}


void time_to_fault(long int size_of_run, long int number_of_runs)
{
    std::map<long int, long int> time_distribution;

    for(long int i=0; i<number_of_runs; i++)
    {
        std::cout << std::endl;
        std::cout << i << std::endl;

        // Setup Simulation
        Simulation current_sim;

        ip_central_Algorithm algorithm;
        current_sim.Algorithm_ptr = & algorithm;
        current_sim.Algorithm_ptr->set_sim(&current_sim);

        unbound_cubic_Lattice_3D_faults_anisotropy lattice;
        lattice.setup_lattice();
        lattice.sim = & current_sim;
        current_sim.Lattice_ptr = & lattice;

        uniform_Strength strength;
        current_sim.Strength_ptr = & strength;

        simple_Bond bond;
        current_sim.Bond_ptr = & bond;

        simple_Site site;
        current_sim.Site_ptr = & site;

        current_sim.Algorithm_ptr->initialize_sim();

        // Run Simulation
        long int j;
        for(j=0; j<size_of_run; j++)
        {
            current_sim.Algorithm_ptr->advance_sim();
            if( current_sim.Algorithm_ptr->check_growth() == true)
            {
                break;
            }
        }

        // Record Simulation
        if(time_distribution.count(j) > 0)
        {
            time_distribution[j] += 1;
        }
        else
        {
            time_distribution.insert(std::make_pair(j, 1));
        }
    }

    // Write Runs to File
    std::ofstream toFile("time_distribution.txt", std::ios::trunc);

    toFile << time_distribution.size() << "\n";

    toFile << "Distribution of Times to Fault\n";

    std::map<long int, long int>::iterator time;

    for(time = time_distribution.begin(); time != time_distribution.end(); time++)
    {
        toFile << time->first << "\t" << time->second << "\n";
    }
}

void largest_strength(long int size_of_run, long int number_of_runs)
{
    std::multiset<double> largest_strengths;

    for(long int i=0; i<number_of_runs; i++)
    {
        std::cout << std::endl;
        std::cout << i << std::endl;

        // Setup Simulation
        Simulation current_sim;

        ip_central_Algorithm algorithm;
        current_sim.Algorithm_ptr = & algorithm;
        current_sim.Algorithm_ptr->set_sim(&current_sim);

        Lattice_1D lattice;
        // unbound_square_Lattice_2D lattice;
        // unbound_cubic_Lattice_3D lattice;
        lattice.setup_lattice();
        lattice.sim = & current_sim;
        current_sim.Lattice_ptr = & lattice;

        uniform_Strength strength;
        current_sim.Strength_ptr = & strength;

        simple_Bond bond;
        current_sim.Bond_ptr = & bond;

        simple_Site site;
        current_sim.Site_ptr = & site;

        current_sim.Algorithm_ptr->initialize_sim();

        double largest_strength = 0;

        // Run Simulation
        long int j;
        double this_strength;
        for(j=0; j<size_of_run; j++)
        {
            current_sim.Algorithm_ptr->advance_sim();
            this_strength = (current_sim.Algorithm_ptr->get_last_invaded())->get_strength();

            if(this_strength > largest_strength)
            {
                largest_strength = this_strength;
            }
        }

        // Record Simulation
        largest_strengths.insert(largest_strength);
    }

    // Write Runs to File
    std::ofstream toFile("strength_distribution.txt", std::ios::trunc);

    toFile.precision(17);

    toFile << largest_strengths.size() << "\n";

    toFile << "Distribution of Largest Strengths for M = " << size_of_run << "\n";

    std::set<double>::iterator strength;

    for(strength = largest_strengths.begin(); strength != largest_strengths.end(); strength++)
    {
        toFile << *strength << "\n";
    }
}

void strength_distribution(long int size_of_run, long int number_of_runs)
{
    int resolution = 100000;

    long int * counts;

    counts = new long int[resolution] { };

    std::multiset<double> largest_strengths;

    for(long int i=0; i<number_of_runs; i++)
    {
        std::cout << std::endl;
        std::cout << i << std::endl;

        // Setup Simulation
        Simulation current_sim;

        ip_central_Algorithm algorithm;
        current_sim.Algorithm_ptr = & algorithm;
        current_sim.Algorithm_ptr->set_sim(&current_sim);

        // Lattice_1D lattice;
        unbound_square_Lattice_2D lattice;
        // unbound_cubic_Lattice_3D lattice;
        lattice.setup_lattice();
        lattice.sim = & current_sim;
        current_sim.Lattice_ptr = & lattice;

        uniform_Strength strength;
        current_sim.Strength_ptr = & strength;

        simple_Bond bond;
        current_sim.Bond_ptr = & bond;

        simple_Site site;
        current_sim.Site_ptr = & site;

        current_sim.Algorithm_ptr->initialize_sim();

        // Run Simulation
        long int j;
        double this_strength;
        for(j=0; j<size_of_run; j++)
        {
            current_sim.Algorithm_ptr->advance_sim();
            this_strength = (current_sim.Algorithm_ptr->get_last_invaded())->get_strength();

            counts[(int)(resolution*this_strength)] += 1;
        }
    }

    // Write Runs to File
    std::ofstream toFile("strength_distribution.txt", std::ios::trunc);

    toFile.precision(17);

    toFile << resolution << "\n";

    toFile << "Distribution of Bond Strengths for M = " << size_of_run << "\n";

    for(int i = 1; i <= resolution; i++)
    {
        toFile << i*double(1)/double(resolution) << "\t" << counts[i] << "\n";
    }

    toFile.close();

    delete[] counts;
}

void natural_fracking(long int size_of_run)
{
    Simulation current_sim;

    natural_fracking_Algorithm algorithm;
    current_sim.Algorithm_ptr = & algorithm;
    current_sim.Algorithm_ptr->set_sim(&current_sim);

    finite_square_Lattice_2D lattice;
    lattice.sim = & current_sim;
    current_sim.Lattice_ptr = & lattice;

    uniform_Strength strength;
    current_sim.Strength_ptr = & strength;

    simple_Bond bond;
    current_sim.Bond_ptr = & bond;

    simple_Site site;
    current_sim.Site_ptr = & site;

    lattice.setup_lattice();

    std::cout << "Setting up Algorithm." << std::endl;

    current_sim.Algorithm_ptr->initialize_sim();

    std::cout << "Running Sim." << std::endl;

    for(long int i=0; i<size_of_run; i++)
    {
        current_sim.Algorithm_ptr->advance_sim();
    }

    current_sim.Algorithm_ptr->write_sim_to_file();
    current_sim.Algorithm_ptr->free_up_memory();

    Statistics * stats;

    stats = new Statistics(&current_sim);

    stats->get_burst_distribution();

    stats->write_burst_distribution_to_file();

    stats->get_mass_r_distribution();

    stats->write_mass_r_distribution_to_file();

    stats->get_mass_l_distribution();

    stats->write_mass_l_distribution_to_file();

    //stats->get_branching_statistics();

    //stats->write_branching_statistics_to_file();

    delete stats;
}
}
