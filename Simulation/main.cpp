#include <time.h>
#include <math.h>
#include "Abstract_Classes.h"
#include "Site.cpp"
#include "Bond.cpp"
#include "Strength.cpp"
#include "Lattice.cpp"
#include "Algorithm.cpp"

class Statistics
{
public:
    Simulation * sim;
    Statistics(Simulation * this_sim)
    {
        sim = this_sim;
    }

    Bond * this_bond;
    double p_c = 0.24885;
    std::map<long int, long int> burst_distribution;
    std::map<long int, long int> mass_r_distribution;
    std::map<long int, long int> mass_l_distribution;
    std::map<std::pair<int, int>, std::map<int, int> > branch_distribution;

    bool network_initialized = false;

    void set_p_c(double value)
    {
        p_c = value;
    }

    void get_burst_distribution(void)
    {
        long int size_of_this_burst = 0;

        for(this_bond = sim->Algorithm_ptr->get_network_begin();
            sim->Algorithm_ptr->more_network(); this_bond = sim->Algorithm_ptr->get_network_next())
        {
            if(this_bond->get_strength() < p_c)
            {
                size_of_this_burst++;
            }
            else
            {
                if(burst_distribution.count(size_of_this_burst) > 0)
                {
                    burst_distribution[size_of_this_burst]++;
                }
                else
                {
                    burst_distribution.insert(std::make_pair(size_of_this_burst, 1));
                    size_of_this_burst = 0;
                }
            }
        }
    }

    void write_burst_distribution_to_file(void)
    {
        std::ofstream toFile("burst_distribution.txt", std::ios::trunc);


        toFile << burst_distribution.size() << "\n";

        toFile << "Example Burst Distribution p_c = " << p_c << "\n";

        std::map<long int, long int>::iterator size;

        for(size = burst_distribution.begin(); size != burst_distribution.end(); size++)
        {
            toFile << size->first << "\t" << size->second << "\n";
        }

    }

    void get_mass_r_distribution(void)
    {
        if(network_initialized == false)
        {
            sim->Lattice_ptr->initialize_network();
            network_initialized = true;
        }

        mass_r_distribution.clear();

        for(this_bond = sim->Algorithm_ptr->get_network_begin();
            sim->Algorithm_ptr->more_network(); this_bond = sim->Algorithm_ptr->get_network_next())
        {
            double distance = ceil(sim->Lattice_ptr->get_euclidian_distance(this_bond->second));


            if(mass_r_distribution.count(distance) > 0)
            {
                mass_r_distribution[distance]++;
            }
            else
            {
                mass_r_distribution.insert(std::make_pair(distance, 1));
            }

        }

    }

    void write_mass_r_distribution_to_file(void)
    {
        std::ofstream toFile("mass_r_distribution.txt", std::ios::trunc);


        toFile << mass_r_distribution.size() << "\n";

        toFile << "Example Mass_r Distribution\n";

        std::map<long int, long int>::iterator radius;

        for(radius = mass_r_distribution.begin(); radius != mass_r_distribution.end(); radius++)
        {
            toFile << radius->first << "\t" << radius->second << "\n";
        }

    }

    void get_mass_l_distribution(void)
    {
        if(network_initialized == false)
        {
            sim->Lattice_ptr->initialize_network();
            network_initialized = true;
        }

        mass_l_distribution.clear();

        for(this_bond = sim->Algorithm_ptr->get_network_begin();
            sim->Algorithm_ptr->more_network(); this_bond = sim->Algorithm_ptr->get_network_next())
        {
            double chem_distance = ceil(sim->Lattice_ptr->get_chemical_distance(this_bond->second));


            if(mass_l_distribution.count(chem_distance) > 0)
            {
                mass_l_distribution[chem_distance]++;
            }
            else
            {
                mass_l_distribution.insert(std::make_pair(chem_distance, 1));
            }

        }

    }

    void write_mass_l_distribution_to_file(void)
    {
        std::ofstream toFile("mass_l_distribution.txt", std::ios::trunc);


        toFile << mass_l_distribution.size() << "\n";

        toFile << "Example Mass_l Distribution\n";

        std::map<long int, long int>::iterator chem_level;

        for(chem_level = mass_l_distribution.begin(); chem_level != mass_l_distribution.end(); chem_level++)
        {
            toFile << chem_level->first << "\t" << chem_level->second << "\n";
        }

    }

    void get_branching_statistics(void)
    {
        if(network_initialized == false)
        {
            sim->Lattice_ptr->initialize_network();
            network_initialized = true;
        }

        // Data structures
        std::multiset<Site *> bond_begins;
        typedef std::pair<int, int> branch;
        typedef std::multiset<branch> branch_set;
        typedef std::multiset<branch>::iterator branch_set_iterator;
        std::map<Site *, branch_set> free_ends;
        std::map<Site *, branch_set>::iterator this_end;
        std::deque<std::map<Site *, branch_set>::iterator> free_ends_to_delete;

        // Get a list of bond_ends in the newtwork
        for(this_bond = sim->Algorithm_ptr->get_network_begin();
            sim->Algorithm_ptr->more_network(); this_bond = sim->Algorithm_ptr->get_network_next())
        {

            // std::cout << this_bond->first << std::endl;
            bond_begins.insert(this_bond->first);
        }

        // Find the current free ends
        for(this_bond = sim->Algorithm_ptr->get_network_begin();
            sim->Algorithm_ptr->more_network(); this_bond = sim->Algorithm_ptr->get_network_next())
        {
            if(bond_begins.count(this_bond->second) == 0)
            {
                if(free_ends.count(this_bond->first) > 0)
                {
                    (free_ends[this_bond->first]).insert(branch(1,1));
                }
                else
                {
                    branch_set temp;
                    temp.insert(branch(1, 1));
                    free_ends[this_bond->first] = temp;
                }
            }
        }

        // std::cout << free_ends.size() << std::endl;

        while(free_ends.size() > 1 || (free_ends.begin())->first != sim->Lattice_ptr->get_origin())
        {
            // Loop over ends and work inward
            for(this_end = free_ends.begin(); this_end != free_ends.end(); this_end++)
            {
                if(this_end->first != sim->Lattice_ptr->get_origin())
                {
                    // std::cout << this_end->first << "\t" << (this_end->second).size() << std::endl;

                    // If all ends are present
                    if((this_end->second).size() == bond_begins.count(this_end->first))
                    {

                        // std::cout << "\t All ends present." << std::endl;

                        // Combine all ends
                        while((this_end->second).size() > 1)
                        {

                            // Get two random ends
                            branch_set_iterator iter_1 = (this_end->second).begin();

                            int random_index = floor((this_end->second).size()*drand48());
                            int current_index = 0;

                            // std::cout << "\t" << random_index;

                            while(random_index != current_index)
                            {
                                iter_1++;
                                current_index++;
                            }

                            while(random_index == current_index)
                            {
                                random_index = floor((this_end->second).size()*drand48());
                            }

                            // std::cout << "\t" << random_index << std::endl;

                            branch_set_iterator iter_2 = (this_end->second).begin();

                            current_index = 0;
                            while(random_index != current_index)
                            {
                                iter_2++;
                                current_index++;
                            }

                            // Get the order and length of the ends
                            int order_1 = iter_1->first;
                            int size_1 = iter_1->second;
                            int order_2 = iter_2->first;
                            int size_2 = iter_2->second;

                            // std::cout << "\t Combining ends: " << order_1 << "\t" << order_2 << std::endl;

                            // If the ends have the same order
                            if(order_1 == order_2)
                            {
                                // std::cout << "\t Ends have same order." << std::endl;

                                std::pair<int, int> branch_order(order_1, order_2);

                                // Add them to the branch distribution
                                if(branch_distribution.count(branch_order) > 0)
                                {
                                    if((branch_distribution[branch_order]).count(size_1) > 0)
                                    {
                                        (branch_distribution[branch_order])[size_1] += 1;
                                    }
                                    else
                                    {
                                        (branch_distribution[branch_order]).insert(std::make_pair(size_1, 1));
                                    }

                                    if((branch_distribution[branch_order]).count(size_2) > 0)
                                    {
                                        (branch_distribution[branch_order])[size_2] += 1;
                                    }
                                    else
                                    {
                                        (branch_distribution[branch_order]).insert(std::make_pair(size_2, 1));
                                    }
                                }
                                else
                                {
                                    std::map<int, int> temp;
                                    branch_distribution.insert(std::make_pair(branch_order, temp));

                                    if(size_1 == size_2)
                                    {
                                        // std::cout << "\t Branch orders stored" << std::endl;

                                        (branch_distribution[branch_order]).insert(std::make_pair(size_1, 2));
                                    }
                                    else
                                    {
                                        (branch_distribution[branch_order]).insert(std::make_pair(size_1, 1));
                                        (branch_distribution[branch_order]).insert(std::make_pair(size_2, 1));
                                    }
                                }

                                // Remove ends from list


                                (this_end->second).erase(iter_1);
                                // std::cout << "\t First end erased.\t " << std::endl;

                                (this_end->second).erase(iter_2);
                                // std::cout << "\t Second end erased." << std::endl;

                                // Add a new end with order one more and length 0
                                (this_end->second).insert(std::make_pair(order_1+1, 0));
                            }
                            else if (order_1 > order_2)
                            {
                                // Add end_2 to branch distribution
                                std::pair<int, int> branch_order(order_2, order_1);

                                if(branch_distribution.count(branch_order) > 0)
                                {
                                    if((branch_distribution[branch_order]).count(size_2) > 0)
                                    {
                                        (branch_distribution[branch_order])[size_2] += 1;
                                    }
                                    else
                                    {
                                        (branch_distribution[branch_order]).insert(std::make_pair(size_2, 1));
                                    }
                                }
                                else
                                {
                                    std::map<int, int> temp;
                                    branch_distribution.insert(std::make_pair(branch_order, temp));
                                    (branch_distribution[branch_order]).insert(std::make_pair(size_2, 1));
                                }

                                // Remove end_2 from list
                                (this_end->second).erase(iter_2);
                            }
                            else
                            {
                                // Add end_1 to branch distribution
                                std::pair<int, int> branch_order(order_1, order_2);

                                if(branch_distribution.count(branch_order) > 0)
                                {
                                    if((branch_distribution[branch_order]).count(size_1) > 0)
                                    {
                                        (branch_distribution[branch_order])[size_1] += 1;
                                    }
                                    else
                                    {
                                        (branch_distribution[branch_order]).insert(std::make_pair(size_1, 1));
                                    }
                                }
                                else
                                {
                                    std::map<int, int> temp;
                                    branch_distribution.insert(std::make_pair(branch_order, temp));
                                    (branch_distribution[branch_order]).insert(std::make_pair(size_1, 1));
                                }

                                // Remove end_2 from list
                                (this_end->second).erase(iter_1);
                            }
                        }

                        Site * new_end = sim->Lattice_ptr->get_upstream(this_end->first);


                        branch this_branch = *((this_end->second).begin());
                        this_branch.second += 1;


                        if(free_ends.count(new_end) > 0)
                        {
                            (free_ends[new_end]).insert(this_branch);
                        }
                        else
                        {
                            branch_set temp;
                            free_ends[new_end] = temp;
                            (free_ends[new_end]).insert(this_branch);
                        }


                        free_ends_to_delete.push_back(this_end);
                    }
                }

                // Delete old ends
                std::deque<std::map<Site *, branch_set>::iterator>::iterator i;

                for(i = free_ends_to_delete.begin(); i != free_ends_to_delete.end(); i++)
                {
                    free_ends.erase(*i);
                }
                free_ends_to_delete.clear();
            }
        }

        branch_set at_origin = (free_ends.begin())->second;

        for(branch_set_iterator i = at_origin.begin(); i != at_origin.end(); i++)
        {
            std::pair<int, int> branch_order(i->first, i->first);

            if(branch_distribution.count(branch_order) > 0)
            {
                if((branch_distribution[branch_order]).count(i->second) > 0)
                {
                    (branch_distribution[branch_order])[i->second] += 1;
                }
                else
                {
                    (branch_distribution[branch_order]).insert(std::make_pair(i->second, 1));
                }
            }
            else
            {
                std::map<int, int> temp;
                branch_distribution.insert(std::make_pair(branch_order, temp));
                (branch_distribution[branch_order]).insert(std::make_pair(i->second, 1));
            }


        }
    }

    void write_branching_statistics_to_file(void)
    {
        std::ofstream toFile("branch_distribution.txt", std::ios::trunc);

        toFile << (branch_distribution.rbegin())->first.first << "\n";

        toFile << "Example Branch Distribution\n";

        std::map<std::pair<int, int>, std::map<int, int> >::iterator this_order;

        for(this_order = branch_distribution.begin(); this_order != branch_distribution.end(); this_order++)
        {
            std::map<int, int>::iterator i;
            for(i = (this_order->second).begin(); i != (this_order->second).end(); i++)
            {
                toFile << this_order->first.first << "\t" << this_order->first.second << "\t" << i->first << "\t" << i->second << "\n";
            }
        }

        toFile.close();

    }

    long int steps_before_fault = 0;
    long int on_fault_count = 0;
    long int off_fault_count = 0;
    double fault_fraction;

    void get_fault_data(void)
    {

        bool at_fault = false;

        for(this_bond = sim->Algorithm_ptr->get_network_begin();
            sim->Algorithm_ptr->more_network(); this_bond = sim->Algorithm_ptr->get_network_next())
        {
            if(! at_fault)
            {
                if(sim->Lattice_ptr->on_any_fault(this_bond))
                {
                    on_fault_count++;
                    at_fault = true;
                    fault_fraction = sim->Lattice_ptr->get_fault_fraction(this_bond);
                }
                else
                {
                    steps_before_fault++;
                }
            }
            else
            {
                if(sim->Lattice_ptr->on_any_fault(this_bond))
                {
                    on_fault_count++;
                }
                else
                {
                    off_fault_count++;
                }
            }
        }
    }

    void write_fault_data_to_file(void)
    {
        std::ofstream toFile("fault_data.txt", std::ios::trunc);

        toFile << "Fault Fraction: " << fault_fraction << "\n";
        toFile << "Steps Before Fault: " << steps_before_fault << "\n";
        toFile << "On Fault Count: " << on_fault_count << "\n";
        toFile << "Off Fault Count: " << off_fault_count << "\n";
    }

};

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

    //current_sim.Algorithm_ptr->write_sim_to_file();
    //current_sim.Algorithm_ptr->free_up_memory();

    Statistics * stats;

    stats = new Statistics(&current_sim);
    stats->set_p_c(p_c);

    stats->get_burst_distribution();

    stats->write_burst_distribution_to_file();

    stats->get_mass_r_distribution();

    stats->write_mass_r_distribution_to_file();

    //stats->get_mass_l_distribution();

    //stats->write_mass_l_distribution_to_file();

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

int main(int argc, char **argv)
{
    long int size_of_run = atoi(argv[1]);

    normal_run(size_of_run, 0.5);

    return 0;
}




