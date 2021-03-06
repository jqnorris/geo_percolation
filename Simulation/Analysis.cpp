#include "Abstract_Classes.h"
#include "Simulation.h"
#include "Analysis.h"
#include <set>
#include <deque>
#include <iostream>
#include <fstream>
#include <math.h>



Statistics::Statistics(Simulation * this_sim)
{
        sim = this_sim;
}

Statistics::Statistics(Simulation * this_sim, tools::params this_params)
{
        sim = this_sim;
        params = this_params;
        parameters_set = true;
}

void Statistics::calculate(void)
{
    if(parameters_set)
    {
        if(params.count("fractures") > 0)
        {
            if(params["fractures"] == "true")
            {
                sim->Algorithm_ptr->write_sim_to_file();
            }
        }
        if(params.count("bursts") > 0)
        {
            if(params["bursts"] == "true")
            {
                set_p_c(std::stof(params["p_c"]));
                get_burst_distribution();
                write_burst_distribution_to_file();
            }
        }
        if(params.count("mass_r") > 0)
        {
            if(params["mass_r"] == "true")
            {
                get_mass_r_distribution();
                write_mass_r_distribution_to_file();
            }
        }
        if(params.count("mass_l") > 0)
        {
            if(params["mass_l"] == "true")
            {
                get_mass_l_distribution();
                write_mass_l_distribution_to_file();
            }
        }
        if(params.count("branching") > 0)
        {
            if(params["branching"] == "true")
            {
                get_branching_statistics();
                write_branching_statistics_to_file();
            }
        }
        if(params.count("faults") > 0)
        {
            if(params["faults"] == "true")
            {
                get_fault_data();
                write_fault_data_to_file();
            }
        }

    }


}

void Statistics::set_p_c(double value)
{
    p_c = value;
}

void Statistics::get_burst_distribution(void)
{
    long int burst_size = 0;

    for(this_bond = sim->Algorithm_ptr->get_network_begin();
        sim->Algorithm_ptr->more_network(); this_bond = sim->Algorithm_ptr->get_network_next())
    {
        if(this_bond->get_strength() < p_c)
        {
            if(burst_size > 0)
            {
                burst_size++;
            }
            else
            {
                burst_size = 1;
            }
        }
        else
        {
            if(burst_distribution.count(burst_size) > 0)
            {
                burst_distribution[burst_size]++;
            }
            else
            {
                burst_distribution.insert(std::make_pair(burst_size, 1));
            }
            burst_size = 0;
        }
    }
}

void Statistics::write_burst_distribution_to_file(void)
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

void Statistics::get_mass_r_distribution(void)
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

void Statistics::write_mass_r_distribution_to_file(void)
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

void Statistics::get_mass_l_distribution(void)
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

void Statistics::write_mass_l_distribution_to_file(void)
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

void Statistics::get_branching_statistics(void)
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


void Statistics::write_branching_statistics_to_file(void)
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

void Statistics::get_fault_data(void)
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

void Statistics::write_fault_data_to_file(void)
{
    std::ofstream toFile("fault_data.txt", std::ios::trunc);

    toFile << "Fault Fraction: " << fault_fraction << "\n";
    toFile << "Steps Before Fault: " << steps_before_fault << "\n";
    toFile << "On Fault Count: " << on_fault_count << "\n";
    toFile << "Off Fault Count: " << off_fault_count << "\n";
}

