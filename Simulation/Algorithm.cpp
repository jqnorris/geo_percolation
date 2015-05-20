#include "Abstract_Classes.h"
#include "Simulation.h"
#include <map>
#include <set>
#include <deque>
#include <fstream>

// Unbound invasion percolation from a central site
class ip_central_Algorithm: public Algorithm
{
private:
    struct bond_comparison
    {
        bool operator()(const Bond * lhs, const Bond * rhs) const
        {
            return lhs->get_strength() < rhs->get_strength();
        }
    };

    std::set<Bond*, bond_comparison> available_bonds;
    std::set<Bond*, bond_comparison>::iterator current_bond;

    std::deque<Bond *> trapped_bonds;
    std::deque<Bond *> invaded_bonds;

    void invade_site(Site * site_ptr)
    {
        sim->Lattice_ptr->set_current_site(site_ptr);

        site_ptr->set_occupied_to(true);

        while(sim->Lattice_ptr->more_neighbors())
        {
            Site * neighbor = sim->Lattice_ptr->get_next_neighbor();

            if(!neighbor->is_occupied())
            {
                Bond * new_bond = sim->Bond_ptr->make_bond(site_ptr, neighbor, sim->Strength_ptr->get_new_strength());
                sim->Lattice_ptr->modify_strength(new_bond);
                modify_strength(new_bond);

                available_bonds.insert(new_bond);
            }
        }

    }

public:
    ~ip_central_Algorithm()
    {
        for(current_bond=available_bonds.begin(); current_bond != available_bonds.end(); current_bond++)
        {
            delete (*current_bond);
        }

        std::deque<Bond *>::iterator to_delete;

        for(to_delete = trapped_bonds.begin(); to_delete != trapped_bonds.end(); to_delete++)
        {
            delete (*to_delete);
        }

        for(to_delete = invaded_bonds.begin(); to_delete != invaded_bonds.end(); to_delete++)
        {
            delete (*to_delete);
        }
    }

    void initialize_sim()
    {
        Site * origin = sim->Lattice_ptr->get_origin();
        invade_site(origin);

    }

    void set_sim(Simulation * to_set)
    {
        sim = to_set;
    }

    void advance_sim()
    {
        current_bond = available_bonds.begin();

        while( (*current_bond)->second->is_occupied())
        {

            trapped_bonds.push_back(*current_bond);
            available_bonds.erase(current_bond);
            current_bond = available_bonds.begin();
        }

        invaded_bonds.push_back(*current_bond);
        invade_site((*current_bond)->second);
        available_bonds.erase(current_bond);

        return;

    }

    void free_up_memory(void)
    {
        available_bonds.clear();
    }

    void write_sim_to_file()
    {
        std::ofstream toFile1("fractures.txt", std::ios::trunc);
        std::ofstream toFile2("trapped.txt", std::ios::trunc);

        toFile1 << invaded_bonds.size() << "\n";
        toFile2 << trapped_bonds.size() << "\n";

        toFile1 << "Example Invasion\n";
        toFile2 << "Example Trapped\n";

        toFile1.precision(17);
        toFile2.precision(17);

        std::deque<Bond *>::iterator to_write;

        for(to_write = invaded_bonds.begin(); to_write != invaded_bonds.end(); to_write++)
        {
            sim->Lattice_ptr->write_bond(toFile1, *to_write);
        }

        for(to_write = trapped_bonds.begin(); to_write != trapped_bonds.end(); to_write++)
        {
            sim->Lattice_ptr->write_bond(toFile2, *to_write);
        }
        toFile1.close();
        toFile2.close();
    }

    std::deque<Bond *>::iterator current_network_iterator;

    Bond * get_network_begin(void)
    {
        current_network_iterator = invaded_bonds.begin();
        return *current_network_iterator;
    }

    Bond * get_network_next(void)
    {
        current_network_iterator++;
        return *current_network_iterator;
    }

    bool more_network(void)
    {
        return current_network_iterator != invaded_bonds.end();
    }

    void modify_strength(Bond * bond){}

    bool check_growth(void)
    {
        return sim->Lattice_ptr->on_any_fault(get_last_invaded());
    }

    Bond * get_last_invaded(void)
    {
        return invaded_bonds.back();
    }

    long int get_breakthrough_count(void)
    {
        return 0;
    }
};

// Natural Fracking
class natural_fracking_Algorithm: public Algorithm
{
private:
    struct bond_comparison
    {
        bool operator()(const Bond * lhs, const Bond * rhs) const
        {
            return lhs->get_strength() < rhs->get_strength();
        }
    };

    long int breakthrough_count = 0;

    std::set<Bond*, bond_comparison> available_bonds;
    std::set<Bond*, bond_comparison>::iterator current_bond;
    std::deque<Bond *> trapped_bonds;
    std::deque<Bond *> invaded_bonds;
    std::map<Site *, long int> cluster_map;
    long int new_cluster_number = 1;
    std::map<long int, long int> cluster_merging;
    std::map<long int, bool> breakout_map;

    long int get_good_cluster_number(long int cluster_number)
    {
        long int good_cluster_number = cluster_number;
        while(good_cluster_number != cluster_merging[good_cluster_number])
        {
            good_cluster_number = cluster_merging[good_cluster_number];
        }

        while(cluster_number != cluster_merging[cluster_number])
        {
            long int next_number = cluster_merging[cluster_number];
            cluster_merging[cluster_number] = good_cluster_number;
            cluster_number = next_number;
        }

        std::cout << good_cluster_number << std::endl;

        return good_cluster_number;
    }

    long int merge_clusters(long int cluster_1, long int cluster_2)
    {
        long int smaller = std::min(get_good_cluster_number(cluster_1),
                get_good_cluster_number(cluster_2));
        long int larger = std::max(get_good_cluster_number(cluster_1),
                get_good_cluster_number(cluster_2));
        cluster_merging[larger] = smaller;

        return smaller;
    }

    void break_bond(Bond * bond_ptr)
    {
        current_bond = available_bonds.begin();
        Site * site_1 = (*current_bond)->first;
        Site * site_2 = (*current_bond)->second;
        long int good_cluster_number;

        if(cluster_map.count(site_1) > 0)
        {
            long int cluster_1 = get_good_cluster_number(cluster_map[site_1]);
            if(cluster_map.count(site_2) > 0)
            {
                long int cluster_2 = get_good_cluster_number(cluster_map[site_2]);

                if(cluster_1 == cluster_2)
                {
                    trapped_bonds.push_back(bond_ptr);
                }
                else
                {
                    invaded_bonds.push_back(bond_ptr);
                    merge_clusters(cluster_1, cluster_2);
                }
            }
            else
            {
                invaded_bonds.push_back(bond_ptr);
                cluster_map.insert(std::make_pair(site_2, cluster_1));
            }
        }
        else
        {
            if(cluster_map.count(site_2) > 0)
            {
                invaded_bonds.push_back(bond_ptr);
                long int cluster_2 = get_good_cluster_number(cluster_map[site_2]);
                cluster_map.insert(std::make_pair(site_1, cluster_2));
            }
            else
            {
                invaded_bonds.push_back(bond_ptr);
                good_cluster_number = new_cluster_number++;
                cluster_merging.insert(std::make_pair(good_cluster_number, good_cluster_number));
                cluster_map.insert(std::make_pair(site_1, good_cluster_number));
                cluster_map.insert(std::make_pair(site_2, good_cluster_number));
            }
        }
    }

public:
    ~natural_fracking_Algorithm()
    {
        for(current_bond=available_bonds.begin(); current_bond != available_bonds.end(); current_bond++)
        {
            delete (*current_bond);
        }

        std::deque<Bond *>::iterator to_delete;


        for(to_delete = trapped_bonds.begin(); to_delete != trapped_bonds.end(); to_delete++)
        {
            delete (*to_delete);
        }

        for(to_delete = invaded_bonds.begin(); to_delete != invaded_bonds.end(); to_delete++)
        {
            delete (*to_delete);
        }
    }

    void initialize_sim()
    {
        std::set<Site *> to_process;
        std::set<Site *> processed;
        std::set<Site *>::iterator this_site;

        to_process.insert(sim->Lattice_ptr->get_origin());

        while(to_process.size() > 0)
        {
            this_site = to_process.begin();

            sim->Lattice_ptr->set_current_site(*this_site);

            while(sim->Lattice_ptr->more_neighbors())
            {
                Site * neighbor = sim->Lattice_ptr->get_next_neighbor();

                if(processed.count(neighbor) < 1)
                {
                    Bond * new_bond = sim->Bond_ptr->make_bond(*this_site, neighbor, sim->Strength_ptr->get_new_strength());
                    sim->Lattice_ptr->modify_strength(new_bond);
                    modify_strength(new_bond);
                    available_bonds.insert(new_bond);

                    if(sim->Lattice_ptr->on_lattice(neighbor))
                    {
                        to_process.insert(neighbor);
                    }
                }
            }

            processed.insert(*this_site);
            to_process.erase(this_site);
        }

    }

    void set_sim(Simulation * to_set)
    {
        sim = to_set;
    }

    void advance_sim()
    {
        if(available_bonds.size() > 0)
        {
            current_bond = available_bonds.begin();

            break_bond(*current_bond);

            //check_breakthrough(*current_bond);

            if(!sim->Lattice_ptr->on_lattice((*current_bond)->second))
            {
                breakthrough_count++;
            }

            available_bonds.erase(current_bond);

            return;
        }

        return;
    }

    void free_up_memory(void)
    {
        available_bonds.clear();
    }

    void write_sim_to_file()
    {
        std::ofstream toFile1("fractures.txt", std::ios::trunc);
        std::ofstream toFile2("trapped.txt", std::ios::trunc);

        toFile1 << invaded_bonds.size() << "\n";
        toFile2 << trapped_bonds.size() << "\n";

        toFile1 << "Example Invasion\n";
        toFile2 << "Example Trapped\n";

        toFile1.precision(17);
        toFile2.precision(17);

        std::deque<Bond *>::iterator to_write;

        for(to_write = invaded_bonds.begin(); to_write != invaded_bonds.end(); to_write++)
        {
            sim->Lattice_ptr->write_bond(toFile1, *to_write);
        }

        for(to_write = trapped_bonds.begin(); to_write != trapped_bonds.end(); to_write++)
        {
            sim->Lattice_ptr->write_bond(toFile2, *to_write);
        }
        toFile1.close();
        toFile2.close();
    }

    std::deque<Bond *>::iterator current_network_iterator;

    Bond * get_network_begin(void)
    {
        current_network_iterator = invaded_bonds.begin();
        return *current_network_iterator;
    }

    Bond * get_network_next(void)
    {
        current_network_iterator++;
        return *current_network_iterator;
    }

    bool more_network(void)
    {
        return current_network_iterator != invaded_bonds.end();
    }

    void modify_strength(Bond * bond){}

    bool check_growth(void)
    {
        return sim->Lattice_ptr->on_any_fault(get_last_invaded());
    }

    Bond * get_last_invaded(void)
    {
        return invaded_bonds.back();
    }

    long int get_breakthrough_count(void)
    {
        return breakthrough_count;
    }
};
