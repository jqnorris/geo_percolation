#include "Abstract_Classes.h"
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
};
