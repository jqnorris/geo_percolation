
#include <map>
#include <set>
#include <deque>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <ctime>
#include <fstream>
#include <cmath>

#define INDEX_TYPE short int

// Forward declare Simulation class
class Simulation;

// Abstract class for sites
class Site
{
public:
    virtual bool is_occupied(void)const=0;
    virtual void set_occupied_to(bool state)=0;
    virtual Site * make_site(Simulation * sim)=0;
    virtual ~Site() {};
};

// Abstract class for bonds
class Bond
{
public:
    Site * first;
    Site * second;
    virtual double get_strength(void)const=0;
    virtual void set_strength_to(double value)=0;
    virtual double get_time_to_failure(void)const=0;
    virtual void set_time_to_failure(double time)=0;
    virtual Bond * make_bond(Simulation * sim, Site *, Site *)=0;
    virtual ~Bond() {};
};


// Abstract class for lattice of sites/bonds
class Lattice
{
public:
    Simulation * sim;
    virtual Site * get_origin(void)=0;
    virtual void set_current_site(Site * this_site)=0;
    virtual bool more_neighbors(void)=0;
    virtual Site * get_next_neighbor(int * =NULL)=0;
    virtual double get_euclidian_distance(Site * site_1, Site * site_2 = NULL)=0;
    virtual void initialize_network(void)=0;
    virtual long int get_chemical_distance(Site * site_1, Site * site_2 = NULL)=0;
    virtual void write_bond(std::ofstream &, Bond * &)=0;
};

// Abstract class for simulation
class Algorithm
{
public:
    Simulation * sim;
    virtual void initialize_sim(void)=0;
    virtual void advance_sim(void)=0;
    virtual void write_sim_to_file()=0;
    virtual Bond * get_network_begin(void)=0;
    virtual Bond * get_network_next(void)=0;
    virtual bool more_network(void)=0;
};

// Abstract class for generating times to failure
class Timing
{
public:
    virtual double new_time_to_failure(void)=0;
};

// Abstract class for generating bond strengths
class Strength
{
public:
    virtual double get_new_strength(void)=0;
};

class Simulation
{
public:
    Algorithm * Algorithm_ptr;
    Lattice * Lattice_ptr;
    Site * Site_ptr;
    Bond * Bond_ptr;
    Timing * Timing_ptr;
    Strength * Strength_ptr;
};


// Simplest type site
class simple_Site: public Site
{
private:
    bool occupied;

public:
    ~simple_Site() {};
    bool is_occupied() const
    {
        return occupied;
    }

    void set_occupied_to(bool state)
    {
        occupied = state;
    }

    Site * make_site(Simulation * sim)
    {
        simple_Site * temp = new simple_Site;

        temp->occupied = false;
        return temp;
    }

};

// Simplest type of bond
class simple_Bond: public Bond
{
private:
    double strength;

public:
    ~simple_Bond() {};
    double get_strength() const
    {
        return strength;
    }
    void set_strength_to(double value)
    {
        strength = value;
    }
    double get_time_to_failure(void)const
    {
        return 0;
    }
    virtual void set_time_to_failure(double time)
    {
        return;
    }
    Bond * make_bond(Simulation * sim, Site * first, Site * second)
    {
        simple_Bond * temp = new simple_Bond;
        temp->first = first;
        temp->second = second;
        temp->set_strength_to(sim->Strength_ptr->get_new_strength());
        return temp;
    }
};

// Uniformly distributed bond strengths
class uniform_Strength: public Strength
{
public:
    uniform_Strength()
    {
        srand48(time(0));
    }

    double get_new_strength(void)
    {
        return drand48();
    }
};

// 2D square lattice
class unbound_square_Lattice_2D: public Lattice
{
private:
    class location
    {
    public:
        int loc[2];

        location() {};
        location(int first, int second)
        {
            loc[0] = first;
            loc[1] =second;
        }
        bool operator< (const location & other) const
        {
            if(loc[0] < other.loc[0]) return true;
            if(other.loc[0] < loc[0]) return false;
            if(loc[1] < other.loc[1]) return true;
            return false;
        }
    };

    location current_site;
    int next_neighbor;
    std::map<location, Site *> sites_by_loc;
    std::map<Site *, location> sites_by_ptr;

    class Parent
    {
    public:
        Site * site;
        long int chemical_level;
        void update(Site* this_site, long int this_chemical_level)
        {
            site = this_site;
            chemical_level = this_chemical_level;

        }
        Parent()
        {
            site = NULL;
            chemical_level = -1;
        }

        Parent(Site* this_site, long int this_chemical_level)
        {
            update(this_site, this_chemical_level);
        }
    };

    std::map<Site *, Parent> network;

    Site * get_site_ptr(location loc, int * new_site = NULL)
    {
        Site * loc_ptr;

        if(sites_by_loc.count(loc) > 0)
        {
            loc_ptr = sites_by_loc[loc];

            if(new_site != NULL)
            {
                *new_site = 0;
            }
        }
        else
        {
            loc_ptr = sim->Site_ptr->make_site(sim);
            sites_by_loc.insert(std::make_pair(loc, loc_ptr));
            sites_by_ptr.insert(std::make_pair(loc_ptr, loc));

            if(new_site != NULL)
            {
                *new_site = 1;
            }
        }

        return loc_ptr;
    }

public:
    ~unbound_square_Lattice_2D()
    {
        std::map<Site *, location>::iterator to_delete;

        for(to_delete = sites_by_ptr.begin(); to_delete != sites_by_ptr.end(); to_delete++)
        {
            delete to_delete->first;
        }
    }

    Site * get_origin(void)
    {
        location origin(0, 0);

        return get_site_ptr(origin);
    }

    void set_current_site(Site * site)
    {
        current_site = sites_by_ptr[site];
        next_neighbor = 0;
    }

    bool more_neighbors()
    {
        if(next_neighbor < 4) return true;

        return false;
    }

    Site * get_next_neighbor(int * new_site = NULL)
    {
        location neighbor = current_site;

        switch(next_neighbor)
        {
        case 0: // Up
            next_neighbor++;
            neighbor.loc[1] += 1;
            return get_site_ptr(neighbor, new_site);
        case 1: // Down
            next_neighbor++;
            neighbor.loc[1] -= 1;
            return get_site_ptr(neighbor, new_site);
        case 2: // Left
            next_neighbor++;
            neighbor.loc[0] -= 1;
            return get_site_ptr(neighbor, new_site);
        case 3: // Right
            next_neighbor++;
            neighbor.loc[0] += 1;
            return get_site_ptr(neighbor, new_site);
        }
    }

    double get_euclidian_distance(Site * site_1, Site * site_2 = NULL)
    {
        double distance;
        location loc_1 = sites_by_ptr[site_1];


        if (site_2 == NULL)
        {

            distance = sqrt(loc_1.loc[0]*loc_1.loc[0] + loc_1.loc[1]*loc_1.loc[1]);
        }
        else
        {
            location loc_2 = sites_by_ptr[site_2];

            double diff_x = loc_1.loc[0] - loc_2.loc[0];
            double diff_y = loc_1.loc[1] - loc_2.loc[1];

            distance = sqrt(diff_x*diff_x + diff_y*diff_y);
        }

        return distance;
    }

    void initialize_network(void)
    {
        Parent temp;

        for(Bond * current_bond = sim->Algorithm_ptr->get_network_begin();
            sim->Algorithm_ptr->more_network();
            current_bond = sim->Algorithm_ptr->get_network_next())
        {
            if (current_bond->first == get_origin())
            {
                temp.update(current_bond->first, 0);
            }
            else
            {
                temp.update(current_bond->first, (network[current_bond->first]).chemical_level+1);
            }
            network.insert(std::make_pair(current_bond->second, temp));
        }

    }

    long int get_chemical_distance(Site * site_1, Site * site_2 = NULL)
    {
        if(site_2 == NULL or site_2 == get_origin())
        {
            return network[site_1].chemical_level;
        }
        else if (site_1 == get_origin())
        {
            return network[site_2].chemical_level;
        }
        else
        {
            Parent ancestor_1 = network[site_1];
            Parent ancestor_2 = network[site_2];


            long int chem_level_1 = ancestor_1.chemical_level+1;
            long int chem_level_2 =  ancestor_2.chemical_level+1;


            while(ancestor_1.chemical_level > ancestor_2.chemical_level)
            {
                ancestor_1 = network[ancestor_1.site];

            }
            while(ancestor_1.chemical_level < ancestor_2.chemical_level)
            {
                ancestor_2 = network[ancestor_2.site];

            }

            while(ancestor_1.site != ancestor_2.site)
            {
                ancestor_1 = network[ancestor_1.site];
                ancestor_2 = network[ancestor_2.site];
            }

            return (chem_level_1 - ancestor_1.chemical_level) + (chem_level_2 - ancestor_2.chemical_level);
        }
    }

    void write_bond(std::ofstream & file, Bond * & bond)
    {
        location site = sites_by_ptr[bond->first];

        file << bond->get_strength() << "\t" << site.loc[0] << "\t" << site.loc[1] << "\t";

        site = sites_by_ptr[bond->second];

        file << site.loc[0] << "\t" << site.loc[1] << "\n";

    }
};

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
                available_bonds.insert(sim->Bond_ptr->make_bond(sim, site_ptr, neighbor));
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
        return *(current_network_iterator++);
    }

    Bond * get_network_next(void)
    {
        return *(current_network_iterator++);
    }

    bool more_network(void)
    {
        return current_network_iterator != invaded_bonds.end();
    }
};


class testing_Strength: public Strength
{
private:
    int next_strength;

public:
    testing_Strength()
    {
        next_strength=0;
    }

    double get_new_strength()
    {
        return next_strength++;
    }

};

class Statistics
{
public:
    Simulation * sim;
    Statistics(Simulation * this_sim)
    {
        sim = this_sim;
    }

    Bond * this_bond;
    double p_c = 0.4975;
    std::map<long int, long int> burst_distribution;
    std::map<long int, long int> mass_r_distribution;
    std::map<long int, long int> mass_l_distribution;

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

        toFile << "Example Burst Distribution\n";

        std::map<long int, long int>::iterator size;

        for(size = burst_distribution.begin(); size != burst_distribution.end(); size++)
        {
            toFile << size->first << "\t" << size->second << "\n";
        }

    }

    void get_mass_r_distribution(void)
    {
        sim->Lattice_ptr->initialize_network();

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

};



int main(int argc, char **argv)
{
    long int N = atoi(argv[1]);

    Simulation current_sim;

    ip_central_Algorithm algorithm;
    algorithm.sim = & current_sim;
    current_sim.Algorithm_ptr = & algorithm;

    unbound_square_Lattice_2D lattice;
    lattice.sim = & current_sim;
    current_sim.Lattice_ptr = & lattice;

    uniform_Strength strength;
    current_sim.Strength_ptr = & strength;

    simple_Bond bond;
    current_sim.Bond_ptr = & bond;

    simple_Site site;
    current_sim.Site_ptr = & site;



    current_sim.Algorithm_ptr->initialize_sim();

    for(long int i=0; i<N; i++)
    {
        current_sim.Algorithm_ptr->advance_sim();
    }

    current_sim.Algorithm_ptr->write_sim_to_file();

    Statistics * stats;

    stats = new Statistics(&current_sim);

    stats->get_burst_distribution();

    stats->write_burst_distribution_to_file();

    stats->get_mass_r_distribution();

    stats->write_mass_r_distribution_to_file();

    stats->get_mass_l_distribution();

    stats->write_mass_l_distribution_to_file();

    delete stats;

    return 0;
}




