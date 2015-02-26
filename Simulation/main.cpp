#include <map>
#include <set>
#include <deque>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <cmath>

#define INDEX_TYPE short int

// Forward declare Simulation class
class Simulation;

#include "Bond.h"

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

    // Would like to eliminate
    virtual bool on_any_fault(Bond * bond)=0;
    virtual double get_fault_fraction(Bond * bond)=0;
    virtual Site * get_upstream(Site *)=0;
};

// Abstract class for simulation
class Algorithm
{
public:
    Simulation * sim;
    virtual void initialize_sim(void)=0;
    virtual void advance_sim(void)=0;
    virtual void free_up_memory(void)=0;
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

// Random number generator seed
unsigned long long rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}

// Uniformly distributed bond strengths
class uniform_Strength: public Strength
{
public:
    uniform_Strength()
    {
    }

    double get_new_strength(void)
    {
        return drand48();
    }
};

// 2D infinite square lattice
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

    Site * get_upstream(Site * site)
    {
        return network[site].site;
    }

    void write_bond(std::ofstream & file, Bond * & bond)
    {
        location site = sites_by_ptr[bond->first];

        file << bond->get_strength() << "\t" << site.loc[0] << "\t" << site.loc[1] << "\t";

        site = sites_by_ptr[bond->second];

        file << site.loc[0] << "\t" << site.loc[1] << "\n";

    }


    bool on_any_fault(Bond * bond)
    {
        return false;
    }
    double get_fault_fraction(Bond * bond)
    {
        return 1;
    }
};

// 3D infinite cubic lattice
class unbound_cubic_Lattice_3D: public Lattice
{
private:
    class location
    {
    public:
        int loc[3];

        location() {};
        location(int first, int second, int third)
        {
            loc[0] = first;
            loc[1] = second;
            loc[2] = third;
        }
        bool operator< (const location & other) const
        {
            if(loc[0] < other.loc[0]) return true;
            if(other.loc[0] < loc[0]) return false;
            if(loc[1] < other.loc[1]) return true;
            if(other.loc[1] < loc[1]) return false;
            if(loc[2] < other.loc[2]) return true;
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
    ~unbound_cubic_Lattice_3D()
    {
        std::map<Site *, location>::iterator to_delete;

        for(to_delete = sites_by_ptr.begin(); to_delete != sites_by_ptr.end(); to_delete++)
        {
            delete to_delete->first;
        }
    }

    Site * get_origin(void)
    {
        location origin(0, 0, 0);

        return get_site_ptr(origin);
    }

    void set_current_site(Site * site)
    {
        current_site = sites_by_ptr[site];
        next_neighbor = 0;
    }

    bool more_neighbors()
    {
        if(next_neighbor < 6) return true;

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
        case 4: // Front
            next_neighbor++;
            neighbor.loc[2] += 1;
            return get_site_ptr(neighbor, new_site);
        case 5: // Front
            next_neighbor++;
            neighbor.loc[2] -= 1;
            return get_site_ptr(neighbor, new_site);
        }
    }

    double get_euclidian_distance(Site * site_1, Site * site_2 = NULL)
    {
        double distance;
        location loc_1 = sites_by_ptr[site_1];


        if (site_2 == NULL)
        {

            distance = sqrt(loc_1.loc[0]*loc_1.loc[0] + loc_1.loc[1]*loc_1.loc[1] + loc_1.loc[2]*loc_1.loc[2]);
        }
        else
        {
            location loc_2 = sites_by_ptr[site_2];

            double diff_x = loc_1.loc[0] - loc_2.loc[0];
            double diff_y = loc_1.loc[1] - loc_2.loc[1];
            double diff_z = loc_1.loc[2] - loc_2.loc[2];

            distance = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);
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

    Site * get_upstream(Site * site)
    {
        return network[site].site;
    }

    void write_bond(std::ofstream & file, Bond * & bond)
    {
        location site = sites_by_ptr[bond->first];

        file << bond->get_strength() << "\t" << site.loc[0] << "\t" << site.loc[1] << "\t" << site.loc[2] << "\t";

        site = sites_by_ptr[bond->second];

        file << site.loc[0] << "\t" << site.loc[1] << "\t" << site.loc[2] << "\n";

    }


    bool on_any_fault(Bond * bond)
    {
        return false;
    }
    double get_fault_fraction(Bond * bond)
    {
        return 1;
    }
};


// 2D infinite square lattice with simple faults
class unbound_square_Lattice_2D_with_faults: public Lattice
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

    class Simple_Fault
    {
    private:
        unbound_square_Lattice_2D_with_faults * lattice;
    public:
        Simple_Fault(double this_ID, unbound_square_Lattice_2D_with_faults * this_lattice,
                     double this_fraction, long int this_distance,
                     char this_orientation)
        {
            ID = this_ID;
            lattice = this_lattice;
            fraction = this_fraction;
            distance = this_distance;
            orientation = this_orientation;
        }

        int ID;
        double fraction;
        long int distance;
        char orientation;

        bool operator< (const Simple_Fault & other) const
        {
            return ID < other.ID;
        }

        bool on_fault(Bond * bond) const
        {
            if(orientation == 'v')
            {
                if(lattice->sites_by_ptr[bond->first].loc[0] ==
                        lattice->sites_by_ptr[bond->second].loc[0])
                {
                    long int x_plane = lattice->sites_by_ptr[lattice->get_origin()].loc[0] + distance;

                    if(lattice->sites_by_ptr[bond->first].loc[0] == x_plane)
                    {
                        return true;
                    }

                }

                return false;
            }
            if(orientation == 'h')
            {
                if(lattice->sites_by_ptr[bond->first].loc[1] ==
                        lattice->sites_by_ptr[bond->second].loc[1])
                {
                    long int y_plane = lattice->sites_by_ptr[lattice->get_origin()].loc[1] + distance;

                    if(lattice->sites_by_ptr[bond->first].loc[1] == y_plane)
                    {
                        return true;
                    }

                }

                return false;
            }

        }
    };

    std::set<Simple_Fault> faults;

    bool on_any_fault(Bond * bond)
    {
        std::set<Simple_Fault>::iterator this_fault;

        for(this_fault = faults.begin(); this_fault != faults.end(); this_fault++)
        {
            if(this_fault->on_fault(bond))
            {
                return true;
            }
        }

        return false;
    }

    double get_fault_fraction(Bond * bond)
    {
        std::set<Simple_Fault>::iterator this_fault;

        // Currently ignores intersections
        for(this_fault = faults.begin(); this_fault != faults.end(); this_fault++)
        {
            if(this_fault->on_fault(bond))
            {
                return this_fault->fraction;
            }
        }

        return 1;
    }


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
    ~unbound_square_Lattice_2D_with_faults()
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

    Site * get_upstream(Site * site)
    {
        return network[site].site;
    }

    void write_bond(std::ofstream & file, Bond * & bond)
    {
        location site = sites_by_ptr[bond->first];

        file << bond->get_strength() << "\t" << site.loc[0] << "\t" << site.loc[1] << "\t";

        site = sites_by_ptr[bond->second];

        file << site.loc[0] << "\t" << site.loc[1] << "\t";

        file << on_any_fault(bond) << "\n";

    }

    int ID_counter = 0;
    void add_fault(double fraction, long int distance, char orientaion)
    {
        faults.insert(Simple_Fault(ID_counter++, this, fraction, distance, orientaion));
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
                Bond * new_bond = sim->Bond_ptr->make_bond(sim, site_ptr, neighbor);

                if(sim->Lattice_ptr->on_any_fault(new_bond))
                {
                    new_bond->set_strength_to(sim->Lattice_ptr->get_fault_fraction(new_bond)*new_bond->get_strength());
                }

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
    double p_c = 0.499;
    std::map<long int, long int> burst_distribution;
    std::map<long int, long int> mass_r_distribution;
    std::map<long int, long int> mass_l_distribution;
    std::map<std::pair<int, int>, std::map<int, int> > branch_distribution;

    bool network_initialized = false;

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

int main(int argc, char **argv)
{
    time_t start, end;

    time(&start);

    // Give random number generator a seed
    srand48(rdtsc());

    long int N = atoi(argv[1]);

    Simulation current_sim;

    ip_central_Algorithm algorithm;
    algorithm.sim = & current_sim;
    current_sim.Algorithm_ptr = & algorithm;

    // unbound_square_Lattice_2D_with_faults lattice;    
    // lattice.add_fault(atof(argv[2]), atoi(argv[3]), 'v');

    unbound_square_Lattice_2D lattice;

    // unbound_cubic_Lattice_3D lattice;
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

    for(long int i=0; i<N; i++)
    {
        current_sim.Algorithm_ptr->advance_sim();
    }

    // current_sim.Algorithm_ptr->write_sim_to_file();
    current_sim.Algorithm_ptr->free_up_memory();

    Statistics * stats;

    stats = new Statistics(&current_sim);

    stats->get_burst_distribution();

    stats->write_burst_distribution_to_file();

    stats->get_mass_r_distribution();

    stats->write_mass_r_distribution_to_file();

    stats->get_mass_l_distribution();

    stats->write_mass_l_distribution_to_file();

    stats->get_branching_statistics();

    stats->write_branching_statistics_to_file();

    delete stats;

    time(&end);

    std::cout << difftime(end,start) << std::endl;

    return 0;
}




