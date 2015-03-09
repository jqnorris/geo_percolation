#include "Abstract_Classes.h"
#include <map>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <set>

// 1D infinite line lattice
class Lattice_1D: public Lattice
{
private:
    typedef long long int location;
    location current_site;
    int next_neighbor;
    std::map<location, Site *> sites_by_loc;
    std::map<Site *, location> sites_by_ptr;

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
    Simulation * sim;
    Site * get_origin(void)
    {
        return get_site_ptr(0);
    }

    void set_current_site(Site * site)
    {
        current_site = sites_by_ptr[site];
        next_neighbor = 0;
    }

    bool more_neighbors()
    {
        if(next_neighbor < 2) return true;

        return false;
    }

    Site * get_next_neighbor(int * new_site = NULL)
    {
        location neighbor = current_site;

        switch(next_neighbor)
        {
        case 0: // left
            next_neighbor++;
            return get_site_ptr(neighbor - 1, new_site);
        case 1: // right
            next_neighbor++;
            return get_site_ptr(neighbor + 1 , new_site);
        }
    }

    double get_euclidian_distance(Site * site_1, Site * site_2 = NULL)
    {
        double distance;
        location loc_1 = sites_by_ptr[site_1];


        if (site_2 == NULL)
        {
            distance = abs(loc_1);
        }
        else
        {
            location loc_2 = sites_by_ptr[site_2];

            distance = abs(loc_2 - loc_1);
        }

        return distance;
    }

    void initialize_network(void)
    {

    }

    long int get_chemical_distance(Site * site_1, Site * site_2 = NULL)
    {
        return get_euclidian_distance(site_1, site_2);
    }

    Site * get_upstream(Site * site)
    {
        location loc = sites_by_ptr[site];

        if (loc < 0) return sites_by_loc[loc + 1];
        if (loc > 0) return sites_by_loc[loc - 1];

        return get_origin();
    }

    void write_bond(std::ofstream & file, Bond * & bond)
    {
        location site = sites_by_ptr[bond->first];

        file << bond->get_strength() << "\t" << site;

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

// 4D infinite hypercubic lattice
class unbound_hypercubic_Lattice_4D: public Lattice
{
private:
    class location
    {
    public:
        int loc[4];

        location() {};

        location(int first, int second, int third, int fourth)
        {
            loc[0] = first;
            loc[1] = second;
            loc[2] = third;
            loc[3] = fourth;
        }
        bool operator< (const location & other) const
        {
            if(loc[0] < other.loc[0]) return true;
            if(other.loc[0] < loc[0]) return false;
            if(loc[1] < other.loc[1]) return true;
            if(other.loc[1] < loc[1]) return false;
            if(loc[2] < other.loc[2]) return true;
            if(other.loc[2] < loc[2]) return false;
            if(loc[3] < other.loc[3]) return true;
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
    ~unbound_hypercubic_Lattice_4D()
    {
        std::map<Site *, location>::iterator to_delete;

        for(to_delete = sites_by_ptr.begin(); to_delete != sites_by_ptr.end(); to_delete++)
        {
            delete to_delete->first;
        }
    }

    Site * get_origin(void)
    {
        location origin(0, 0, 0, 0);

        return get_site_ptr(origin);
    }

    void set_current_site(Site * site)
    {
        current_site = sites_by_ptr[site];
        next_neighbor = 0;
    }

    bool more_neighbors()
    {
        if(next_neighbor < 8) return true;

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
        case 6: // Out
            next_neighbor++;
            neighbor.loc[3] += 1;
            return get_site_ptr(neighbor, new_site);
        case 7: // In
            next_neighbor++;
            neighbor.loc[3] -= 1;
            return get_site_ptr(neighbor, new_site);
        }
    }

    double get_euclidian_distance(Site * site_1, Site * site_2 = NULL)
    {
        double distance;
        location loc_1 = sites_by_ptr[site_1];


        if (site_2 == NULL)
        {

            distance = sqrt(loc_1.loc[0]*loc_1.loc[0] + loc_1.loc[1]*loc_1.loc[1]
                    + loc_1.loc[2]*loc_1.loc[2] + loc_1.loc[3]*loc_1.loc[3]);
        }
        else
        {
            location loc_2 = sites_by_ptr[site_2];

            double diff_x = loc_1.loc[0] - loc_2.loc[0];
            double diff_y = loc_1.loc[1] - loc_2.loc[1];
            double diff_z = loc_1.loc[2] - loc_2.loc[2];
            double diff_a = loc_1.loc[3] - loc_2.loc[3];

            distance = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z + diff_a*diff_a);
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

// 5D infinite hypercubic lattice
class unbound_hypercubic_Lattice_5D: public Lattice
{
private:
    class location
    {
    public:
        int loc[5];

        location() {};

        location(int first, int second, int third, int fourth, int fifth)
        {
            loc[0] = first;
            loc[1] = second;
            loc[2] = third;
            loc[3] = fourth;
            loc[4] = fifth;
        }
        bool operator< (const location & other) const
        {
            if(loc[0] < other.loc[0]) return true;
            if(other.loc[0] < loc[0]) return false;
            if(loc[1] < other.loc[1]) return true;
            if(other.loc[1] < loc[1]) return false;
            if(loc[2] < other.loc[2]) return true;
            if(other.loc[2] < loc[2]) return false;
            if(loc[3] < other.loc[3]) return true;
            if(other.loc[3] < loc[3]) return false;
            if(loc[4] < other.loc[4]) return true;
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
    ~unbound_hypercubic_Lattice_5D()
    {
        std::map<Site *, location>::iterator to_delete;

        for(to_delete = sites_by_ptr.begin(); to_delete != sites_by_ptr.end(); to_delete++)
        {
            delete to_delete->first;
        }
    }

    Site * get_origin(void)
    {
        location origin(0, 0, 0, 0, 0);

        return get_site_ptr(origin);
    }

    void set_current_site(Site * site)
    {
        current_site = sites_by_ptr[site];
        next_neighbor = 0;
    }

    bool more_neighbors()
    {
        if(next_neighbor < 10) return true;

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
        case 6: // Out
            next_neighbor++;
            neighbor.loc[3] += 1;
            return get_site_ptr(neighbor, new_site);
        case 7: // In
            next_neighbor++;
            neighbor.loc[3] -= 1;
            return get_site_ptr(neighbor, new_site);
        case 8: // Top
            next_neighbor++;
            neighbor.loc[4] += 1;
            return get_site_ptr(neighbor, new_site);
        case 9: // Bottom
            next_neighbor++;
            neighbor.loc[4] -= 1;
            return get_site_ptr(neighbor, new_site);
        }
    }

    double get_euclidian_distance(Site * site_1, Site * site_2 = NULL)
    {
        double distance;
        location loc_1 = sites_by_ptr[site_1];


        if (site_2 == NULL)
        {

            distance = sqrt(loc_1.loc[0]*loc_1.loc[0] + loc_1.loc[1]*loc_1.loc[1]
                    + loc_1.loc[2]*loc_1.loc[2] + loc_1.loc[3]*loc_1.loc[3]+ loc_1.loc[4]*loc_1.loc[4]);
        }
        else
        {
            location loc_2 = sites_by_ptr[site_2];

            double diff_x = loc_1.loc[0] - loc_2.loc[0];
            double diff_y = loc_1.loc[1] - loc_2.loc[1];
            double diff_z = loc_1.loc[2] - loc_2.loc[2];
            double diff_a = loc_1.loc[3] - loc_2.loc[3];
            double diff_b = loc_1.loc[4] - loc_2.loc[4];

            distance = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z + diff_a*diff_a + diff_b*diff_b);
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

// 6D infinite hypercubic lattice
class unbound_hypercubic_Lattice_6D: public Lattice
{
private:
    class location
    {
    public:
        int loc[6];

        location() {};

        location(int first, int second, int third, int fourth, int fifth, int sixth)
        {
            loc[0] = first;
            loc[1] = second;
            loc[2] = third;
            loc[3] = fourth;
            loc[4] = fifth;
            loc[5] = sixth;
        }
        bool operator< (const location & other) const
        {
            if(loc[0] < other.loc[0]) return true;
            if(other.loc[0] < loc[0]) return false;
            if(loc[1] < other.loc[1]) return true;
            if(other.loc[1] < loc[1]) return false;
            if(loc[2] < other.loc[2]) return true;
            if(other.loc[2] < loc[2]) return false;
            if(loc[3] < other.loc[3]) return true;
            if(other.loc[3] < loc[3]) return false;
            if(loc[4] < other.loc[4]) return true;
            if(other.loc[4] < loc[4]) return false;
            if(loc[5] < other.loc[5]) return true;
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
    ~unbound_hypercubic_Lattice_6D()
    {
        std::map<Site *, location>::iterator to_delete;

        for(to_delete = sites_by_ptr.begin(); to_delete != sites_by_ptr.end(); to_delete++)
        {
            delete to_delete->first;
        }
    }

    Site * get_origin(void)
    {
        location origin(0, 0, 0, 0, 0, 0);

        return get_site_ptr(origin);
    }

    void set_current_site(Site * site)
    {
        current_site = sites_by_ptr[site];
        next_neighbor = 0;
    }

    bool more_neighbors()
    {
        if(next_neighbor < 12) return true;

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
        case 6: // Out
            next_neighbor++;
            neighbor.loc[3] += 1;
            return get_site_ptr(neighbor, new_site);
        case 7: // In
            next_neighbor++;
            neighbor.loc[3] -= 1;
            return get_site_ptr(neighbor, new_site);
        case 8: // Top
            next_neighbor++;
            neighbor.loc[4] += 1;
            return get_site_ptr(neighbor, new_site);
        case 9: // Bottom
            next_neighbor++;
            neighbor.loc[4] -= 1;
            return get_site_ptr(neighbor, new_site);
        case 10: // Left
            next_neighbor++;
            neighbor.loc[5] += 1;
            return get_site_ptr(neighbor, new_site);
        case 11: // Right
            next_neighbor++;
            neighbor.loc[5] -= 1;
            return get_site_ptr(neighbor, new_site);
        }
    }

    double get_euclidian_distance(Site * site_1, Site * site_2 = NULL)
    {
        double distance;
        location loc_1 = sites_by_ptr[site_1];


        if (site_2 == NULL)
        {

            distance = sqrt(loc_1.loc[0]*loc_1.loc[0] + loc_1.loc[1]*loc_1.loc[1]
                    + loc_1.loc[2]*loc_1.loc[2] + loc_1.loc[3]*loc_1.loc[3]
                    + loc_1.loc[4]*loc_1.loc[4] + loc_1.loc[5]*loc_1.loc[5]);
        }
        else
        {
            location loc_2 = sites_by_ptr[site_2];

            double diff_x = loc_1.loc[0] - loc_2.loc[0];
            double diff_y = loc_1.loc[1] - loc_2.loc[1];
            double diff_z = loc_1.loc[2] - loc_2.loc[2];
            double diff_a = loc_1.loc[3] - loc_2.loc[3];
            double diff_b = loc_1.loc[4] - loc_2.loc[4];
            double diff_c = loc_1.loc[5] - loc_2.loc[5];

            distance = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z + diff_a*diff_a
                            + diff_b*diff_b + diff_c*diff_c);
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

