#ifndef ABSTRACT_CLASSES_H
#define ABSTRACT_CLASSES_H

#include <iostream>

// List of abstract classes
class Site;
class Bond;
class Lattice;
class Algorithm;
class Strength;
class Timing;

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
};

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
    virtual Bond * make_bond(Site * site_1, Site * site_2, double strength)=0;
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

    // Would like to eliminate
    virtual bool on_any_fault(Bond * bond)=0;
    virtual double get_fault_fraction(Bond * bond)=0;
    virtual Site * get_upstream(Site *)=0;

    // Developing
    virtual void modify_strength(Bond * bond)=0;
    virtual void setup_lattice(void)=0;
    virtual bool on_lattice(Site * site)=0;
};

// Abstract class for simulation
class Algorithm
{
public:
    Simulation * sim;
    virtual void initialize_sim(void)=0;
    virtual void set_sim(Simulation * to_set)=0;
    virtual void advance_sim(void)=0;
    virtual void free_up_memory(void)=0;
    virtual void write_sim_to_file()=0;
    virtual Bond * get_network_begin(void)=0;
    virtual Bond * get_network_next(void)=0;
    virtual bool more_network(void)=0;

    // Developing
    virtual void modify_strength(Bond * bond)=0;
    virtual bool check_growth(void)=0;
    virtual Bond * get_last_invaded(void)=0;
    virtual long int get_breakthrough_count(void)=0;
};

// Abstract class for generating bond strengths
class Strength
{
public:
    virtual double get_new_strength(void)=0;
};

// Abstract class for generating times to failure
class Timing
{
public:
    virtual double new_time_to_failure(void)=0;
};

#endif // ABSTRACT_CLASSES_H

/*
 * NOTES
 *
 * Sites are managed by Lattice.
 * Bonds are managed by Algorithm.
 *
*/
