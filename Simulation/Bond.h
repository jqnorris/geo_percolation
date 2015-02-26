#include "Site.h"

#ifndef BOND_H
#define BOND_H

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
    virtual Bond * make_bond(Simulation * sim, Site * site_1, Site * site_2)=0;
    virtual ~Bond() {};
};


#endif // BOND_H
