
#include <map>
#include <set>
#include <iostream>

#define INDEX_TYPE short int


// Abstract class for lattice of sites/bonds
class Lattice
{
public:
    virtual void
};

// Abstract class for sites
class Site
{
public:
    virtual bool is_occupied()=0;
    virtual void set_occupied_to()=0;
};

class simple_Site: public Site
{
private:
    bool occupied;

public:
    bool is_occupied() return occupied;
    void set_occupied_to(bool state) occupied=state return void;
};


// Abstract class for bonds
class Bond
{
public:
    virtual double get_strength()=0;
    virtual void set_strength_to()=0
};


class Square_Lattice
{
private:
    Site;
    Bonds;
};



// Abstract class for generating times to failure
class Timing
{
public:
    virtual double new_time_to_failure()=0;
};

// Abstract class for generating bond properties
class Strength
{
public:
    virtual double net_strength()=0;
};

int main(int argc, char **argv)
{
    return 0;
}
