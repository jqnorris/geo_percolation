#include "Abstract_Classes.h"

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
    Bond * make_bond(Site * first, Site * second, double strength)
    {
        simple_Bond * temp = new simple_Bond;
        temp->first = first;
        temp->second = second;
        temp->set_strength_to(strength);
        return temp;
    }
};
