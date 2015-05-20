#include "Abstract_Classes.h"
#include "Simulation.h"

// Simple Site: Simplest type site
class simple_Site: public Site
{
private:
    bool occupied;
public:
    ~simple_Site()
    {
    }

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
