#ifndef SITE_H
#define SITE_H

#include "Forward_Declarations.h"
#include "Abstract_Classes.h"

// Simplest Site
class simple_Site: public Site
{
private:
    bool occupied;
public:
    bool is_occupied() const;
    void set_occupied_to(bool state);
    Site * make_site(Simulation * sim);
    ~simple_Site();
};

#endif // SITE_H
