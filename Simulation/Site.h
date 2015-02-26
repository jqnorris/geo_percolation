#ifndef SITE_H
#define SITE_H

// Abstract class for sites
class Site
{
public:
    virtual bool is_occupied(void)const=0;
    virtual void set_occupied_to(bool state)=0;
    virtual Site * make_site(Simulation * sim)=0;
    virtual ~Site() {};
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

#endif // SITE_H
