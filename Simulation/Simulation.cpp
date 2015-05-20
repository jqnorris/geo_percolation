#include "Simulation.h"
#include "Tools.h"
#include "Algorithm.cpp"
#include "Lattice.cpp"
#include "Strength.cpp"
#include "Bond.cpp"
#include "Site.cpp"
#include <map>
#include <fstream>

void Simulation::configure_from_file(void)
{
    // Store all configuration parameters
    std::map<std::string, std::string> config_params;

    // Load lattice parameters from file
    std::ifstream from_file("sim_config.txt");

    // Load header
    std::getline(from_file, name);

    // Load parameters
    std::string temp;
    while(std::getline(from_file, temp))
    {
        std::cout << tools::split_string(temp, "=").first << "\t" << tools::split_string(temp, "=").second << "\n";
        config_params.insert(tools::split_string(temp, "="));
    }

    // Configure simulation
    if(config_params.count("size_of_run") > 0)
    {
        size_of_run = std::stoi(config_params["size_of_run"]);
    }
    else
    {
        size_of_run = 0;
    }

    if(config_params.count("number_of_runs") > 0)
    {
        number_of_runs = std::stoi(config_params["number_of_runs"]);
    }
    else
    {
        number_of_runs = 0;
    }

    if(config_params.count("p_c") > 0)
    {
        p_c = std::stof(config_params["number_of_runs"]);
    }
    else
    {
        p_c = 0.5;
    }

    if(config_params.count("algorithm") > 0)
    {
        if(config_params["algorithm"] == "ip_central")
        {
            Algorithm_ptr = new ip_central_Algorithm;
            Algorithm_ptr->set_sim(this);
        }
        else
        {
            std::cout << "WARNING! Unknown Algorithm." << std::endl;
            exit(1);
        }
    }
    else
    {
        Algorithm_ptr = new ip_central_Algorithm;
        Algorithm_ptr->set_sim(this);
    }

    if(config_params.count("lattice_type") > 0)
    {
        if(config_params["lattice_type"] == "unbound_square")
        {
            Lattice_ptr = new unbound_square_Lattice_2D;
            Lattice_ptr->setup_lattice();
            Lattice_ptr->sim = this;
        }
        else if(config_params["lattice_type"] == "unbound_cubic")
        {
            Lattice_ptr = new unbound_cubic_Lattice_3D;
            Lattice_ptr->setup_lattice();
            Lattice_ptr->sim = this;
        }
        else
        {
            std::cout << "WARNING! Lattice type unknown." << std::endl;
            exit(1);
        }
    }
    else if (config_params.count("dimensions") > 0)
    {
        if(config_params["dimensions"] == "1")
        {
            Lattice_ptr = new Lattice_1D;
            Lattice_ptr->setup_lattice();
            Lattice_ptr->sim = this;
        }
        else if(config_params["dimensions"] == "2")
        {
            Lattice_ptr = new unbound_square_Lattice_2D;
            Lattice_ptr->setup_lattice();
            Lattice_ptr->sim = this;
        }
        else if(config_params["dimensions"] == "3")
        {
            Lattice_ptr = new unbound_cubic_Lattice_3D;
            Lattice_ptr->setup_lattice();
            Lattice_ptr->sim = this;
        }
        else if(config_params["dimensions"] == "4")
        {
            Lattice_ptr = new unbound_hypercubic_Lattice_4D;
            Lattice_ptr->setup_lattice();
            Lattice_ptr->sim = this;
        }
        else if(config_params["dimensions"] == "5")
        {
            Lattice_ptr = new unbound_hypercubic_Lattice_5D;
            Lattice_ptr->setup_lattice();
            Lattice_ptr->sim = this;
        }
        else if(config_params["dimensions"] == "6")
        {
            Lattice_ptr = new unbound_hypercubic_Lattice_6D;
            Lattice_ptr->setup_lattice();
            Lattice_ptr->sim = this;
        }
        else
        {
            std::cout << "WARNING! Lattice type unknown." << std::endl;
            exit(1);
        }
    }
    else
    {
        Lattice_ptr = new unbound_square_Lattice_2D;
        Lattice_ptr->setup_lattice();
        Lattice_ptr->sim = this;
    }

    if(config_params.count("type_of_strength") > 0)
    {
           if(config_params["type_of_strength"] == "uniform")
           {
               Strength_ptr = new uniform_Strength;
           }
           else if(config_params["type_of_strength"] == "testing")
           {
               Strength_ptr = new testing_Strength;
           }
           else
           {
               std::cout << "WARNING! Type of Strength unknown." << std::endl;
               exit(1);
           }
    }
    else
    {
        Strength_ptr = new uniform_Strength;
    }

    Bond_ptr = new simple_Bond;
    Site_ptr = new simple_Site;

    return;
}

void Simulation::run_sim(void)
{
    Algorithm_ptr->initialize_sim();

    for(long int i=0; i<size_of_run; i++)
    {
        Algorithm_ptr->advance_sim();
    }

    Algorithm_ptr->write_sim_to_file();
    Algorithm_ptr->free_up_memory();

    return;
}

Simulation::~Simulation()
{
    delete Algorithm_ptr;
    delete Lattice_ptr;
    delete Site_ptr;
    delete Bond_ptr;
    delete Strength_ptr;
}
