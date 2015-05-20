#include "Common_Runs.h"
#include "Simulation.h"

int main(int argc, char **argv)
{
    Simulation sim;

    sim.configure_from_file();
    sim.run_sim();

    return 0;
}




