#include "Common_Runs.h"
#include "Simulation.h"
#include "Tools.h"
#include "Analysis.h"

int main(int argc, char **argv)
{
    tools::params params = tools::load_cofiguration_from_file();

    Simulation sim;
    sim.configure_sim(params);
    sim.run_sim();

    Statistics * stats = new Statistics(&sim, params);
    stats->calculate();

    return 0;
}




