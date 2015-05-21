#ifndef FORWARD_DECLARATIONS_H
#define FORWARD_DECLARATIONS_H

// Simulation Class
class Simulation;

// Abstract Classes
class Site;
class Bond;
class Lattice;
class Algorithm;
class Strength;
class Timing;

// Sites
class simple_Site;

// Bonds
class simple_Bond;

// Lattices
class Lattice_1D;
class unbound_square_Lattice_2D;
class unbound_cubic_Lattice_3D;
class unbound_hypercubic_Lattice_4D;
class unbound_hypercubic_Lattice_5D;
class unbound_hypercubic_Lattice_6D;
class unbound_square_Lattice_2D_with_faults;

// Algorithms
class ip_central_Algorithm;

// Timings


// Strengths
class uniform_Strength;
class testing_Strength;

#endif // FORWARD_DECLARATIONS_H
