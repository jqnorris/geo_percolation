 Currently this is a set of notes about general design of this program.
 
 
 1. Simulation
 The simulation is written using interfaces (abstract classes). This allows for the details of the simulation to be determined at run time. The subclasses of these interfaces are implemented using forward declaration. The class names are contained in a single header file, and the implementation details are contained in separate files for each interface. Current naming conventions include the name of the interface in the name of the subclass. The simulation files work together to create a network of connected bonds. In the simulation Sites are created and destroyed within the Lattice interface, and Bonds are created and destroyed within the Algorithm interface. The network must: 1. Be connected, and 2. Have well ordered bonds.
 
 2. Analysis/Statisitcs
 Originally the analysis was written to be performed on a single network. The simulation would run and grow a network. The simulation would then be passed to the Statistics class and various statistical properties of the network could be calculated, and written to file.
 
 3. Load Simulation Parameters
 Currently anytime I want to change one of the parameters, I have to recompile the program. I want to take advantage of the use of interfaces and allow the simulation details to be determined at runtime using a configuration file.
