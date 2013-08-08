/// \file shy.cpp Main file. Its only responsibility is to
/// instantiate a simulation and call its main methods.
#include <iostream>
#ifdef HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI
#include <string>

#include "simulation.hpp"

int Usage()
{
  std::cerr << "Usage: shy [input file]\n";
  return 0;
}

int main(int argc, char* argv[]) {

  if (argc != 2)
    return Usage();

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif // HAVE_MPI

  char const * const filename_cstr = argv[1];
  std::string filename(filename_cstr);
  
  Simulation sim;

  sim.Load(filename);

  sim.Init();

  sim.Run();

  sim.Finalize();
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif // HAVE_MPI

  return 0;

}
