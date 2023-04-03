#ifndef MPI_SUPPORT
#define MPI_SUPPORT

#include <mpi.h>
#include "simulation_configuration.h"
#include "simulation_support.h"

struct ProcConfig{
  int size;
  int rank;
  MPI_Comm comm;
  int neutron_num;
  // int neutron_start;
  // int neutron_end;
  // int reactor_x_start;
  // int reactor_x_end;
};


void MPI_Initialize(struct simulation_configuration_struct*, struct ProcConfig*);

extern unsigned long int getTotalNumberActiveNeutrons(struct neutron_struct*, struct ProcConfig);

extern unsigned long int getTotalNumFissions(struct channel_struct**, struct simulation_configuration_struct*, struct ProcConfig proc_config);

extern double ***returnInitialAtomQuantities(struct channel_struct**, struct simulation_configuration_struct*, struct ProcConfig);

extern double ***returnTotalAtomQuantities(struct channel_struct**, struct simulation_configuration_struct*, struct ProcConfig, double***);

// extern void getAtomQuantities(struct fuel_assembly_struct*, struct simulation_configuration_struct*, struct ProcConfig, double**);

extern void free_matrix(double***, struct channel_struct**, struct simulation_configuration_struct*);

#endif