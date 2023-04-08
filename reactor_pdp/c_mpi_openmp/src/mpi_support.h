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
};


void MPI_Initialize(struct simulation_configuration_struct*, struct ProcConfig*);

extern unsigned long int getTotalNumberActiveNeutrons(struct neutron_struct*, struct ProcConfig);

extern unsigned long int getTotalNumFissions(struct channel_struct**, struct simulation_configuration_struct*, struct ProcConfig proc_config);

extern double ****returnInitialAtomQuantities(struct channel_struct**, struct simulation_configuration_struct*, struct ProcConfig);

extern void synchronize_reactor_core(struct channel_struct *, struct simulation_configuration_struct* , struct ProcConfig proc_config, double**);

extern double ***returnTotalAtomQuantities(struct channel_struct**, struct simulation_configuration_struct*, struct ProcConfig, double****);

// extern void getAtomQuantities(struct fuel_assembly_struct*, struct simulation_configuration_struct*, struct ProcConfig, double**);

extern double ****returnEmptyAtomQuantities(struct channel_struct**, struct simulation_configuration_struct*, struct ProcConfig);

extern void synchronize_quantities(struct channel_struct **, struct simulation_configuration_struct*, struct ProcConfig, double****, int, int);

extern void free3dMatrix(double***, struct channel_struct**, struct simulation_configuration_struct*);

extern void free4dMatrix(double ****, struct channel_struct **, struct simulation_configuration_struct*);

#endif