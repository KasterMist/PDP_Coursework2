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


void MPI_Initialize(struct simulation_configuration_struct *, struct ProcConfig *);

extern unsigned long int getTotalNumberActiveNeutrons(struct neutron_struct *, struct ProcConfig);

extern unsigned long int getTotalNumFissions(struct channel_struct **, struct simulation_configuration_struct *, struct ProcConfig proc_config);

extern double **** returnEmptyAtomQuantities(struct channel_struct **, struct simulation_configuration_struct *, struct ProcConfig);

extern double **** returnInitialAtomQuantities(struct channel_struct **, struct simulation_configuration_struct *, struct ProcConfig);

extern double *** setInitialModifiedU235(struct channel_struct **, struct simulation_configuration_struct *, struct ProcConfig);

extern double *** returnTotalAtomQuantities(struct channel_struct **, struct simulation_configuration_struct *, struct ProcConfig, double ****);

extern void synchronize_and_update(struct channel_struct **, struct simulation_configuration_struct *, struct ProcConfig, double ****, double ****, int, int);

extern void free3dMatrix(double ***, struct channel_struct **, struct simulation_configuration_struct *);

extern void free4dMatrix(double ****, struct channel_struct **, struct simulation_configuration_struct *);

#endif