#ifndef MPI_SUPPORT
#define MPI_SUPPORT

#include <mpi.h>
#include "simulation_configuration.h"
#include "simulation_support.h"

struct ProcConfig{
  int size;
  int rank;
  MPI_Comm comm;
  int neutron_start;
  int neutron_end;
  int neutron_num;
  // int reactor_x_start;
  // int reactor_x_end;
};

// struct BasicInformation{
//   unsigned long int current_total_active_neutrons;
//   unsigned long int total_num_fissions;
//   double total_mev;
//   double total_joules;
// };


void MPI_Initialize(struct simulation_configuration_struct*, struct ProcConfig*);

// void getBasicInformation(struct channel_struct **, struct simulation_configuration_struct*, struct ProcConfig*, struct BasicInformation*);

extern unsigned long int getTotalNumberActiveNeutrons(struct neutron_struct*, struct ProcConfig);

extern unsigned long int getTotalNumFissions(struct channel_struct**, struct simulation_configuration_struct*, struct ProcConfig proc_config);

extern double ***getAtomQuantities(struct channel_struct**, struct simulation_configuration_struct*, struct ProcConfig);
// extern void getAtomQuantities(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config, unsigned long int (*total_proc_quantities)[configuration->channels_x][configuration->channels_y][11]);

extern void free_matrix(unsigned long int***, struct channel_struct**, struct simulation_configuration_struct*);
#endif