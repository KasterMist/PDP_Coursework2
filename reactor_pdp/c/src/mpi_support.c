#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>

#include "mpi_support.h"
#include "simulation_configuration.h"
#include "simulation_support.h"



void MPI_Initialize(struct simulation_configuration_struct* configuration, struct ProcConfig* proc_config){

  // struct ProcConfig proc_config;
  proc_config->comm = MPI_COMM_WORLD;
  MPI_Comm_size(proc_config->comm, &proc_config->size);
  MPI_Comm_rank(proc_config->comm, &proc_config->rank);
  

  proc_config->neutron_start = configuration->max_neutrons / proc_config->size * proc_config->rank;
  if(proc_config->rank != proc_config->size - 1){
    proc_config->neutron_end = proc_config->neutron_start + (configuration->max_neutrons / proc_config->size);
  }else{
    proc_config->neutron_end = configuration->max_neutrons;
  }
  proc_config->neutron_num = proc_config->neutron_end - proc_config->neutron_start;
  
}

unsigned long int getTotalNumberActiveNeutrons(struct neutron_struct *neutrons, struct ProcConfig proc_config) {
  unsigned long int activeNeutrons=0;
  for (unsigned long int i=0;i<proc_config.neutron_num;i++) {
    if (neutrons[i].active) activeNeutrons++;
  }
  MPI_Reduce(&activeNeutrons, &activeNeutrons, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, proc_config.comm);
  return activeNeutrons; 
}

extern unsigned long int getTotalNumFissions(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config){ 
  unsigned long int proc_num_fissions = 0;
  unsigned long int total_num_fissions = 0;
  for (int i = 0; i < configuration->channels_x; i++) {
    for (int j = 0; j < configuration->channels_y; j++) { 
      proc_num_fissions += reactor_core[i][j].contents.fuel_assembly.num_fissions; 
    }
  } 
  MPI_Reduce(&proc_num_fissions, &total_num_fissions, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, proc_config.comm);
  return total_num_fissions;
}


extern double ***getAtomQuantities(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config){
  unsigned long int a;
  size_t size = sizeof(a);
  // unsigned long int total_num_fissions = 0;
  
  unsigned long int proc_quantities[configuration->channels_x][configuration->channels_y][11];
  // unsigned long int total_quantities[configuration->channels_x][configuration->channels_y][11];
  // unsigned long int ***total_proc_quantities;
  // total_proc_quantities = (unsigned long int ***)malloc(configuration->channels_x * sizeof(unsigned long int**));
  // for (int i = 0; i < configuration->channels_x; i++) {
  //   total_proc_quantities[i] = (unsigned long int **)malloc(configuration->channels_y * sizeof(unsigned long int*));
  //   for (int j = 0; j < configuration->channels_y; j++) { 
  //     total_proc_quantities[i][j] = (unsigned long int *)malloc(11 * sizeof(unsigned long int));
  //   }
  // }
  double ***total_proc_quantities;
    total_proc_quantities = (double ***)malloc(configuration->channels_x * sizeof(double **));
    for (int i = 0; i < configuration->channels_x; i++) {
        total_proc_quantities[i] = (double **)malloc(configuration->channels_y * sizeof(double *));
        for (int j = 0; j < configuration->channels_y; j++) {
            total_proc_quantities[i][j] = (double *)malloc(11 * sizeof(double));
            for (int k = 0; k < 11; k++) {
                total_proc_quantities[i][j][k] = 1;
            }
        }
    }

  for (int i = 0; i < configuration->channels_x; i++) {
    // total_proc_quantities[i] = (unsigned long int **)malloc(configuration->channels_y * size);
    for (int j = 0; j < configuration->channels_y; j++) {  
      // total_proc_quantities[i][j] = (unsigned long int *)malloc(11 * size);
      if (reactor_core[i][j].type == FUEL_ASSEMBLY) {
      for(int k = 0; k < reactor_core[i][j].contents.fuel_assembly.num_pellets; k++){
        for(int l = 0; l < 11; l++){
          total_proc_quantities[i][j][l] += reactor_core[i][j].contents.fuel_assembly.quantities[k][l]; 
        }
        // total_proc_quantities[i][j][0] += reactor_core[i][j].contents.fuel_assembly.quantities[k][0];
        // total_proc_quantities[i][j][1] += reactor_core[i][j].contents.fuel_assembly.quantities[k][1]; 
        // total_proc_quantities[i][j][2] += reactor_core[i][j].contents.fuel_assembly.quantities[k][2];
        // total_proc_quantities[i][j][3] += reactor_core[i][j].contents.fuel_assembly.quantities[k][3];
        // total_proc_quantities[i][j][4] += reactor_core[i][j].contents.fuel_assembly.quantities[k][4];
        // total_proc_quantities[i][j][5] += reactor_core[i][j].contents.fuel_assembly.quantities[k][5];
        // total_proc_quantities[i][j][6] += reactor_core[i][j].contents.fuel_assembly.quantities[k][6];
        // total_proc_quantities[i][j][7] += reactor_core[i][j].contents.fuel_assembly.quantities[k][7];
        // total_proc_quantities[i][j][8] += reactor_core[i][j].contents.fuel_assembly.quantities[k][8];
        // total_proc_quantities[i][j][9] += reactor_core[i][j].contents.fuel_assembly.quantities[k][9];
        // total_proc_quantities[i][j][10] += reactor_core[i][j].contents.fuel_assembly.quantities[k][10];
      }
    }
      MPI_Reduce(total_proc_quantities[i][j], total_proc_quantities[i][j], 11, MPI_DOUBLE, MPI_SUM, 0, proc_config.comm);
    }
  }
  if(proc_config.rank == 0){
    printf("!!!!!!!!!!!!%ld\n", total_proc_quantities[0][0][0]);
  } 
  return total_proc_quantities;
}


// extern void getAtomQuantities(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config, unsigned long int (*total_proc_quantities)[configuration->channels_x][configuration->channels_y][11]){
//   // unsigned long int a;
//   // size_t size = sizeof(a);
//   // unsigned long int total_num_fissions = 0;
//   unsigned long int proc_quantities[configuration->channels_x][configuration->channels_y][11];
//   // unsigned long int total_proc_quantities[configuration->channels_x][configuration->channels_y][11];
  
//   for (int i = 0; i < configuration->channels_x; i++) {
//     for (int j = 0; j < configuration->channels_y; j++) { 
//     // proc_num_fissions += reactor_core[i][j].contents.fuel_assembly.num_fissions; 
//     // MPI_Reduce(&reactor_core[i][j].contents.fuel_assembly.num_fissions, &reactor_core[i][j].contents.fuel_assembly.num_fissions, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, proc_config.comm);
//       for(int k = 0; k < reactor_core[i][j].contents.fuel_assembly.num_pellets; k++){
//         proc_quantities[i][j][0] = reactor_core[i][j].contents.fuel_assembly.quantities[k][U235];
//         proc_quantities[i][j][1] = reactor_core[i][j].contents.fuel_assembly.quantities[k][U238]; 
//         proc_quantities[i][j][2] = reactor_core[i][j].contents.fuel_assembly.quantities[k][Pu239];
//         proc_quantities[i][j][3] = reactor_core[i][j].contents.fuel_assembly.quantities[k][U236];
//         proc_quantities[i][j][4] = reactor_core[i][j].contents.fuel_assembly.quantities[k][Ba141];
//         proc_quantities[i][j][5] = reactor_core[i][j].contents.fuel_assembly.quantities[k][Kr92];
//         proc_quantities[i][j][6] = reactor_core[i][j].contents.fuel_assembly.quantities[k][Xe140];
//         proc_quantities[i][j][7] = reactor_core[i][j].contents.fuel_assembly.quantities[k][Sr94];
//         proc_quantities[i][j][8] = reactor_core[i][j].contents.fuel_assembly.quantities[k][Xe134];
//         proc_quantities[i][j][9] = reactor_core[i][j].contents.fuel_assembly.quantities[k][Zr103];
//         proc_quantities[i][j][10] = reactor_core[i][j].contents.fuel_assembly.quantities[k][Pu240];
//       }
//       MPI_Reduce(proc_quantities[i][j], total_proc_quantities[i][j], 11, MPI_UNSIGNED_LONG, MPI_SUM, 0, proc_config.comm);
//     }
//   }
// }

extern void free_matrix(unsigned long int ***matrix, struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration){
  for (int i = 0; i < configuration->channels_x; i++) {
    for (int j = 0; j < configuration->channels_y; j++) {
      free(matrix[i][j]);
    } 
    free(matrix[i]);
  }
  free(matrix);
}

// void getBasicInformation(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig* proc_config, struct BasicInformation* basic_information){
//   unsigned long int current_active_neutrons;
//   unsigned long int proc_num_fissions = 0;
//   current_active_neutrons = getNumberActiveNeutrons(proc_config);
//   for (int i = 0; i < configuration->channels_x; i++) {
//     for (int j = 0; j < configuration->channels_y; j++) { 
//       proc_num_fissions += reactor_core[i][j]->contents->fuel_assembly.num_fissions; 
//     }
//   } 
//   MPI_Reduce(&current_active_neutrons, basic_information->current_total_active_neutrons, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, proc_config.comm);
//   MPI_Reduce(&proc_num_fissions, basic_information->total_num_fissions, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, proc_config.comm);
  
//   basic_information->total_mev=getMeVFromFissions(basic_information->total_num_fissions);
//   basic_information->total_joules=getJoulesFromMeV(basic_information->total_mev);
// }