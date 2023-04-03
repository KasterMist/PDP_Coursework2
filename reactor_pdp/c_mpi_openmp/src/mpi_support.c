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


/* Initialize the MPI configuration and send it to proc_config. 
 * It will contain the information of the MPI environment and each process's allocated neutrons length.
 */ 
void MPI_Initialize(struct simulation_configuration_struct* configuration, struct ProcConfig* proc_config){

  // struct ProcConfig proc_config;
  proc_config->comm = MPI_COMM_WORLD;
  MPI_Comm_size(proc_config->comm, &proc_config->size);
  MPI_Comm_rank(proc_config->comm, &proc_config->rank);
  
  // proc_config->neutron_start = configuration->max_neutrons / proc_config->size * proc_config->rank;
  // if(proc_config->rank != proc_config->size - 1){
  //   proc_config->neutron_end = proc_config->neutron_start + (configuration->max_neutrons / proc_config->size);
  // }else{
  //   proc_config->neutron_end = configuration->max_neutrons;
  // }
  // proc_config->neutron_num = proc_config->neutron_end - proc_config->neutron_start;

  proc_config->neutron_num = configuration->max_neutrons / proc_config->size;
  if(proc_config->rank == proc_config->size - 1){
    proc_config->neutron_num += configuration->max_neutrons % proc_config->size;
  }
  
}

/* This will replace the original function "getNumberActiveNeutrons" in main.c
 * Based on the previous function "getNumberActiveNeutrons", this function add MPI_Reduce to 
 * merge all the active neutrons together.
 */
unsigned long int getTotalNumberActiveNeutrons(struct neutron_struct *neutrons, struct ProcConfig proc_config) {
  unsigned long int activeNeutrons=0;

  #pragma omp parallel for default(none) firstprivate(proc_config, neutrons) reduction(+:activeNeutrons)
  for (unsigned long int i=0;i<proc_config.neutron_num;i++) {
    if (neutrons[i].active) activeNeutrons++;
  }

  MPI_Reduce(&activeNeutrons, &activeNeutrons, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, proc_config.comm);
  return activeNeutrons; 
}

/*
 * This function replaces the original function "getTotalNumberFissions" in the main.c.
 * Based on the previous function, this function add MPI_Reduce to merge the fission number together.
 */
extern unsigned long int getTotalNumFissions(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config){ 
  unsigned long int proc_num_fissions = 0;
  unsigned long int total_num_fissions = 0;
  #pragma omp parallel for shared(reactor_core) reduction(+:proc_num_fissions)
  for (int i = 0; i < configuration->channels_x; i++) {
    for (int j = 0; j < configuration->channels_y; j++) { 
      proc_num_fissions += reactor_core[i][j].contents.fuel_assembly.num_fissions; 
    }
  } 
  MPI_Reduce(&proc_num_fissions, &total_num_fissions, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, proc_config.comm);
  return total_num_fissions;
}


/*
 * This function is to get the atom quantities of the initial reactors.
 * It is called after running initialiseReactorCore function.
 */
extern double ***returnInitialAtomQuantities(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config){

  double ***initial_quantities;
  initial_quantities = (double ***)malloc(configuration->channels_x * sizeof(double **));
  #pragma omp parallel for default(none) shared(initial_quantities) firstprivate(configuration)
  for (int i = 0; i < configuration->channels_x; i++) {
    initial_quantities[i] = (double **)malloc(configuration->channels_y * sizeof(double *));
    for (int j = 0; j < configuration->channels_y; j++) {
      initial_quantities[i][j] = (double *)malloc(11 * sizeof(double));
      for (int k = 0; k < 11; k++) {
        initial_quantities[i][j][k] = 0.0;
      }
    }
  }

  #pragma omp parallel for default(none) firstprivate(configuration, reactor_core) shared(initial_quantities)
  for (int i = 0; i < configuration->channels_x; i++) {
    
    for (int j = 0; j < configuration->channels_y; j++) {  
      
      if (reactor_core[i][j].type == FUEL_ASSEMBLY) {
        // getAtomQuantities(&reactor_core[i][j].contents.fuel_assembly, configuration, proc_config, &initial_quantities[i][j]);
        for(int k = 0; k < reactor_core[i][j].contents.fuel_assembly.num_pellets; k++){
          for(int l = 0; l < 11; l++){
            initial_quantities[i][j][l] += reactor_core[i][j].contents.fuel_assembly.quantities[k][l]; 
          }
        }
      }
    }
  }
  return initial_quantities;
}

/*
 * This function aims to get all the atom quantities and return the data to a 3 dimensional array.
 * First two dimension determines the reactor core.
 * The third dimension determines the atom type.
 */
extern double ***returnTotalAtomQuantities(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config, double*** initial_quantities){
  unsigned long int a;
  size_t size = sizeof(a);
  // double minus;
  // unsigned long int total_num_fissions = 0;
  
  // double proc_quantities[configuration->channels_x][configuration->channels_y][11];
  // unsigned long int total_quantities[configuration->channels_x][configuration->channels_y][11];
  // unsigned long int ***total_proc_quantities;
  double ***total_proc_quantities;
  total_proc_quantities = (double ***)malloc(configuration->channels_x * sizeof(double **));
  for (int i = 0; i < configuration->channels_x; i++) {
    total_proc_quantities[i] = (double **)malloc(configuration->channels_y * sizeof(double *));
    for (int j = 0; j < configuration->channels_y; j++) {
      total_proc_quantities[i][j] = (double *)malloc(11 * sizeof(double));
      for (int k = 0; k < 11; k++) {
        total_proc_quantities[i][j][k] = 0.0;
        // proc_quantities[i][j][k] = 0.0;
      }
    }
  }
  
  for (int i = 0; i < configuration->channels_x; i++) {
    // total_proc_quantities[i] = (double **)malloc(configuration->channels_y * sizeof(double*));
    for (int j = 0; j < configuration->channels_y; j++) {  
      // total_proc_quantities[i][j] = (double *)malloc(11 * sizeof(double));
      if (reactor_core[i][j].type == FUEL_ASSEMBLY) {
        // getAtomQuantities(&reactor_core[i][j].contents.fuel_assembly, configuration, proc_config, &total_proc_quantities[i][j]);
        for(int k = 0; k < reactor_core[i][j].contents.fuel_assembly.num_pellets; k++){
          for(int l = 0; l < 11; l++){
            total_proc_quantities[i][j][l] += reactor_core[i][j].contents.fuel_assembly.quantities[k][l]; 
          }
        }
      }
      MPI_Reduce(total_proc_quantities[i][j], total_proc_quantities[i][j], 11, MPI_DOUBLE, MPI_SUM, 0, proc_config.comm);

      for(int m = 0; m < 11; m++){
        // minus = (proc_config.size - 1) * initial_quantities[i][j][m];
        // total_proc_quantities[i][j][m] = total_proc_quantities[i][j][m] - minus;
        total_proc_quantities[i][j][m] = total_proc_quantities[i][j][m] - (proc_config.size - 1) * initial_quantities[i][j][m];
      }
      // MPI_Barrier(proc_config.comm);

    }
  }
  
  return total_proc_quantities;
}

// extern void getAtomQuantities(struct fuel_assembly_struct * fuel_rod, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config, double** total_proc_quantities){
//   for(int k = 0; k < fuel_rod->num_pellets; k++){
//     for(int l = 0; l < 11; l++){
//       (*total_proc_quantities)[l] += fuel_rod->quantities[k][l]; 
//     }
//   }
// }


/*
 * After finishing the total atom quantities calculation, free the created array.
 */
extern void free_matrix(double ***matrix, struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration){
  for (int i = 0; i < configuration->channels_x; i++) {
    for (int j = 0; j < configuration->channels_y; j++) {
      free(matrix[i][j]);
    } 
    free(matrix[i]);
  }
  free(matrix);
}
