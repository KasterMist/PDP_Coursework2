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
  
  for (int i = 0; i < configuration->channels_x; i++) {
    for (int j = 0; j < configuration->channels_y; j++) { 
      proc_num_fissions += reactor_core[i][j].contents.fuel_assembly.num_fissions; 
    }
  } 
  MPI_Reduce(&proc_num_fissions, &total_num_fissions, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, proc_config.comm);
  return total_num_fissions;
}

/*
 * This function is to create an empty 4-dimensional array which is used to store atom quantities.
 * The value of each atom quantities is set 0.
 */
extern double **** returnEmptyAtomQuantities(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config){
  
  double **** empty_quantities = (double ****)malloc(configuration->channels_x * sizeof(double ***));
  
  for (int i = 0; i < configuration->channels_x; i++) {
    empty_quantities[i] = (double ***)malloc(configuration->channels_y * sizeof(double **));
    for (int j = 0; j < configuration->channels_y; j++) {

      if (reactor_core[i][j].type == FUEL_ASSEMBLY) {

        empty_quantities[i][j] = (double **)malloc(reactor_core[i][j].contents.fuel_assembly.num_pellets * sizeof(double *));

        #pragma omp parallel for default(none) shared(empty_quantities) firstprivate(reactor_core, i, j)
        for (int k = 0; k < reactor_core[i][j].contents.fuel_assembly.num_pellets; k++) {
          empty_quantities[i][j][k] = (double *)malloc(NUM_CHEMICALS * sizeof(double)); 
          for(int l = 0; l < NUM_CHEMICALS; l++){
            empty_quantities[i][j][k][l] = 0.0;
          }
        }
      }
      
      
    }
  }
  return empty_quantities;
}


/*
 * This function is to get the atom quantities of the initial reactors.
 * It is called after running initialiseReactorCore function.
 */
extern double **** returnInitialAtomQuantities(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config){

  double **** initial_quantities = returnEmptyAtomQuantities(reactor_core, configuration, proc_config);

  
  for (int i = 0; i < configuration->channels_x; i++) {
    
    for (int j = 0; j < configuration->channels_y; j++) {  
      
      if (reactor_core[i][j].type == FUEL_ASSEMBLY) {
        #pragma omp parallel for default(none) firstprivate(reactor_core, i, j) shared(initial_quantities)
        for(int k = 0; k < reactor_core[i][j].contents.fuel_assembly.num_pellets; k++){
          for(int l = 0; l < 11; l++){
            initial_quantities[i][j][k][l] = reactor_core[i][j].contents.fuel_assembly.quantities[k][l]; 
          }
        }
      }
    }
  }
  return initial_quantities;
}


/*
 * This function aims to create a 3 dimensional array which stores the modified U235.
 * It is not used in the main function.
 */
extern double *** setInitialModifiedU235(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config){


  double *** modified_U235 = (double ***)malloc(configuration->channels_x * sizeof(double **));
  for (int i = 0; i < configuration->channels_x; i++) {
    modified_U235[i] = (double **)malloc(configuration->channels_y * sizeof(double *));
    for (int j = 0; j < configuration->channels_y; j++) {

      if (reactor_core[i][j].type == FUEL_ASSEMBLY) {

        modified_U235[i][j] = (double *)malloc(reactor_core[i][j].contents.fuel_assembly.num_pellets * sizeof(double));
        #pragma omp parallel for default(none) firstprivate(reactor_core, i, j) shared(modified_U235)
        for (int k = 0; k < reactor_core[i][j].contents.fuel_assembly.num_pellets; k++) {
          modified_U235[i][j][k] = 0.0;
        }
      }
    }
  }
  return modified_U235;
}

/*
 * This function is to synchronize U235 and update the other atoms.
 * The U235 will be synchronized to each process. Other atoms will be updated to the reactor core in each process.
 */
extern void synchronize_and_update(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config, double**** initial_quantities, double**** modified_quantities, int frequency, int current_step){
  
  if(current_step % frequency == 0){
    for (int i=0;i<configuration->channels_x;i++) {
      for (int j=0;j<configuration->channels_y;j++) {
        if(reactor_core[i][j].type == FUEL_ASSEMBLY){
          double modified_U235[reactor_core[i][j].contents.fuel_assembly.num_pellets];
          double total_modified_U235[reactor_core[i][j].contents.fuel_assembly.num_pellets];

          // Update the modified atoms to the reactor core
          #pragma omp parallel for default(none) firstprivate(i, j) shared(modified_U235, modified_quantities, reactor_core)
          for(int k = 0; k < reactor_core[i][j].contents.fuel_assembly.num_pellets; k++){
            for(int l = 0; l < 11; l++){
              if(l == Ba141 | l == Kr92 | l == Xe140 | l == Sr94 | l == Xe134 | l == Zr103){
                reactor_core[i][j].contents.fuel_assembly.quantities[k][l] += modified_quantities[i][j][k][l];
                modified_quantities[i][j][k][l] = 0.0;
                
              }else if(l == U235){
                modified_U235[k] =  modified_quantities[i][j][k][l];
                modified_quantities[i][j][k][l] = 0.0;
              }

            }
          }
          // Synchronize the U235 quantities to all of the processes.
          MPI_Allreduce(&modified_U235, &total_modified_U235, reactor_core[i][j].contents.fuel_assembly.num_pellets, MPI_DOUBLE, MPI_SUM, proc_config.comm);  
          #pragma omp parallel for default(none) firstprivate(modified_U235, i, j) shared(reactor_core)
          for(int k = 0; k < reactor_core[i][j].contents.fuel_assembly.num_pellets; k++){
            reactor_core[i][j].contents.fuel_assembly.quantities[k][U235] += modified_U235[k];
          }
        }
        
      }
    }
  }

}


/*
 * This function aims to get all the atom quantities and return the data to a 3 dimensional array.
 * First two dimension determines the reactor core.
 * The third dimension determines the atom type.
 */
extern double *** returnTotalAtomQuantities(struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration, struct ProcConfig proc_config, double**** initial_quantities){
  double modified;
  double initial_reactor_core_atoms[configuration->channels_x][configuration->channels_y][11];
  
  double *** total_proc_quantities;
  total_proc_quantities = (double ***)malloc(configuration->channels_x * sizeof(double **));
  for (int i = 0; i < configuration->channels_x; i++) {
    total_proc_quantities[i] = (double **)malloc(configuration->channels_y * sizeof(double *));
    for (int j = 0; j < configuration->channels_y; j++) {
      total_proc_quantities[i][j] = (double *)malloc(11 * sizeof(double));
      for (int k = 0; k < 11; k++) {
        total_proc_quantities[i][j][k] = 0.0;
        initial_reactor_core_atoms[i][j][k] = 0.0;
      }
    }
  }
  
  for (int i = 0; i < configuration->channels_x; i++) {
    
    for (int j = 0; j < configuration->channels_y; j++) {  
      
      if (reactor_core[i][j].type == FUEL_ASSEMBLY) {

        for(int k = 0; k < reactor_core[i][j].contents.fuel_assembly.num_pellets; k++){
          for(int l = 0; l < 11; l++){
            modified = reactor_core[i][j].contents.fuel_assembly.quantities[k][l] - initial_quantities[i][j][k][l];
            total_proc_quantities[i][j][l] += modified;

            initial_reactor_core_atoms[i][j][l] += initial_quantities[i][j][k][l]; 
          }
        }
        
      }
      MPI_Reduce(total_proc_quantities[i][j], total_proc_quantities[i][j], 11, MPI_DOUBLE, MPI_SUM, 0, proc_config.comm);
      
      for(int m = 0; m < 11; m++){
        total_proc_quantities[i][j][m] += initial_reactor_core_atoms[i][j][m];
      }
      
    }
  }
  
  return total_proc_quantities;
}



/*
 * After finishing the total atom quantities calculation, free the created 3d array.
 * It is not used in the main function
 */
extern void free3dMatrix(double *** matrix, struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration){
  for (int i = 0; i < configuration->channels_x; i++) {
    for (int j = 0; j < configuration->channels_y; j++) {
      free(matrix[i][j]);
    } 
    free(matrix[i]);
  }
  free(matrix);
}

/*
 * After finishing the total atom quantities calculation, free the created 4d array.
 * It is not used in the main function.
 */
extern void free4dMatrix(double **** matrix, struct channel_struct ** reactor_core, struct simulation_configuration_struct* configuration){
  for(int i = 0; i < configuration->channels_x; i++){
    for(int j = 0; j < configuration->channels_y; j++){
      if (reactor_core[i][j].type == FUEL_ASSEMBLY) {
        for (int k = 0; k < reactor_core[i][j].contents.fuel_assembly.num_pellets; k++) {
          free(matrix[i][j][k]);
        }
      }
      free(matrix[i][j]);
    }
    free(matrix[i]);
  }
  free(matrix);
}
