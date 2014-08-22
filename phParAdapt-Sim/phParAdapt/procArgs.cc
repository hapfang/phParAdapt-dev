#include <stdio.h>
#include <string.h>
#include "func.h"
#ifdef SIM
#include "SimPartitionedMesh.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

char options[255];
void procArgs(int argc, char *argv[])
{
  int c;

  /* first copy the entire argv into out dumpstring */
  if ( argc > 1 ){
     strcpy(options,argv[1]);
     for(c=2; c< argc; c++){
       strcat(options," ");
       strcat(options,argv[c]);
     }
   }

  /* look for command line arguments */
  while (--argc > 0 && (*++argv)[0] == '-')  {
    while ((argc > 0) && (c = *++argv[0]))  {
      switch (c)  {

      case 'h':
        if (PMU_rank() == 0)
         printf("\nYou have to specify input options in the following format in the file adapt.inp\n"
"\n=============================================================================================== \n"
"\nnumberErrorVars #Error Indicators {- usually 10}\n"

"\nrefWeights Weights for each Error Indicator type (default is 1 1 1 1 1 1 1 1 1 1)\n"

"\nglobalP Global Polynomial Order (default is 1)\n"

"\ntimeStepNumber #Time Step Number (default is 0)\n"

"\nnumelX #of elements in x direction (default is 0)\n" 

"\nNSFTAG -1 (default is  -1)\n" 

"\nensa_dof #of Solution variables to be write (default is 5, for the two-phase models -7)\n" 

"\nattributeFileName geom.spj (deafult for the preprocessing. For Adapt you have to put attribute file name without IC - geomNOIC.spj or other)\n"

"\nmeshFileName geom_p.sms (default name of the partitioned mesh directory. Change it, if different)\n"

"\nmodelFileName Parasolid Model File name \n"

"\nIdirection do sliced partitioning in specified direction, needs also numelX and BYPASS=1. By default is 0\n" 

"\nBYPASS 0 (by default)\n"

"\nnoldPhastaStyle 0 (using old PHASTA style. default is 0)\n" 

"\nzScale 0 prompt user for a x-y-z-coordinate scaling factor(default is 0)\n"   

"\nadaptFlag 0 ( Flag for adaptation. defalut(no adapt)=0, with adapt=1)\n"

"\nerrorName 0 not using  \n"

"\nSONFATH 0 (sets the number of homogenous directions for dynamic model averaging. Default is 0)\n"  

"\nlStart 0 (lin-2-quad. By default is 0)\n"   

"\nrRead #of Solution variables to be read (default is 5,for the two-phase models -7 ) Needs also rStart=1\n"    

"\nrStart 1 restart mode (read ICs from restart files) (default=0), \n"

"\nin case of using restart file and converting from density to pressure rStart=2)\n"

"\n in case of using linear portion of restart file rStart=3\n"

"\nin case of using quadratic portion of restart file rStart=4\n"

"\n AdaptStrategy #of Adapt Strategy mode\n"

"\n Strategy 1,2 anisotropic adaptation ()\n"

"\n Strategy 3,4 tag -driven isotropic adaptation (3 -read EI from error.step#.proc#, 4-from restart files)\n"

"\n Strategy 5,6 size-field driven isotropic adaptation (5 -read EI from error.step#.proc#, 6-from restart files)\n"

"\n Startegy 7 - octree splitting only for generating BIG meshes\n"

"\nAdaptFactor Value of adapt Factor (common value range is 0.6-0.7)\n"  

"\nAdaptOption 1 \n"  

"\nhmax #value of the Maximal Coarsening factor  \n"

"\nhmin 1 #value of the Maximal Refining factor \n"

"\n multipleRestarts 0 use files 'restartc.inp.ProcNumber+1' for initial conditions\n"

"\n (expects an integer argument for the number of procs)  only works together with the -rRead -option \n"

"\nPeriodic 0 (periodic case. Default is 0)\n"

"\ntiming 0 (print detailed timing statistics. Default is 0)\n"

"\nwGraph 0 (default is 0)\n"

"\nphastaVersion #PHASTA version\n"

"\nold_format 0 (flag for old format)\n"

"\nFortFormFlag 0 (default is 0)\n"

"\noutputFormat binary/ascii (default is binary)\n"

"\nCUBES 0 (Use multi-directional Inertial Partitioning, partition.in must be provided. Default is 0)\n"

"\ninternalBCNodes 0 (default is 0)\n"

"\nversion #version (not using)\n" 

"\nWRITEASC 0 (write in ASCII format. Default is 0)\n"

"\nphastaIO 1 (Use phastaIO libs for output file generation. Default is 0)\n"

"\nparted 1 (use already partitioned mesh-parted=1, parted=0)\n"  
          "\n\n");

        SimPartitionedMesh_stop();
        exit(1);
        break;

      default:
        printf(" Unknown option %c , Please use NSpre -h for the list of available options\n", c);
      }
    }
  }
}

#ifdef __cplusplus
}
#endif
