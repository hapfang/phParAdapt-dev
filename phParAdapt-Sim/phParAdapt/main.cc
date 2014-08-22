/////////////////////////////////////////////////////////////////////////////////////////////////////
// phParAdapt main
// 
// parallel adaptive pre/post processor for phasta
// switch for adptation added
//
// J. Mueller 2005/2006
/////////////////////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <strstream>
#include <unistd.h>

#include "phParAdapt.h"
#ifdef SIM
#include "SimPartitionedMesh.h"
#include "SimAdvMeshing.h"
#endif
#include "func.h"
#include "ccfunc.h"


//  #if ( defined SGI )
//  #include <stream.h>
//  #endif

#if ( defined SIM_PARASOLID )
#include "SimParasolidKrnl.h"
#define MkPschema(x) "P_SCHEMA=" #x
#define XMkPschema(x) MkPschema(x)
#elif ( defined DISCRETE )
#include "SimDiscrete.h"
#endif

using std::strstream;


// catches the debugger
// to run debugger: setenv (=export) catchDebugger=1
void
catchDebugger() {
    static volatile int debuggerPresent =0;
    while (!debuggerPresent ); // assign debuggerPresent=1
}

int main(int argc, char *argv[])
{

// depending on cluster setup nodes sometimes don't get these --
// even if exported through job script
// Make sure PATHS are CORRECT !
//     putenv("PARASOLID=/usr/local/parasolid/16.0");
//     putenv("P_LIST=/usr/local/parasolid/16.0/lispdata");
//     putenv("P_SCHEMA=/usr/local/parasolid/16.0/schema");

//SimPartitioneMesh_start HAS to be the first call!!
    SimPartitionedMesh_start(&argc, &argv);
    Sim_readLicenseFile(NULL);

    char log_file[128];
    sprintf(log_file,"MeshSim.%i.log",PMU_gid(PMU_rank(),0)+1);
    Sim_logOn(log_file);

    MS_init();

    SimAdvMeshing_start();

    SimMeshing_start();

    SimModel_start();

    // read the model and mesh
#if  ( defined SIM_PARASOLID )
    if(SimParasolid_start(1) != 0 ){
        if(PMU_rank()==0){
            printf("\nerror in SimParasolid_start\n");
        }
        exit(0);
    }
#elif ( defined DISCRETE )
    // first initialize !!!
    SimDiscrete_start(0);
#endif

    
    // to run debugger: setenv (=export) catchDebugger=1
    // 
    if ( getenv( "catchDebugger" ) ) {

        int parent_pid = getpid();
        int gdb_child = fork();

        if( gdb_child == 0 ) {
     
            printf("Debugger Process initiating\n");
            strstream exec_string;

#if ( defined decalp )
            exec_string <<"xterm -e idb " 
                        << " -pid "<< parent_pid <<" "<< argv[0] << "\n";
#endif
#if ( defined LINUX )
            exec_string <<"xterm -e gdb " 
                        << " -pid "<< parent_pid <<" "<< argv[0] << "\n";
#endif
#if ( defined SUN4 )
            exec_string <<"xterm -e dbx " 
                        << " - "<< parent_pid <<" "<< argv[0] << "\n";
#endif
#if ( defined SGI )
            exec_string <<"xterm -e dbx " 
                        << " -p "<< parent_pid <<" "<< argv[0] << "\n";
#endif
            system( exec_string.str() );
            exit(0);
        }
        catchDebugger();
    }

//     // readLicenseFile(NULL);
//     MS_registerKey("sci001 parasolid 20060115 0 sAWcdPBlTDijvfHnrSi++A==");
//     MS_registerKey("sci001 adv 20060115 0 z+LrlbPSqcz5WKfgbLUXwA==");
//     MS_registerKey("sci001 surface 20060115 0 NAsdjvJQ4Uycb7w9Ep3sHQ==");
//     MS_registerKey("sci001 volume 20060115 0 YyhekWwXgOTGHgZurZHWzQ==");
//     MS_registerKey("sci001 adapt 20060115 0 FaIVzwkqy+6efRD6ISG9rQ==");
//     MS_registerKey("sci001 discrete 20060115 0 ldO5LJHunJSdHc/2J3j+aQ==");
//     MS_registerKey("sci001 pmesh 20060115 0 M4Hqwj4RThWeK8OoXrMfPg==");    
//     MS_registerKey("sci001 attributes 20060115 0 KQKXuy6FrB94b2iu44VS+g==");

#if defined(PARALLEL) && defined(DEBUG)
  int t;
  if (sscanf(argv[argc-1], "%d", &t) == 1)
    sleep(t);
#endif

    switchAdapt_Preproc(argc, argv);

#ifdef DISCRETE 
    SimDiscrete_stop(0);
#endif
            
#ifdef SIM_PARASOLID
    int sflag ;
    if(( sflag = SimParasolid_stop(1) ) != 0){
        if(PMU_rank()==0){
            printf("\nerror in SimParasolid_stop(1)\n"); 
            printf("\nencountered error %d\n",sflag);
        }   
        exit(0);
         
    }   
#endif

    SimModel_stop(); 

    SimMeshing_stop();

    SimPartitionedMesh_stop();

    SimAdvMeshing_stop();
    MS_exit();

    Sim_logOff();

    return 0;
}
