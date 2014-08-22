#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "phParAdapt.h"

int
phParAdaptProcSize() {

    char filename[128];
    sprintf(filename,"procsize_file_tmp.%d",PMU_rank());

    char sss[128];
    sprintf(sss,"ps -p %d -o vsz > %s",getpid(),filename);

    system(sss);

    char dummy[128];

    // read the file
    FILE* fin=fopen(filename,"r");
    fscanf(fin,"%s",dummy);
    int size;
    fscanf(fin,"%d",&size);
    fclose(fin);

    sprintf(sss,"rm %s",filename);
    system(sss);

    return size;
}
