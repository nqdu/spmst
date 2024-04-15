#ifndef _SPMST_PARALLEL_H
#define _SPMST_PARALLEL_H

void allocate_tasks(int ntasks,int nprocs,int myrank, int &start,int &end);
void print_progressbar(float percentage);

#endif