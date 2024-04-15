#include <iostream>

void allocate_tasks(int ntasks,int nprocs,int myrank, int &start,int &end)
{
    // allocate jobs to each rank
    int sub_n = ntasks / nprocs;
    int num_larger_procs = ntasks - nprocs * sub_n;
    if (myrank < num_larger_procs){ sub_n = sub_n + 1;
        start = 0 + myrank * sub_n;
    }
    else if (sub_n > 0){ 
        start = 0 + num_larger_procs + myrank * sub_n;
    }
    else { // this process has only zero elements
        start = -1;
        sub_n = 0;
    }
    end = start + sub_n - 1;
}

/**
 * print progress bar to the screen
 * @param percentage percentage of current progress
*/
void print_progressbar(float percentage) 
{
    // Ensure percentage is within [0.0, 100.0]
    percentage = (percentage < 0.0) ? 0.0 : (percentage > 100.0) ? 100.0 : percentage;

    // Determine the number of characters to represent the progress bar
    int numChars = (int)(percentage / 2.0);

    // Print the progress bar
    printf("[");
    for (int i = 0; i < 50; ++i) {
        if (i < numChars) {
            printf("=");
        } else if (i == numChars) {
            printf(">");
        } else {
            printf(" ");
        }
    }
    printf("] %.1f%%\r", percentage);  // Use carriage return to overwrite the line

    // Flush the output to ensure immediate display
    fflush(stdout);
}