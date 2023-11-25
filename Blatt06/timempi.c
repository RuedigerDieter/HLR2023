#define _DEFAULT_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

int main(void) {

    struct timeval tv;
    time_t time;
    int micro_sec;
    char time_string[30];
    char output[80];
    char hostname[30];

    int proc_id, proc_num;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

    gettimeofday(&tv, NULL);
    gethostname(hostname, 30);

    time = tv.tv_sec;
    micro_sec = tv.tv_usec;
    
    // strftime(time_string, 30, "%Y-%m-%d %T", localtime(&time));
    // snprintf(output, 80, "%s : %s.%d", hostname, time_string, (int)micro_sec);

    // printf("%s\n", output);
    // printf("%d\n", (int)micro_sec);

    if(proc_id == proc_num - 1) 
    {
        char* proc_output = malloc(sizeof(char) * 80);
        int proc_time = 0;

        int us_min = micro_sec;
        int us_max = micro_sec;
        for (int i = 0; i < proc_num; i++)
        {
            MPI_Recv(proc_output, 80, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // TODO receive time

            for (int i = 0; i < proc_num; i++)
            {
                if (us_min > proc_time])
                {
                    us_min = proc_time;
                }
                else if (us_max < proc_time)
                {
                    us_max = proc_time;
                }
            }
            // TODO: print string
        }
       
        printf("Kleinster uS-Anteil: %d\n", us_min);
        printf("Größte Differenz: %d\n", us_max - us_min);

        free(output);

        MPI_Bcast(NULL,0,MPI_INT,proc_num - 1, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Send(output, 80, MPI_CHAR, proc_num - 1, 0, MPI_COMM_WORLD);
        // TODO send time as int to - 1
        MPI_Bcast(NULL,0,MPI_INT,proc_num - 1, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
