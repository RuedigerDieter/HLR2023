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
    
    strftime(time_string, 30, "%Y-%m-%d %T", localtime(&time));
    snprintf(output, 80, "%s : %s.%d", hostname, time_string, (int)micro_sec);

    if(proc_id != proc_num - 1) 
    {
        MPI_Send(output, 80, MPI_CHAR, proc_num - 1, 0, MPI_COMM_WORLD);
        MPI_Send(&micro_sec, 1, MPI_LONG, proc_num - 1, 1, MPI_COMM_WORLD);
        // printf("Sent Proc %d\n", proc_id);

        MPI_Bcast(NULL,0,MPI_INT,proc_num - 1, MPI_COMM_WORLD);
    }
    else
    {
        char proc_output[80];
        long proc_time = 0;

        long us_min = micro_sec;
        long us_max = micro_sec;
        if(proc_num == 1)
        {
            printf("%s\n", output);
            //printf("%d\n", (int)micro_sec);

            printf("Kleinster uS-Anteil: %ld\n", us_min);
            printf("Größte Differenz: %ld\n", us_max - us_min);

            return 0;
        }
        for (int i = 0; i < proc_num - 1; i++)
        {
            MPI_Recv(&proc_output, 80, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&proc_time, 1, MPI_LONG, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (us_min > proc_time)
            {
                us_min = proc_time;
            }
            else if (us_max < proc_time)
            {
                us_max = proc_time;
            }

            printf("[%d] %s\n", i, proc_output);
        }

        printf("Kleinster uS-Anteil: %ld\n", us_min);
        printf("Größte Differenz: %ld\n", us_max - us_min);

        MPI_Bcast(NULL,0,MPI_INT,proc_num - 1, MPI_COMM_WORLD);
    }
    printf("Beendet jetzt!\n");
    MPI_Finalize();
    return 0;
}
