#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>


int*
init (int N, int rank)
{
	// TODO
	int* buf = (int*)malloc(sizeof(int) * N);

	srand(time(NULL) + rank * 100);

	for (int i = 0; i < N; i++)
	{
		// Do not modify "% 25"
		buf[i] = rand() % 25;
	}

	return buf;
}

int
circle (int* buf, int array_size, int world_size, int rank)
{
	int* tmp;
	int* firstbuf;
	int b_iterate = 0;

	int src = (rank == 0)? world_size : rank - 1;
	int dest = (rank == world_size)? 0 : rank + 1;

	tmp = malloc(sizeof(int) * array_size);
	firstbuf = (rank == 0)?  NULL : malloc(sizeof(int) * array_size);

	if (rank == 0)
	{
		printf("rank %d: Broadcast startet\n", rank);
		MPI_Bcast(buf, array_size, MPI_INT, 0, MPI_COMM_WORLD);
		printf("rank %d: Broadcast gesendet\n", rank);
	}
	else
	{
		MPI_Bcast(firstbuf, array_size, MPI_INT, 0, MPI_COMM_WORLD);
		printf("rank %d: Broadcast empfangen\n", rank);
	}

	printf("rank %d: Broadcast beendet\n", rank);

	while(!b_iterate || 1)
	{
		printf("rank %d: startet circle", rank);
		tmp = buf;
		

		if (rank == world_size - 1)
		{
			printf("rank %d: sendet", rank);
			MPI_Ssend(buf, array_size, MPI_INT, dest, 0, MPI_COMM_WORLD);
			printf("rank %d: empfängt", rank);
			MPI_Recv(buf, array_size, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else
		{
			printf("rank %d: empfängt", rank);
			MPI_Recv(buf, array_size, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("rank %d: sendet", rank);
			MPI_Ssend(tmp, array_size, MPI_INT, dest, 0, MPI_COMM_WORLD);
		}

		/**
		 * Check Termination
		*/
		if (rank == world_size - 1)
		{
			b_iterate = (buf[0] == firstbuf[0])? 1 : 0;
			MPI_Bcast(&b_iterate, 1, MPI_INT, world_size - 1, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Bcast(&b_iterate, 1, MPI_INT, world_size - 1, MPI_COMM_WORLD);
		}
	}

	printf("rank %d: beendet circle", rank);
	return 0;
}

int
main (int argc, char** argv)
{
	int N;
	int rank;
	int* buf;
  	int ret;
	int* buf_ret;

	int world_size;
	int proc_array_size;
	int rest;

	if (argc < 2)
	{
		printf("Arguments error!\nPlease specify a buffer size.\n");
		return EXIT_FAILURE;
	}

	// Array length
	N = atoi(argv[1]);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	if (N >= world_size)
	{
		proc_array_size = N / world_size;
	}
	else
	{
		proc_array_size = 1;
	}
	rest = N % world_size;
	buf_ret = (rank == 0)? (int*)malloc(sizeof(int) * proc_array_size) : NULL;

	// TODO: Prozesse müssen wissen, wie groß ihr Array ist, bzw. Reste verteilen

	buf = init(proc_array_size, rank);

	/*
	 *	Alle Prozesse haben die richtigen Elemente, die Prozesse mit "Lücken" kriegen illegale Werte,
	 *	um die Kommunikation sauberer zu gestalten
	 */
	if (rank >= rest)
	{
		buf[proc_array_size - 1] = 25;
	}


	if (rank == 0)
	{
		printf("\nBEFORE\n");
		
		printf("rank 0: ");
		for (int i = 0; i < proc_array_size; i++)
		{
			if (!(buf[i] == 25))
				{
					printf("%d ", buf[i]);
				}		
		}
		printf("\n");

		for (int i = 1; i < world_size; i++)
		{
			MPI_Recv(buf_ret, proc_array_size, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("rank %d: ", i);
			for (int j = 0; j < proc_array_size; j++)
			{
				if (!(buf_ret[j] == 25))
				{
					printf("%d ", buf_ret[j]);
				}
			}
			printf("\n");
		}
	}
	else
	{
		MPI_Ssend(buf, proc_array_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	ret = circle(buf, proc_array_size, world_size, rank);

	if (rank == 0)
	{
		printf("\nAFTER\n");

		printf("rank 0: ");
		for (int i = 0; i < proc_array_size; i++)
		{
			printf("%d ", buf[i]);
		}
		printf("\n");

		for (int i = 1; i < world_size; i++)
		{
			MPI_Recv(buf, proc_array_size, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("rank %d: ", rank);
			for (int j = 0; j < proc_array_size; j++)
			{
				printf("%d", buf[j]);
			}
			printf("\n");
		}
	}
	else
	{
		MPI_Ssend(buf, proc_array_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	free(buf);
	MPI_Finalize();

	return EXIT_SUCCESS;
}
