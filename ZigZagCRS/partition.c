//#include <iostream>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include <stdlib.h>

#define row_block_size 8
#define colVal_block_size 20
#define init_block_size 5

int main(int argc, char *argv[]) {

	clock_t startTime, endTime;
	startTime = clock();

	int p, my_rank;

	/* initialize MPI stuff */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);	

	int num_threads = atoi(argv[1]); //= argv - command line arg;
	int row_num, nz, col_num, i;
	int *partition = (int *)malloc(sizeof(int)*(p+1));
	int *local_row_array;
	int entries_per_proc;// = nz/p;
	int np, lower_index, ySize;
	float *Y, *X;


	FILE *fp;
	fp = fopen("crs48x48.txt", "r");
	fscanf(fp, "%d", &nz);
	while (fgetc(fp) != '\n');
	fscanf(fp, "%d", &row_num);
	while (fgetc(fp) != '\n');
	fscanf(fp, "%d", &col_num);
	while (fgetc(fp) != '\n');

	//%d,\n%d,\n", &nz, &row_num, &col_num);

	printf("%d => NZ = %d\n",my_rank, nz);
	entries_per_proc = nz/p;
	X = (float *)malloc(sizeof(float)*col_num);
	

	if (my_rank == 0) {
		int *row_arr;
		
		row_arr = (int *)malloc(sizeof(int)*row_num);

		for (i=0; i<row_num; i++) {
			fscanf(fp, "%d", &row_arr[i]);
			//printf("rownum[%d]: %d\n", i, row_arr[i]);
			while (fgetc(fp) != '\n');
		}

		partition[0] =0;

		#pragma omp parallel for num_threads(p-1) private(lower_index, i)
		for (i=1; i<p; ++i) {
			lower_index = entries_per_proc * i;
			partition[i] = bSearch(row_arr, row_num, lower_index);
			printf("searching for np: %d got %d\n", lower_index, partition[i]);
			//nz_arr[i] = row_arr[indices[i]];
		}
		partition[p] = row_num;
		
		printf("Partition array:\n");
		for(i=0; i<=p; i++)
		{
			printf("%d\n", partition[i]);
		}
  	//	printf("sending in bcat, recv\n");

  		//READ INPUT VECTOR X
  		FILE *fp2;
		fp2 = fopen("Xvector.txt", "r");
  		for (i=0; i<col_num; i++) {
			fscanf(fp2, "%f,", &X[i]);
		}
		fclose(fp2);
	}
	
	MPI_Bcast(partition, p+1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(X, col_num, MPI_FLOAT, 0, MPI_COMM_WORLD);

  	//printf("success in bcat, recv\n");

  	int num_rows_local;
  	num_rows_local = partition[my_rank+1] - partition[my_rank];
  	if (num_rows_local == 0) {
 		//can't exit!!
 		//should send and exit? 
 		exit(0); 		
   	}

  	int *col_arr;
  	float *val_arr;
  	int local_nz;

	local_row_array = (int *)malloc(sizeof(int) * num_rows_local);

	// SEEK TO THE CORRECT ROW POSITION
	fseek(fp, row_block_size*(partition[my_rank]) + init_block_size*3, SEEK_SET);
	//printf("fp for rank %d --  %d, jumped %d\n", my_rank, ftell(fp), row_block_size*(partition[my_rank]) + init_block_size*3);
	//printf("fp for rank %d --  curpos %d\n", my_rank, ftell(fp));
	for (i = 0; i < num_rows_local; i++) {
		fscanf(fp, "%d", &local_row_array[i]);
		//printf("local row arr[%d, %d]: %d  curpos %d\n", my_rank, i, local_row_array[i], ftell(fp));
		
		//printf("fp for rank %d --  %d, \n", my_rank, ftell(fp));
		while (fgetc(fp) != '\n');		
	}
	//printf("local: %d\n", local_row_array[0]);
	/*for (i=0; i<num_rows_local; i++) {
		printf("local row arr[%d, %d]: %d\n", my_rank, i, local_row_array[i]);
	}*/
	//local_nz = (nz_recv_arr[my_rank] - nz_recv_arr[my_rank-1]);
	
	int t2;
	if (my_rank != p-1) {
		fscanf(fp, "%d", &t2);
		//printf("rank %d read char %d\n",my_rank,t2);
		while (fgetc(fp) != '\n');
	}
	else {
		t2 = nz;
	}
	local_nz = t2 - local_row_array[0];

	col_arr = (int *)malloc(sizeof(int) * local_nz);
  	val_arr = (float *)malloc(sizeof(float) * local_nz);

	// SEEK TO CORRECT COLUMN AND VALUE
	fseek(fp, row_block_size*(row_num+1) + init_block_size*3 + colVal_block_size*(local_row_array[0]), SEEK_SET);

	//printf("col val array\n");
	for (i = 0; i < local_nz; i++) {
		fscanf(fp, "%d, %f", &col_arr[i], &val_arr[i]);
		while (fgetc(fp) != '\n');
		//printf("rank: %d, %d, %f,\n", my_rank,col_arr[i], val_arr[i]);
	}
	//printf("val check %d\n", val_arr[0]);
	int nz_per_thread = (local_nz) / num_threads;
	int *partition_threads = (int *)malloc(sizeof(int)*(num_threads+1));

	printf("\n********\nRank: %d, local nz: %d\n*******\n", my_rank, local_nz);
		

	partition_threads[0] = 0;
	#pragma omp parallel for num_threads(num_threads)
	for (i=1; i<num_threads; i++) {
		lower_index = local_row_array[0] + i*nz_per_thread;
		//printf("Thread searching for: %d\n", lower_index);

		partition_threads[i] = bSearch(local_row_array, num_rows_local, lower_index);
		//printf("part threads[%d] from proc:%d = %d\n", i, my_rank, partition_threads[i]);
	}
	partition_threads[num_threads] = num_rows_local;


	#pragma omp single
	{
		printf("Partiton for rank %d \n",my_rank);	
		for (i=0; i<=num_threads; i++) {
			printf("Thread[%d] = %d\n", i, partition_threads[i] + partition[my_rank]);
		}
	}
	/*
		Array partition_threads holds partition of indices over local_row_array
	*/

	//Allocate output array	and initialise to zero => use calloc()
	if(my_rank == 0){
		ySize = row_num;
	}
	else {
		ySize = num_rows_local;
	}
	Y = (float*)calloc(ySize, sizeof(float));

	//printf("Rank %d, Calloc done\n", my_rank);

	#pragma omp parallel num_threads(num_threads)
	{
		int id = omp_get_thread_num();
		int local_start = partition_threads[id];
		int local_end = partition_threads[id+1]-1;
		int thread_nz;
		if(local_start < local_end)
		{
			int local_index = local_start;
			int col_num, val, cnt, n;
			//printf("Rank %d, Thread %d, \n", my_rank, id);
			while(local_index <= local_end)
			{
				if(id == num_threads-1 && local_index == local_end){
					thread_nz = local_nz - local_row_array[local_index] + local_row_array[0];
				}
				else{
					thread_nz = local_row_array[local_index+1] - local_row_array[local_index];
				}
				
				//*******************************ZIG - ZAG - ZOOM *******************************

				if((local_index - local_start) % 2 == 0)
				{	
					//ZIG
					cnt = 0;
					while( cnt < thread_nz)
					{
						//non-zeroes for this thread start from number of nonzeroes processed by all other threads uptil now.
						n = local_row_array[local_index] - local_row_array[0]+cnt;
						//printf("\n***CHECK  %d\n", col_arr[n]);
						//X[col_arr[n]] += 2;
						Y[local_index] += val_arr[n] * X[col_arr[n]];
						cnt++;
					}
				}
				else
				{	
					//ZAG
					cnt = thread_nz-1;
					while(cnt >= 0)	
					{
						//non-zeroes for this thread start from number of nonzeroes processed by all other threads uptil now.
						n = local_row_array[local_index] - local_row_array[0]+cnt;
						//Y[local_index] += 2;
						Y[local_index] += val_arr[n] * X[col_arr[n]];
						cnt--;
					}	
				}

				local_index++;
			}
			
			/*#pragma omp critical
			{
				printf("\n\n*******");
				printf("Rank %d, Thread %d, Thread nz = %d\n", my_rank, id, thread_nz);
				for(i=0; i<local_end- local_start + 1; i++)
				{
				printf("Y[%d] = %f\n", i, Y[i]);
				}	
				printf("*******\n\n");
			}
*/
		}			
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	int base = partition[my_rank];
	
	/*printf("\n**********\nRank %d\n", my_rank);	
	for(i=0; i<num_rows_local; i++)
	{
		printf("Y[%d+%d] = %f\n", base,i, Y[i]);
	}	
*/
	//Concatenate output vector Y from all processes.
	if(my_rank!=0)
	{
		//All processes send
		MPI_Send(Y, num_rows_local, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
		//exit(0);
	}
	else{
		//Receive from all others
		int base, num_rows_to_receive; 
		for(i=1; i < p; i++){
			base = partition[i];
			num_rows_to_receive  = partition[i+1] - partition[i];
			if(num_rows_to_receive>0){
				MPI_Recv(&Y[base], num_rows_to_receive, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);			
			}
		}
	}

	if(my_rank == 0)
	{
		for(i=0; i<row_num; i++)
		{
			printf("Y[%d] = %f\n", i, Y[i]);
		}	
	}
		
	fclose(fp);
	endTime = clock();
	printf("\nNP = %d, Time taken: %ld\n",p,endTime-startTime);
	MPI_Finalize();
	return 0;
}

int bSearch(int arr[], int len, int val) {
	int i, start = 0, end = len-1, mid;
		
	while (start < end) {
		mid = (start+end)/2;	
		
		if (mid > 0 && arr[mid-1] <= val && val <= arr[mid]) {
			if(arr[mid]-val >= val-arr[mid-1]){
				return mid-1;
			}
			else{
				return mid;
			}
		}
		else if (mid < len-1 && arr[mid+1] >= val && val >= arr[mid]) {
			if(arr[mid+1]-val >= val-arr[mid]){
				return mid;
			}
			else{
				return mid+1;
			}
		}
		else if (val > arr[mid]) {
			start = mid + 1;
		}
		else if(val < arr[mid]) {
			end = mid - 1;
		}
	}
	
	if(start == end)
	{
		//printf("ERROR\n");
		return start;
		//exit(0);
	}	
}