#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include "hashtable.h"
#include <time.h>
#include <string.h>
#include <math.h>

#define col_block_size 8
#define rowVal_block_size 20
#define init_block_size 5
#define prime_arr_len 100
#define hash_vector_size 10
#define vector_block_size 6
#define NUM_SEEDS 10

void genHashWeightsArr(int hash_weights[], int number_of_rows, int prime_arr[]) {
	int r, i;

	for (i=0; i<number_of_rows; i++) {
		r = rand() % prime_arr_len;
		hash_weights[i] = prime_arr[r];
	}
}

HASHTABLE createHashtable(int num_partitions, int row_num) {
	HASHTABLE hashtbl = (HASHTABLE)malloc(sizeof(struct hashtable));
	//hashtbl->buckets = (NODE)malloc(sizeof(struct node)*50);
	//hashtbl->count = 0;
	//printf("allocated inside\n");
	int i;
	hashtbl->buckets = (NODE *)malloc(sizeof(NODE) * num_partitions);
	hashtbl->tails = (NODE *)malloc(sizeof(NODE) * num_partitions);
	hashtbl->col_counts = (int *)malloc(sizeof(int) * num_partitions);
	hashtbl->bucket_nz = (int *)malloc(sizeof(int) * num_partitions);
	hashtbl->avg_vector = (float *)calloc(hash_vector_size * num_partitions, sizeof(float));
	hashtbl->row_indices = (int **)malloc(sizeof(int *) * num_partitions);

	for (i=0; i<num_partitions; i++) {
		hashtbl->buckets[i] = NULL;
		hashtbl->tails[i] = NULL;
		hashtbl->col_counts[i] = 0;
		hashtbl->bucket_nz[i] = 0;
		hashtbl->row_indices[i] = (int *)calloc(row_num, sizeof(int));
	}
	
	return hashtbl;
}

HASHTABLE insert_list(HASHTABLE hashtbl, int insert_at_index, int hashVector[], int column_index, int row_indices[], int nz) {
	NODE cur = hashtbl->buckets[insert_at_index];
	//hashtbl->col_counts[insert_at_index]++;
	
	NODE newNode = (NODE)malloc(sizeof(struct node));
	newNode->col_index = column_index;
	newNode->col_nz = nz;
	newNode->next = NULL;
	
	if (cur == NULL) {
		int i;

		hashtbl->buckets[insert_at_index] = newNode;
		hashtbl->tails[insert_at_index] = newNode;
		//printf("inserted at bucket: %d, id: %d, val: %s\n", hashedVal, newNode->type, newNode->str);
		for (i=0; i<nz; i++) {
			hashtbl->row_indices[insert_at_index][row_indices[i]] += 1;
		}
		for (i=0; i<hash_vector_size; i++) {
			hashtbl->avg_vector[insert_at_index*hash_vector_size + i] = hashVector[i];
		}
		hashtbl->col_counts[insert_at_index] = 1;
		//printf("NULL: %d\n", insert_at_index);
		return hashtbl;
	}

	int i;
	for (i=0; i<nz; i++) {
		hashtbl->row_indices[insert_at_index][row_indices[i]] += 1;
	}
	
	//printf("not empty");
	/*
	while(cur->next != NULL) {
		cur = cur->next;
	}
	cur->next = newNode;
	*/
	float temp;
	int temp_index, cur_cnt =hashtbl->col_counts[insert_at_index];
	for (i=0; i<hash_vector_size; i++) {
		temp_index = insert_at_index * hash_vector_size + i;

		temp = cur_cnt * hashtbl->avg_vector[temp_index];
		hashtbl->avg_vector[temp_index] += hashVector[i];
		hashtbl->avg_vector[temp_index] /= cur_cnt+1;
	}
	hashtbl->col_counts[insert_at_index] += 1;

	hashtbl->tails[insert_at_index]->next = newNode;
	hashtbl->tails[insert_at_index] = newNode;
	//printf("inserted at bucket: %d, id: %d, val: %s\n", hashedVal, newNode->type, newNode->str);

	return hashtbl;
}

void hash(int column_nz_indices[], int column_nz_indices_num, int hash_weights[], int num_partitions, int row_num, int vector[hash_vector_size]) {
	int i, j, sum, wt;
	int batchSize = row_num/hash_vector_size;
	int c = 0;
	
	for (i=0; i<column_nz_indices_num; i++) {
		wt = hash_weights[column_nz_indices[i]];
		vector[column_nz_indices[i]/batchSize]+=wt;
	}
	/*
	for (i=0; i<column_nz_indices_num; i++) {
		while (c<column_nz_indices_num &&  column_nz_indices[c] < hash_vector_size*i) {
			sum += hash_weights[column_nz_indices[i*batchSize + j]];
			c+=1;
		vector[i] = sum;
	}*/
}

HASHTABLE insert_hash(HASHTABLE hashtbl, int column_nz_indices[], int column_nz_indices_num, int column_index, int num_partitions, int hash_weights[], int number_of_rows, int isSeedForP) {
	
	int hashedVector[hash_vector_size] = {0}; 
	hash(column_nz_indices, column_nz_indices_num, hash_weights, num_partitions, number_of_rows, hashedVector);
	//printf("inside insert_hash\n");
	int i, j;
	/*for (i=0; i<hash_vector_size; i++) {
		printf("hashed vector for col sthg %d\n", hashedVector[i]);
	}*/
	//printf("hashedVal: %d\n", hashedVal);
	float diff;
	float min_sum = 65534;
	int min_diff_index;

	if(isSeedForP == -1)
	{
		for (i=0; i<num_partitions; i++) {
			diff = 0;
			for (j=0; j<hash_vector_size; j++) {
				diff += abs(hashtbl->avg_vector[i*hash_vector_size + j] - hashedVector[j]);
			}

			//printf("diff with partition: %d = %f\n", i, diff);
			
			if (diff < min_sum) {
				min_sum = diff;
				min_diff_index = i;
			}
		}
		//printf("min diff with partition: %d = %f\n", min_diff_index, min_sum);

	}
	else
	{
		min_diff_index = isSeedForP;
	}

	hashtbl = insert_list(hashtbl, min_diff_index, hashedVector, column_index, column_nz_indices, column_nz_indices_num);
	
	return hashtbl;	
}

void print_hash(HASHTABLE hashtbl, int num_partitions) {
	int i;
	for(i=0; i<num_partitions; i++) {
		NODE cur = hashtbl->buckets[i];
		printf("hashtable row: %d, count: %d\n", i, hashtbl->col_counts[i]);
		//printf("%d:\t", i);
		while (cur != NULL) {
			printf("%d  ==>  ", cur->col_index);
			cur = cur->next;
		}
		printf("\n***\n");
	}
	printf("--------------------\n");
}

/*HASHTABLE delete_hash(HASHTABLE hashtbl, char str[]) {
	int hashValue = hash(str);
	NODE cur = hashtbl->buckets[hashValue];
	NODE prev = hashtbl->buckets[hashValue];
	if (strcmp(cur->str, str) == 0) {
		hashtbl->tails[hashValue] = hashtbl->buckets[hashValue]->next;	
		hashtbl->buckets[hashValue] = hashtbl->buckets[hashValue]->next;
		return hashtbl;
	}
	
	while (cur != NULL) {
		if (strcmp(cur->str, str) == 0) {
			if (cur->next != NULL) {
				hashtbl->tails[hashValue] = cur->next;
			}
			else {
				hashtbl->tails[hashValue] = prev;	
			}

			prev->next = cur->next;
			return hashtbl;
		}
		cur = cur->next;
	}
	return hashtbl;
}

int find_hash(HASHTABLE hashtbl, char str[]) {
	int hashValue = hash(str);
	//printf("hash of inside find_hash: %s = %d\n", str, hashValue);
	NODE cur = hashtbl->buckets[hashValue];
	while(cur != NULL) {
		if (strcmp(cur->str, str) == 0) {
			return cur->type;
		}
		cur = cur->next;
	}
	return -1;
}*/


int main(int argc, char *argv[]) {
	
	clock_t startTime, endTime;
	startTime = clock();

	int p, my_rank;

	/* initialize MPI stuff */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);	

	srand(time(NULL));
	
	int row_num, nz, col_num, i;

	FILE *fp;
	fp = fopen("crs48x48.txt", "r");
	fscanf(fp, "%d", &nz);
	while (fgetc(fp) != '\n');

	fscanf(fp, "%d", &row_num);
	while (fgetc(fp) != '\n');

	fscanf(fp, "%d", &col_num);
	while (fgetc(fp) != '\n');

	printf("%d => NZ = %d\n",my_rank, nz);

	FILE *fpseed;
	int seed[p];
	

	//int *column_partition = (int *)malloc(sizeof(int)*col_num);
	int *column_ptr;
	int *hash_weights;

	int num_cols_per_process[p];
	int *YPartition = (int*)calloc(row_num, sizeof(int));
	int *YmaxPartition = (int*)calloc(row_num, sizeof(int));
	
	const int nitems 		=  2;
    int blocklengths[2] 	= {1, 1};
    MPI_Datatype types[2] 	= {MPI_INT, MPI_INT};
    MPI_Datatype mpi_pair;
    MPI_Aint offsets[2];

    offsets[0] = offsetof(pair, col);
    offsets[1] = offsetof(pair, nz);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_pair);
    MPI_Type_commit(&mpi_pair);

    //printf("datatype created\n");
	pair A_partition_column[p];
	pair *A_partition[p];

	pair *my_columns;

	column_ptr = (int *)malloc(sizeof(int) * (col_num+1));
	// I need how many non-zeros in each column in the matrix data
	for (i=0; i <= col_num; i++) {
		fscanf(fp, "%d", &column_ptr[i]);
		while (fgetc(fp) != '\n');
	}
	//column_ptr[i] = nz;

	if (my_rank == 0) {	
		fpseed = fopen("seed48x48.txt", "r");
		for(i=0; i<p; i++)
		{
			fscanf(fpseed, "%d\n", &seed[i]);
			printf("seed[%d]: %d\n", i, seed[i]);
		}
		fclose(fpseed);

		int i;
		int prime_arr[prime_arr_len] = {  2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
											101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
											191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
											283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
											401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
											509, 521, 523, 541 
										 };

	
		hash_weights = (int *)malloc(sizeof(int)*row_num);
		genHashWeightsArr(hash_weights, row_num, prime_arr);

		/*for (i=0; i<row_num; i++) {
			printf("hashweights[%d]: %d\n",i,  hash_weights[i]);
		}*/

		int *current_column_rows = (int *)malloc(sizeof(int)*row_num);	
		//printf("check 1\n");
		HASHTABLE hash_columns;
		hash_columns = createHashtable(p, row_num);

		// read row_arr and insert in the hashtable for each column
		int j,c, flag;

		//insert seed cols
		for (c=0; c<p; c++) {
			j = seed[c];
			int nz_in_current_col = column_ptr[j+1] - column_ptr[j];
		
			fseek(fp, col_block_size*(col_num+1) + init_block_size*3 + rowVal_block_size*(column_ptr[j]), SEEK_SET);
			//printf("inserting\n");
			for (i=0; i<nz_in_current_col; i++) {
				fscanf(fp, "%d,", &current_column_rows[i]);
				while (fgetc(fp) != '\n');
			}
			
			hash_columns = insert_hash(hash_columns, current_column_rows, nz_in_current_col, j, p, hash_weights, row_num, c);
		}

		printf("\nSeeds:\n");
		print_hash(hash_columns, p);

		fseek(fp, col_block_size*(col_num+1) + init_block_size*3, SEEK_SET);

		//#pragma omp parallel for private(fp, j) num_threads(8)
		for (j=0; j<col_num; j++) {
			int nz_in_current_col = column_ptr[j+1] - column_ptr[j];
			flag =1;
			for (i=0; i<nz_in_current_col; i++) {
				fscanf(fp, "%d,", &current_column_rows[i]);
				while (fgetc(fp) != '\n');
			}
			//current_column_rows[i] = -1;
			/*if (j==0)
			{
				for (i=0; i<nz_in_current_col; i++) {
					printf("cur col[%d]: %d\n",i,  current_column_rows[i]);
				}
			*/	
			for(c =0 ; c<p; c++)
			{
				if(seed[c] == j)
				{
					flag = 0;
				}
			}
			
			if(flag == 1)
			{
				hash_columns = insert_hash(hash_columns, current_column_rows, nz_in_current_col, j, p, hash_weights, row_num, -1);
			}
				
			//}
		}
		// Load balancing
		//printf("inserted in hash\n");
		print_hash(hash_columns, p);

		// Generate a column-wise index storing the partition alloted to each column
		NODE temp;
		int max;

		#pragma omp parallel for num_threads(p)
		for (i=0; i<p; i++) {
			max = 0;
			A_partition_column[i].col = hash_columns->col_counts[i];
			A_partition[i] = (pair *)malloc(sizeof(pair)*A_partition_column[i].col);

			temp = hash_columns->buckets[i];
			for (j = 0; j < A_partition_column[i].col; j++) {
				A_partition[i][j].col = temp->col_index;
				A_partition[i][j].nz = temp->col_nz;
				if (temp->col_nz > max) {
					max = temp->col_nz;
				}
				temp = temp->next;
			}

			for (j=0; j<row_num; j++) {
				if(hash_columns->row_indices[i][j] > YmaxPartition[j]) {
					YmaxPartition[j] = hash_columns->row_indices[i][j];
					YPartition[j] = i;
				}
			}
			
			A_partition_column[i].nz = max;
		}		
	}

	// Broadcast the column-wise partition array
	MPI_Bcast(A_partition_column, p, mpi_pair, 0, MPI_COMM_WORLD);

	if (my_rank == 0) {
		my_columns = *A_partition;
	}
	else {
		my_columns = (pair *)malloc(sizeof(struct _pair)*A_partition_column[my_rank].col);
	}

	if (my_rank == 0) {
		for (i=1; i<p; i++) {
			MPI_Send(A_partition[i], A_partition_column[i].col, mpi_pair, i, 0, MPI_COMM_WORLD);
		}
	}
	else {
		MPI_Recv(my_columns, A_partition_column[my_rank].col, mpi_pair, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);			
	}

	MPI_Bcast(YPartition, row_num, MPI_INT, 0, MPI_COMM_WORLD);
	

	//check what recvd in mycolumns
	/*for(i=0; i<A_partition_column[my_rank].col; i++)
	{
		printf("Rank %d , Col no: %d, myNz : %d\n", my_rank, my_columns[i].col, my_columns[i].nz);
	}*/

	//partition_fp = (FILE **)malloc(sizeof(FILE*)*p);
	FILE *my_output;
	char f_name[20];


	int colIndex, myNz, rowIndex, j;
	float val;
	char *buffer;
	float *Y;
	
	Y = (float*)calloc(row_num, sizeof(float));

	
	//Read X
	FILE *fp2;
	fp2 = fopen("Xvector_algo2.txt", "r");
	int *X;
	X = (int*)malloc(sizeof(int)*A_partition_column[my_rank].col);
	
	printf("Rank %d recvd %d columns\n", my_rank, A_partition_column[my_rank].col);

	#pragma omp parallel for private(fp2, colIndex, i)
	for(i=0; i<A_partition_column[my_rank].col; i++)
	{
		colIndex = my_columns[i].col;
		fseek(fp2, colIndex*vector_block_size, SEEK_SET);
		fscanf(fp2, "%d\n", &X[i]);	
	}
  	fclose(fp2);

  	/*for(i=0; i<A_partition_column[my_rank].col; i++)
	{
		printf("Rank %d ::, X[%d] = %d\n",my_rank, i, X[i]);	
	}*/

	//for each column in A_partition_column[my_rank]...
	//Read non zeroes and multiply (computing local Y)...
	#pragma omp parallel for private(fp, colIndex, myNz, rowIndex, val)
	for(i=0; i<A_partition_column[my_rank].col; i++)
	{
		//printf("proc: %d, Operating on col %d \n", my_rank, colIndex);
		colIndex = my_columns[i].col;
		myNz = my_columns[i].nz;
		
		
		//seek to non-zeroes corresponding to this column in file
		fseek(fp, col_block_size*(col_num+1) + init_block_size*3 + rowVal_block_size*(column_ptr[colIndex]), SEEK_SET);
		
		//fread(buffer, myNz*rowVal_block_size,1,fp);
		//for each non zero...
		for(j=0; j<myNz; j++)
		{
			fscanf(fp, "%d, %f", &rowIndex, &val);
			while (fgetc(fp) != '\n');
			if(rowIndex>=row_num)
			{
				//printf("\n\n***********ERROR %d\n\n\n\n", rowIndex);
			}
			#pragma omp atomic
			Y[rowIndex]+= X[i]*val;
		}
	}
	//printf("end of loop: %d\n", my_rank);

	pairF *sendOthers[p];
	int numRowsInPartition[p], part;
	//numRowsInPartition = (int*) malloc(sizeof(int)*p);
	
	#pragma omp parallel for 
	for(i=0; i<row_num; i++)
	{
		numRowsInPartition[i] = 0;
	}

/*	for(i=0; i<row_num; i++)
	{
		printf("YPartition[%d] = %d\n", i, YPartition[i]);
	}
*/
	#pragma omp parallel for
	for(i=0; i<row_num; i++)
	{
		part = YPartition[i];
		#pragma omp atomic
		numRowsInPartition[part]++;
	}
	
	if (my_rank == 0) {
		for(i=0;i < p; i++)
		{
			printf("Rank %d got %d rows of Y vector\n", i, numRowsInPartition[i]);
		}
	}

	//make the arrays that have to be sent to other processes. 
	//pair arrays that store rowIndex and val.
	//allocate!
	for(i=0; i<p;i++)
	{
		//if(i!=my_rank)
		//{
			sendOthers[i] = (pairF*)malloc(sizeof(pairF)*numRowsInPartition[i]);
		//}
	}
	
	int *current = (int*) calloc(p, sizeof(int));
	int other, other_pos;

	//populate!
	for(i=0; i<row_num; i++)
	{
		other = YPartition[i];
		//if(other!=my_rank)
		//{
			other_pos = current[other];
			sendOthers[other][other_pos].row = i;
			sendOthers[other][other_pos].val = Y[i];
			current[other]++;
		//}
	}

	//write to respective files
	FILE *partition_fp[p];
	//open output files
	for (i=0; i< p; i++)
	{
		sprintf(f_name, "%d", i);
		//printf("open file %d\n", i);
		partition_fp[i] = fopen(f_name, "a");
	}

	//FILE *fp21 = fopen("hehe.txt", "a");
	for(i=0; i<p; i++)
	{
		if(i!=my_rank)
		{	
			other = i;
			for(j=0; j< numRowsInPartition[other]; j++)
			{
				if(sendOthers[other][j].val!=0){
					fprintf(partition_fp[other], "%d, %f, process %d\n",sendOthers[other][j].row, sendOthers[other][j].val, my_rank);
				}
			}
		}
	}

	//read from respective files and add! 
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (i = 0; i < p; ++i)
	{
		fclose(partition_fp[i]);
	}
	//printf("all files closed by rank %d\n", my_rank);
	sprintf(f_name, "%d", my_rank);
	partition_fp[my_rank] = fopen(f_name, "r");
	strcat(f_name, "_output.txt");
	my_output = fopen(f_name, "w");

	while (fscanf(partition_fp[my_rank], "%d, %f", &rowIndex, &val) > 0) // expect 1 successful conversion
	{
		while (fgetc(partition_fp[my_rank]) != '\n');
		//update local y
		//printf("\n****\nRank %d read value %f\n\n", my_rank, val);
		Y[rowIndex]+=val;
	}

	for(i=0; i<numRowsInPartition[my_rank]; i++)
	{
		rowIndex = sendOthers[my_rank][i].row;
		sendOthers[my_rank][i].val = Y[rowIndex];
		//these are the final values!
		fprintf(my_output, "%d, %f, process %d\n",sendOthers[my_rank][i].row, sendOthers[my_rank][i].val, my_rank);
	}

	fclose(partition_fp[my_rank]);
	fclose(my_output);
	
	endTime = clock();

	printf("\nrank = %d, Time taken: %lf\n", my_rank, (double)(endTime - startTime)/CLOCKS_PER_SEC);
	MPI_Finalize();
	return 0;

}