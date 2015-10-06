/***************************************************
SIDDHANT TULI			2012A7PS077P
SHUCHITA BANTHIA		2012A7PS011P
***************************************************/
#define hash_vector_size 10

struct _pair {
	int col;
	int nz;
};
typedef struct _pair pair;		//col_index, numnz 
										//numCols, maxnz
struct _pairF {
	int row;
	float val;
};
typedef struct _pairF pairF;

struct node {
	int col_index;
	struct node *next;
	int col_nz;
};
typedef struct node *NODE;

struct hashtable {
	NODE *buckets;
	NODE *tails;
	int *bucket_nz;
	int *col_counts;
	float *avg_vector;
	int **row_indices;
};
typedef struct hashtable *HASHTABLE;
