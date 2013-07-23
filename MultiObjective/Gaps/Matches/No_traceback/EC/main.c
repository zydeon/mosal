#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <errno.h>

#define MAX_LENGTH			(4000)
#define NEG_INF				(INT_MIN)

void read_sequence(char seq[] , char *filename);
void init_dynamic_tables();
void init_base_cases();
void align( );
int max(int R, int S, int T);
void add_solution( int *last_match, int gaps );
void swap_tables( int *** t1, int *** t2 );
void print_table(int **T);
int max2(int n1, int n2);
int lcs_length();
void remove_dynamic_tables();


/* read sequences */
char seq1[MAX_LENGTH];
char seq2[MAX_LENGTH];


// 'X1' -> k-1 gaps,    'X2' -> k gaps
int ** R1, ** R2;		/* (ai , bj ) */
int ** S1, ** S2;		/* ('-', bj ) */
int ** T1, ** T2;		/* (ai , '-') */

int M, N, K;

FILE *f_out;

int main( int argc, char *argv[]) {

	if(argc != 3){
		printf("Usage: %s <seq1_file> <seq2_file>\n", argv[0]);
		return 0;
	}
	
	read_sequence(seq1, argv[1]);
	read_sequence(seq2, argv[2]);

	M = strlen(seq1); N = strlen(seq2);
	K = N+M-2*lcs_length();

	init_dynamic_tables();
	align();
	remove_dynamic_tables();

	return 0;
}

void read_sequence(char seq[], char *filename){
	char fasta_header[300];
	char line[300];

	FILE *f = fopen(filename, "r");
	if(f){
		if( fscanf(f, ">%[^\n]\n", fasta_header)>0 ){	// determine if file is valid
			while( fscanf(f, "%s", line)!=EOF )
				strcat(seq, line);			
		}
		else{
			printf("Sequence file format is not correct!\n");
			exit(-1);
		}
		fclose(f);
	}
	else{
		fprintf(stderr, "Failed to open sequence file: %s\n", strerror(errno));
		exit(-1);
	}
}

void init_dynamic_tables(){
	int i, j, match;

	R1 = (int **) malloc( (M+1)*sizeof(int *));
	R2 = (int **) malloc( (M+1)*sizeof(int *));

	T1 = (int **) malloc( (M+1)*sizeof(int *));
	T2 = (int **) malloc( (M+1)*sizeof(int *));

	S1 = (int **) malloc( (M+1)*sizeof(int *));
	S2 = (int **) malloc( (M+1)*sizeof(int *));

	for( i = 0; i <= M; ++i ) {
		R1[i] = (int *) malloc( (N+1)*sizeof(int) );
		R2[i] = (int *) malloc( (N+1)*sizeof(int) );

		T1[i] = (int *) malloc( (N+1)*sizeof(int) );
		T2[i] = (int *) malloc( (N+1)*sizeof(int) );

		S1[i] = (int *) malloc( (N+1)*sizeof(int) );
		S2[i] = (int *) malloc( (N+1)*sizeof(int) );
	}

	// initialize 3 tables for num_gaps = 0
	R1[0][0] = 0;
	T1[0][0] = NEG_INF;
	S1[0][0] = NEG_INF;

	// first columns
	for( i = 1; i <= M; ++i )
		R1[i][0] = S1[i][0] = T1[i][0] = NEG_INF;	
	
	// first lines
	for( j = 1; j <= N; ++j )
		R1[0][j] = S1[0][j] = T1[0][j] = NEG_INF;				
	
	for( i = 1; i <= M; ++i ){
		for( j = 1; j <= N; ++j ){
			if( i==j ){
				match = seq1[i-1] == seq2[j-1] ;
				R1[i][j] = R1[i-1][j-1] + match;
			}
			else
				R1[i][j] = NEG_INF;

			T1[i][j] = S1[i][j] = NEG_INF;	
		}
	}
}

void init_base_cases(int k){
	// 'k' - number of constant indels
	// 'k' >= 1
	int i;

	R2[0][0] = T2[0][0] = S2[0][0] = NEG_INF;

	if(k==1){
		for( i = 1; i < M+1; ++i ){
			S2[i][0] = NEG_INF;			// against table definition	(wrong position of gap)
			R2[i][0] = NEG_INF;			// out of recursive definition
			T2[i][0] = 0;			 
			
		}
		
		for( i = 1; i < N+1; ++i ){
			R2[0][i] = NEG_INF;			// against table definition	(wrong position of gap)
			S2[0][i] = NEG_INF;			// out of recursive definition
			T2[0][i] = 0;
		}
	}
	else {
		for( i = 1; i < M+1; ++i )
			R2[i][0] = S2[i][0] = T2[i][0] = NEG_INF;		
		
		for( i = 1; i < N+1; ++i )
			R2[0][i] = T2[0][i] = S2[0][i] = NEG_INF;
		
	}	
}

void align(){

	int k, i, j, match, last_match = NEG_INF, r, s, t;

	if( M==N ){			// possible solution with no indels
		printf("%d %d\n", R1[M][N], 0);
		last_match = R1[M][N];
	}


	for( k = 1; k <= K; ++k ){
		init_base_cases(k);

		for( i = 1; i <= M; ++i ){
			for( j = 1; j <= N; ++j ){

				match = seq1[i-1] == seq2[j-1];	

				// do not allow the value of 'match' to make it possible from NEG_INF
				r = R2[i-1][j-1] + match *(R2[i-1][j-1]!=NEG_INF);
				s = S2[i-1][j-1] + match *(S2[i-1][j-1]!=NEG_INF);
				t = T2[i-1][j-1] + match *(T2[i-1][j-1]!=NEG_INF);
				
				R2[i][j] = max( r, s, t );
				S2[i][j] = max( R1[i][j-1], S2[i][j-1], T1[i][j-1] );
				T2[i][j] = max( R1[i-1][j], S1[i-1][j], T2[i-1][j] );		
			}
		}
		
		// add to the solution set
		add_solution( &last_match , k);

		swap_tables( &R1, &R2 );
		swap_tables( &S1, &S2 );
		swap_tables( &T1, &T2 );
	}
	fclose(f_out);
}

int max( int n1, int n2, int n3 ){
	if(n1 >= n2)
		return n1 >= n3 ? n1 : n3;
	else
		return n2 >= n3 ? n2 : n3;
}

void add_solution( int *last_match, int gaps ){
	int solution;

	if(R2[M][N] >= T2[M][N]){
		if(R2[M][N] >= S2[M][N]){
			solution = R2[M][N];
		}
		else{
			solution = S2[M][N];			
		}
	}
	else{
		if( T2[M][N] >= S2[M][N] ){
			solution = T2[M][N];
		}
		else{
			solution = S2[M][N];			
		}
	}

	if( solution > *last_match ){
		printf("%d %d\n", solution, gaps);
		*last_match = solution;
	}
}

void swap_tables( int *** t1, int *** t2 ){
	int ** tmp = *t1;
	*t1 = *t2;
	*t2 = tmp;
}

int lcs_length(){

	int i, j, **L, result; 

	// init L
	L = (int **) malloc( (M+1)*sizeof(int *) );
	for( i = 0; i <= M; ++i ){
		L[i] = (int *) malloc( (N+1)*sizeof(int) );
		L[i][0] = 0;
	}
	for( j = 1; j <= N; ++j ){
		L[0][j] = 0;
	}	
	for( i = 1; i <= M; ++i ){
		for( j = 1; j <= N; ++j ){
			if( seq1[i-1] == seq2[j-1] )
				L[i][j] = L[i-1][j-1] + 1;
			else
				L[i][j] = max2(L[i-1][j] , L[i][j-1]);
		}
	}

	result = L[M][N];

	// remove table
	for( i = 0; i < M+1; ++i ){
		free(L[i]);
	}
	free(L);

	return result;
}


int max2(int n1, int n2){
	return n1 >= n2 ? n1 : n2;
}

void remove_dynamic_tables(){
	int i;

	for( i = 0; i < M+1; ++i ){
		free( R1[i] );
		free( R2[i] );
		free( T1[i] );
		free( T2[i] );
		free( S1[i] );
		free( S2[i] );
	}
	free( R1 );
	free( R2 );
	free( T1 );
	free( T2 );
	free( S1 );
	free( S2 );

}