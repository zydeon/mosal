#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <errno.h>

#define MAX_LENGTH		(4000)
#define NEG_INF			(INT_MIN)

void read_sequence(char seq[] , char *filename);
void init_dynamic_tables();
void align();
void init_base_cases(int k );
int max(int n1, int n2, int n3);
void swap_dynamic_tables();
int lcs_length();
int max2(int n1, int n2);
void remove_dynamic_tables();

/* read sequences */
char seq1[MAX_LENGTH];
char seq2[MAX_LENGTH];

int ** P1, **P2;		// P1 -> k-1	  P2 -> k
int N, M, K;

int main(int argc, char *argv[]) {

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

	P1 = (int **) malloc( (M+1)*sizeof(int *));
	P2 = (int **) malloc( (M+1)*sizeof(int *));

	for( i = 0; i <= M; ++i ) {
		P1[i] = (int *) malloc( (N+1)*sizeof(int) );
		P2[i] = (int *) malloc( (N+1)*sizeof(int) );
	}

	// initialize P1 as num_indels=0, traceback is only in diagonal
	P1[0][0] = 0;
	// first column
	for( i = 1; i <= M; ++i ){
		P1[i][0] = NEG_INF;
	}

	// first line
	for( j = 1; j <= N; ++j ){
		P1[0][j] = NEG_INF;
	}	
	for( i = 1; i <= M; ++i ){
		for( j = 1; j <= N; ++j ){
			if( i==j ){
				match = seq1[i-1] == seq2[j-1];
				P1[i][j] = P1[i-1][j-1] + match;
			}
			else
				P1[i][j] = NEG_INF;
		}
	}

}

void init_base_cases(int k){
	// 'k' - number of constant indels

	int i;
	// first column
	for( i = 0; i <= M; ++i ){
		if(i==k)
			P2[i][0] = 0;
		else
			P2[i][0] = NEG_INF;
	}

	// first line
	for( i = 0; i <= N; ++i ){
		if(i==k)
			P2[0][i] = 0;
		else
			P2[0][i] = NEG_INF;
	}	
}


void align(){

	int k, i, j, match, last_match = NEG_INF, d;

	//printf("k=%d\n", 0);
	if( M==N ){			// possible solution with no indels
		printf("%d %d\n", P1[M][N], 0);
		last_match = P1[M][N];		
	}


	for( k = 1; k <= K; ++k ){
		init_base_cases(k);

		for( i = 1; i <= M; ++i ){
			for( j = 1; j <= N; ++j ){

				match = seq1[i-1] == seq2[j-1];

				// do not allow the value of 'match' to make it possible from NEG_INF
				d = P2[i-1][j-1] + match *(P2[i-1][j-1]!=NEG_INF);
				
				P2[i][j] = max( d, P1[i-1][j], P1[i][j-1] );

			}
		}
		
		// add to the solution set
		if( P2[M][N] > last_match ){
			printf("%d %d\n", P2[M][N], k);
			last_match = P2[M][N];
		}

		swap_dynamic_tables();
	}
}


void swap_dynamic_tables(){
	int ** tmp = P1;
	P1 = P2;
	P2 = tmp;
}


int max( int n1, int n2, int n3 ){
	if(n1 >= n2)
		return n1 >= n3 ? n1 : n3;
	else
		return n2 >= n3 ? n2 : n3;
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
		free( P1[i] );
		free( P2[i] );
	}
	free( P1 );
	free( P2 );

}
