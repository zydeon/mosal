#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <errno.h>

#define MAX_LENGTH			(4000)
#define LEN_ALPHABET		(26)


void read_sequence(char seq[] , char *filename);
void init_dynamic_tables();
void init_subs_table( char * filename );
void align();
int get_score( int i, int j );
void init_base_cases(int k );
int max(int n1, int n2, int n3);
void swap_dynamic_tables();
int lcs_length();
int max2(int n1, int n2);
void print_table();
void remove_dynamic_tables();

/* read sequences */
char seq1[MAX_LENGTH];
char seq2[MAX_LENGTH];

int ** P1, **P2;		// P1 -> k-1	  P2 -> k
int N, M, K;

int SS[LEN_ALPHABET][LEN_ALPHABET];	/* subsitution score table */

int NEG_INF;

int main(int argc, char *argv[]) {

	if(argc != 4){
		printf("Usage: %s <seq1_file> <seq2_file> <subs_file>\n", argv[0]);
		return 0;
	}

	read_sequence(seq1, argv[1]);
	read_sequence(seq2, argv[2]);

	M = strlen(seq1); N = strlen(seq2);
	K = N+M-2*lcs_length();

	init_subs_table( argv[3] );
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
	int i, j, score;

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

	// rest of the table	
	for( i = 1; i <= M; ++i ){
		for( j = 1; j <= N; ++j ){
			if( i==j ){
				score = get_score(i-1,j-1);
				P1[i][j] = P1[i-1][j-1] + score;
			}
			else
				P1[i][j] = NEG_INF;
		}
	}

}

void init_subs_table( char * filename ){
	/*
	*	

	builds substitution table from files that have a table with the format:
		      A1  A2 ... An
		  A1   -  -   -   -
		  A2   -  -   -   -
		  ...  -  -   -   -
		  An   -  -   -   -

		  where Ai is the ith letter of the alphabet and each in each cell is an
		  integer with the correspnding substitution score of the letters crossing

	*
	*/
		
	FILE *f;
	char line[200];
	char c;
	int i, j, val, min_penalty = INT_MAX;
	int len_alphabet = 0;
	char alphabet[LEN_ALPHABET];

	f = fopen(filename, "r");

	while( fscanf(f,"#%[^\n]\n", line) ); 	/* read comments */

	while( fscanf(f, "%c", &c) && c !='\n' ){
		if( c >= 'A' && c <= 'Z' ){
			alphabet[ len_alphabet++ ] = c;
		}
	}
	
	for (i = 0; i < len_alphabet; ++i){
		fscanf(f,"\n%c", &c);	

		for (j = 0; j < len_alphabet; ++j){
			fscanf(f,"%d", &val);
			SS[ alphabet[i]-'A' ][ alphabet[j]-'A' ] = val;									/* subs table */
		}
		
		if(val < min_penalty)
			min_penalty = val;
	}

	NEG_INF = INT_MIN - min_penalty;	// do not overflow limits of an integer when adding the subs_score
	
	fclose(f);
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

	int k, i, j, score, last_match = NEG_INF, d;

	//printf("k=%d\n", 0);
	if( M==N ){			// possible solution with no indels
		printf("%d %d\n", P1[M][N], 0);
		last_match = P1[M][N];		
	}


	for( k = 1; k <= K; ++k ){
		init_base_cases(k);

		//printf("k=%d\n", k);
		for( i = 1; i <= M; ++i ){
			for( j = 1; j <= N; ++j ){

				score = get_score(i-1,j-1);				

				// do not allow the value of 'match' to make it possible from NEG_INF
				d = P2[i-1][j-1] + score *(P2[i-1][j-1]!=NEG_INF);
				
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

int get_score( int i, int j ){
	return SS[ seq1[i]-'A' ][ seq2[j]-'A' ];
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
