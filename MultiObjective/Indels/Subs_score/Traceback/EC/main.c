#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <errno.h>

#define MAX_LENGTH			(4000)
#define LEN_ALPHABET		(26)


void read_sequence(char seq[] , char *filename);
void init_dynamic_tables();
void remove_dynamic_tables();
void init_traceback_tables();
void remove_traceback_tables();
void init_subs_table( char * filename );
void align();
int get_score( int i, int j );
void init_base_cases(int k );
int max(int n1, int n2, int n3, int i, int j, int k, int score);
void swap_dynamic_tables();
int lcs_length();
int max2(int n1, int n2);
void print_table();
void traceback( int k );
void reverse(char s[]);

/* read sequences */
char seq1[MAX_LENGTH];
char seq2[MAX_LENGTH];

int ** P1, **P2;		// P1 -> k-1	  P2 -> k
char *** T;				// traceback

int SS[LEN_ALPHABET][LEN_ALPHABET];		/* subsitution score table */

int N, M, K;

int NEG_INF;	//defined in init_subs_table()

int main(int argc, char *argv[]) {

	if(argc <= 3){
		printf("Usage: %s <seq1_file> <seq2_file> <subs_file>\n", argv[0]);
		return 0;
	}
	
	read_sequence(seq1, argv[1]);
	read_sequence(seq2, argv[2]);

	M = strlen(seq1); N = strlen(seq2);
	K = N+M-2*lcs_length();

	init_subs_table( argv[3] );
	init_dynamic_tables();
	init_traceback_tables();

	align();			// traceback inside align()

	remove_dynamic_tables();
	remove_traceback_tables();

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

void init_traceback_tables(){
	int i, j, k;

	T = (char ***) malloc( (M+1)*sizeof(char **) );

	for( i = 0; i <= M; ++i ){
		T[i] = (char **) malloc((N+1)*sizeof(char *));
	}
	for( i = 0; i < M+1; ++i ){
		for( j = 0; j < N+1; ++j ){
			T[i][j] = (char *) malloc((K+1)*sizeof(char));			
		}
	}

	// first cell
	for( k = 0; k < K+1; ++k )
		T[0][0][k] = '-';

	// first columns of each table for each value of k
	for( k = 1; k < K+1; ++k ){
		for( i = 0; i < M+1; ++i ){
			T[i][0][k] = 'u';
		}
	}
	// first line of each table for each value of k
	for( k = 1; k < K+1; ++k ){
		for( j = 0; j < N+1; ++j ){
			T[0][j][k] = 'l';
		}
	}	

	// first traceback tables with no indels
	for( i = 1; i < M+1; ++i ){
		for( j = 1; j < N+1; ++j ){
			if( i==j )	T[i][j][0] = 'd';				
			else		T[i][j][0] = '-';
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

void init_subs_table( char * filename ){
	/* le tabela de ficheiro com formato */
	FILE *f;
	char line[200];
	char c;
	int i, j, val, min_penalty = INT_MAX;
	int len_alphabet = 0;
	char alphabet[LEN_ALPHABET];

	f = fopen(filename, "r");

	while( fscanf(f,"#%[^\n]\n", line) ); 	/* ler comentarios */

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

	if( min_penalty < 0)
		NEG_INF = INT_MIN - min_penalty;	// nao ultrapassar menos infinito (ao adicionar subs_score) ou valor passa a positivo (valor maximo)
	else
		NEG_INF = INT_MIN;
	

	fclose(f);
}

void align(){

	int k, i, j, score, last_match = NEG_INF;
	FILE *f_out = fopen("out", "w");


	if( M==N ){			// possible solution with no indels
		fprintf(f_out, "%d %d\n", P1[M][N], 0);
		printf("matches=%d indels=%d\n", P1[M][N], 0);
		traceback(0);
		last_match = P1[M][N];		
	}


	for( k = 1; k <= K; ++k ){
		init_base_cases(k);

		for( i = 1; i <= M; ++i ){
			for( j = 1; j <= N; ++j ){
				score = get_score(i-1,j-1);
				P2[i][j] = max( P2[i-1][j-1], P1[i-1][j], P1[i][j-1], i, j, k, score);
			}
		}
		
		// add to the solution set
		if( P2[M][N] > last_match ){
			fprintf(f_out, "%d %d\n", P2[M][N], k);
			printf("matches=%d indels=%d\n", P2[M][N], k);
			traceback(k);
			last_match = P2[M][N];
		}

		swap_dynamic_tables();
	}
	fclose(f_out);
}

int get_score( int i, int j ){
	return SS[ seq1[i]-'A' ][ seq2[j]-'A' ];
}

void swap_dynamic_tables(){
	int ** tmp = P1;
	P1 = P2;
	P2 = tmp;
}

int max(int diagonal, int up, int left, int i, int j, int k, int score){

	// if 3 options are impossible do not allow the value of 'match' to make it possible
	if( diagonal == NEG_INF && up==NEG_INF && left == NEG_INF ){
		T[i][j][k] = '-';
		return NEG_INF;
	}

	diagonal += score;

	if(diagonal >= up){
		if(diagonal >= left){
			T[i][j][k] = 'd';
			return diagonal;
		}
		else{
			T[i][j][k] = 'l';
			return left;			
		}
	}
	else{
		if( up >= left ){
			T[i][j][k] = 'u';
			return up;
		}
		else{
			T[i][j][k] = 'l';
			return left;			
		}
	}
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

void traceback( int k ){
	int i, j, len1, len2;
	char res1[MAX_LENGTH], res2[MAX_LENGTH];

	len1 = len2 = 0;
	i = M; 	j = N;	

	while( i!=0 || j !=0 ){

		switch( T[i][j][k] ){			
			case 'u':
				k = k-1;			// previous table
				i--;
				res2[len2++] = '-';
				res1[len1++] = seq1[i];		
			break;
			case 'l':	
				k = k-1;			// previous table
				j--;
				res1[len1++] = '-';
				res2[len2++] = seq2[j];				
			break;
			case 'd':
				i--; j--;
				res1[len1++] = seq1[i];
				res2[len2++] = seq2[j];	
			break;			
		}		

	}

	res1[len1] = '\0';
	res2[len2] = '\0';		
	reverse(res1);
	reverse(res2);	
	printf("%s\n%s\n\n", res1, res2);	
		
}

void reverse(char s[]){
	int i, len;
	char tmp;
	len = strlen(s);

	for(i=0; i<len/2; i++){
		tmp = s[i];
		s[i] = s[len-i-1];
		s[len-i-1] = tmp;
	}
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

void remove_traceback_tables(){
	int i, j;

	for( i = 0; i < M+1; ++i ){
		for( j = 0; j < N+1; ++j )
			free( T[i][j] );	
		free( T[i] );
	}
	free( T );

}