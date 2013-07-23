#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <errno.h>

#define MAX_LENGTH			(4000)
#define NEG_INF				(INT_MIN)

void read_sequence(char seq[] , char *filename);
void init_dynamic_tables();
void remove_dynamic_tables();
void init_traceback_tables();
void remove_traceback_tables();
void init_base_cases();
void align( );
int maxR(int R , int T, int S, int i, int j, int k, int match);
int maxS(int R , int T, int S, int i, int j, int k);
int maxT(int R , int T, int S, int i, int j, int k);
void add_solution( int *last_match, int gaps );
void swap_tables( int *** t1, int *** t2 );
void print_table(int **T);
void traceback( char *** matrix, int k );
void reverse(char s[]);
int max2(int n1, int n2);
int lcs_length();

/* read sequences */
char seq1[MAX_LENGTH];
char seq2[MAX_LENGTH];


// 'X1' -> k-1 gaps,    'X2' -> k gaps
int ** R1, ** R2;		/* (ai , bj ) */
int ** S1, ** S2;		/* ('-', bj ) */
int ** T1, ** T2;		/* (ai , '-') */

//tracebacks
char *** TR, ***TT, ***TS;

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
	init_traceback_tables();
	
	align();		// traceback inside align()

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
	for( i = 1; i <= M; ++i ){
		R1[i][0] = S1[i][0] = T1[i][0] = NEG_INF;		
	}
	// first lines
	for( j = 1; j <= N; ++j ){
		R1[0][j] = S1[0][j] = T1[0][j] = NEG_INF;		
	}	

	// rest of the table
	for( i = 1; i <= M; ++i ){
		for( j = 1; j <= N; ++j ){
			if( i==j ){
				match = seq1[i-1] == seq2[j-1];
				R1[i][j] = R1[i-1][j-1] + match;
			}
			else
				R1[i][j] = NEG_INF;

			T1[i][j] = S1[i][j] = NEG_INF;
		}
	}
}

void init_traceback_tables(){
	int i, j, k;

	TR = (char ***) malloc( (M+1)*sizeof(char **) );
	TT = (char ***) malloc( (M+1)*sizeof(char **) );
	TS = (char ***) malloc( (M+1)*sizeof(char **) );

	for( i = 0; i <= M; ++i ){
		TR[i] = (char **) malloc((N+1)*sizeof(char *));
		TT[i] = (char **) malloc((N+1)*sizeof(char *));
		TS[i] = (char **) malloc((N+1)*sizeof(char *));
	}
	for( i = 0; i < M+1; ++i ){
		for( j = 0; j < N+1; ++j ){
			TR[i][j] = (char *) malloc((K+1)*sizeof(char));			
			TT[i][j] = (char *) malloc((K+1)*sizeof(char));			
			TS[i][j] = (char *) malloc((K+1)*sizeof(char));			
		}
	}

	// first columns of each table for each value of k
	for( k = 2; k < K+1; ++k ){
		for( i = 0; i < M+1; ++i ){
			TR[i][0][k] = '-';
			TT[i][0][k] = '-';
			TS[i][0][k] = '-';
		}
	}
	// first line of each table for each value of k
	for( k = 2; k < K+1; ++k ){
		for( j = 0; j < N+1; ++j ){
			TR[0][j][k] = '-';
			TT[0][j][k] = '-';
			TS[0][j][k] = '-';
		}
	}	

	// num_gaps = 1
	for( i = 1; i < M+1; ++i ){		// first column
		TR[i][0][1] = '-';
		TT[i][0][1] = '5';
	}
	for( j = 1; j < N+1; ++j ){		// first line
		TR[0][j][1] = '-';
		TS[0][j][1] = '9';
	}


	// first traceback tables with no gaps
	for( i = 1; i < M+1; ++i ){
		for( j = 1; j < N+1; ++j ){
			if( i==j )	TR[i][j][0] = '1';				
			else		TR[i][j][0] = '-';
			TT[i][j][0] = '-';
			TS[i][j][0] = '-';
		}
	}
	
}

void init_base_cases(int k){
	// 'k' - number of constant gaps
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
			T2[0][i] = NEG_INF;			// against table definition	(wrong position of gap)
			R2[0][i] = NEG_INF;			// out of recursive definition
			S2[0][i] = 0;
		}
	}
	else {
		for( i = 1; i < M+1; ++i )
			R2[i][0] = S2[i][0] = T2[i][0] = NEG_INF;		
		
		for( i = 1; i < N+1; ++i )
			R2[0][i] = T2[0][i] = S2[0][i] = NEG_INF;
		
	}		
}

FILE *fout;
void align(){

	int k, i, j, match, last_match = NEG_INF;
	fout = fopen("out", "w");

	if( M==N ){			// possible solution with no gaps
		fprintf(fout, "%d %d\n", R1[M][N], 0);
		printf("matches=%d gaps=%d\n", R1[M][N], 0);
		traceback( TR, 0 );
		last_match = R1[M][N];
	}


	for( k = 1; k <= K; ++k ){
		init_base_cases(k);

		for( i = 1; i <= M; ++i ){
			for( j = 1; j <= N; ++j ){

				match = seq1[i-1] == seq2[j-1];				

				R2[i][j] = maxR( R2[i-1][j-1], T2[i-1][j-1], S2[i-1][j-1], i,j,k, match);
				T2[i][j] = maxT( R1[i-1][j], T2[i-1][j], S1[i-1][j], i,j,k );
				S2[i][j] = maxS( R1[i][j-1], T1[i][j-1], S2[i][j-1], i,j,k );		

			}
		}
		
		// add to the solution set
		add_solution( &last_match , k);

		swap_tables( &R1, &R2 );
		swap_tables( &T1, &T2 );
		swap_tables( &S1, &S2 );
	}
	fclose(f_out);
}

int maxR(int R, int T, int S, int i, int j, int k, int match){

	// if 3 options are impossible do not allow the value of 'match' to make it possible
	if( R == NEG_INF && S == NEG_INF && T == NEG_INF ){
		TR[i][j][k] = '-';
		return NEG_INF;
	}

	R += match;			
	S += match;			
	T += match;			

	if(R >= T){
		if(R >= S){
			TR[i][j][k] = '1';
			return R;
		}
		else{
			TR[i][j][k] = '3';
			return S;			
		}
	}
	else{
		if( T >= S ){
			TR[i][j][k] = '2';
			return T;
		}
		else{
			TR[i][j][k] = '3';
			return S;			
		}
	}
}

int maxT(int R, int T, int S, int i, int j, int k){

	// if 3 options are impossible do not allow the value of 'match' to make it possible
	if( R == NEG_INF && S == NEG_INF && T == NEG_INF ){
		TT[i][j][k] = '-';
		return NEG_INF;
	}	

	if(R >= T){
		if(R >= S){
			TT[i][j][k] = '4';
			return R;
		}
		else{
			TT[i][j][k] = '6';
			return S;			
		}
	}
	else{
		if( T >= S ){
			TT[i][j][k] = '5';
			return T;
		}
		else{
			TT[i][j][k] = '6';
			return S;			
		}
	}
}

int maxS(int R, int T, int S, int i, int j, int k){

	// if 3 options are impossible do not allow the value of 'match' to make it possible
	if( R == NEG_INF && S == NEG_INF && T == NEG_INF ){
		TS[i][j][k] = '-';
		return NEG_INF;
	}	

	if(R >= T){
		if(R >= S){
			TS[i][j][k] = '7';
			return R;
		}
		else{
			TS[i][j][k] = '9';
			return S;			
		}
	}
	else{
		if( T >= S ){
			TS[i][j][k] = '8';
			return T;
		}
		else{
			TS[i][j][k] = '9';
			return S;			
		}
	}
}


void add_solution( int *last_match, int gaps ){
	int solution;
	char ***matrix;

	if(R2[M][N] >= T2[M][N]){
		if(R2[M][N] >= S2[M][N]){
			solution = R2[M][N];
			matrix = TR;
		}
		else{
			solution = S2[M][N];
			matrix = TS;			
		}
	}
	else{
		if( T2[M][N] >= S2[M][N] ){
			solution = T2[M][N];
			matrix = TT;
		}
		else{
			solution = S2[M][N];
			matrix = TS;			
		}
	}

	if( solution > *last_match ){
		fprintf(fout, "%d %d\n", solution, gaps);
		printf("matches=%d gaps=%d\n", solution, gaps);
		traceback( matrix, gaps );
		*last_match = solution;
	}
}

void swap_tables( int *** t1, int *** t2 ){
	int ** tmp = *t1;
	*t1 = *t2;
	*t2 = tmp;
}

void traceback( char *** matrix, int k ){
	int i, j, len1, len2;
	char res1[MAX_LENGTH], res2[MAX_LENGTH];

	len1 = len2 = 0;
	i = M; 	j = N;	

	while( i!=0 || j !=0 ){

		switch( matrix[i][j][k] ){

			// TR
			case '1':
				i--; j--;
				res1[len1++] = seq1[i];
				res2[len2++] = seq2[j];				
			break;
			case '2':
				i--; j--;
				res1[len1++] = seq1[i];
				res2[len2++] = seq2[j];
				matrix = TT;
			break;
			case '3':
				i--; j--;
				res1[len1++] = seq1[i];
				res2[len2++] = seq2[j];
				matrix = TS;
			break;			

			// TT
			case '4':
				k-=1;
				i--;
				res2[len2++] = '-';
				res1[len1++] = seq1[i];	
				matrix = TR;			
			break;
			case '5':
				i--;
				res2[len2++] = '-';
				res1[len1++] = seq1[i];					
			break;
			case '6':
				k-=1;
				i--;
				res2[len2++] = '-';
				res1[len1++] = seq1[i];	
				matrix = TS;
			break;	

			// TS
			case '7':
				k-=1;
				j--;
				res1[len1++] = '-';
				res2[len2++] = seq2[j];	
				matrix = TR;			
			break;
			case '8':
				k-=1;
				j--;
				res1[len1++] = '-';
				res2[len2++] = seq2[j];	
				matrix = TT;				
			break;
			case '9':
				j--;
				res1[len1++] = '-';
				res2[len2++] = seq2[j];		
			break;	

			case 'u':
				i--;
				res2[len2++] = '-';
				res1[len1++] = seq1[i];		
			break;
			case 'l':	
				j--;
				res1[len1++] = '-';
				res2[len2++] = seq2[j];				
			break;		

			default:
				if( matrix == TR )
					printf("ERROR in traceback on table R -> i=%d j=%d k=%d !\n", i,j,k);	
				if( matrix == TT )
					printf("ERROR in traceback on table S -> i=%d j=%d k=%d !\n", i,j,k);	
				if( matrix == TS )
					printf("ERROR in traceback on table T -> i=%d j=%d k=%d !\n", i,j,k);	

				exit(-1);								
			break;
		}		

	}

	res1[len1] = '\0';
	res2[len2] = '\0';		
	reverse(res1);
	reverse(res2);	
	printf("%s\n%s\n\n", res1, res2);	
		
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

void remove_traceback_tables(){
	int i, j;

	for( i = 0; i < M+1; ++i ){
		for( j = 0; j < N+1; ++j ){
			free( TR[i][j] );
			free( TS[i][j] );
			free( TT[i][j] );
		}

		free( TR[i] );
		free( TS[i] );
		free( TT[i] );
	}

	free( TR );
	free( TS );
	free( TT );

}