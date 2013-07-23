#include <stdio.h>
#include <stdlib.h>
#include <string.h>	

typedef int ** tabel;

void init_subs_lookup_tabel( int *** S, int *** A, char * filename );
void init_dynamic_tabel(int *** t, int ** A, int **S, char *seq1, char *seq2 );
void init_traceback_tabel( char *** T, int n, int m );

void align( int ** M, int **S, int ** A, char **T, char *s1, char *s2, int n, int m);
void traceback( int **M, char **T, char *s1, char *s2, char *s1_, char *s2_ );

void print_tabel(tabel M, int n, int m);
void print_tabel_(char ** M, int n, int m);
void reverse(char s[]);

int compare(const int ** a,const int ** b);
int get_translation(int **A, int c);

int len_alphabet;

int main(int argc, char const *argv[]) {

	if( argc <= 3 ){
		printf("Usage: %s <seq1> <seq2> <table filename>\n", argv[0]);
		return 0;
	}

	/* Tabela de programacao dinamica */
	int ** M;
	/* Tabela para traceback */
	char ** T;		/* each cell: 'u'-up   'l'-left   'd'-diagonal */
	/* Tabela de subsituicao */
	int ** S;
	/* Tabela lookup de traducao de caracteres do alfabeto em indices na tabela de substituicao */
	int ** A;

	char subs_filename[40];
	strcpy(subs_filename, argv[3]);

	char *seq1 = argv[1];
	char *seq2 = argv[2];
	char seq1_[50];		/* resultado do alinhamento de seq1 */
	char seq2_[50];		/* resultado do alinhamento de seq2 */
	int len_seq1, len_seq2;

	len_seq1 = strlen(seq1);
	len_seq2 = strlen(seq2);

	init_subs_lookup_tabel( &S, &A, subs_filename );
	init_dynamic_tabel( &M, A, S, seq1, seq2 );
	init_traceback_tabel( &T, len_seq1, len_seq2 );

	align(M, S, A, T, seq1, seq2, len_seq1 , len_seq2);
	traceback(M, T, seq1, seq2, seq1_, seq2_);

	reverse(seq1_);
	reverse(seq2_);

	printf("%s\n%s\n", seq1_, seq2_);


	//print_tabel(S, len_alphabet, len_alphabet);
	//print_tabel(A, len_alphabet, 2);
	//print_tabel_( T, len_seq1+1, len_seq2+1);

	return 0;
}


void init_dynamic_tabel(int *** t, int ** A, int **S, char *seq1, char *seq2 ){
	int i, indel_index = len_alphabet - 1;
	int n = strlen(seq1), m = strlen(seq2);

	*t = (int **) malloc( (n+1) * sizeof(int *) );
	for(i=0; i<n+1 ; i++)
		(*t)[i] = (int *) malloc( (m+1) * sizeof(int) );

	// condicoes base 
	for( i = 1; i < n + 1; i++ )
		(*t)[i][0] = S[ i-1 ] [ indel_index ] ;	

	for( i = 1; i < m + 1; i++ )
		(*t)[0][i] = S[ indel_index ] [ i-1 ] ;

	(*t)[0][0] = 0;		
}

void init_traceback_tabel( char *** T, int n, int m ){
	int i;

	*T = (char **) malloc( (n+1) * sizeof(char *) );
	for(i=0; i<n+1 ; i++)
		(*T)[i] = (char *) malloc( (m+1) * sizeof(char) );	

	// condicoes base 
	for( i = 1; i < n + 1; i++ )
		(*T)[i][0] = 'u' ;	

	for( i = 1; i < m + 1; i++ )
		(*T)[0][i] = 'l' ;

	(*T)[0][0] = 'o';		
}

void init_subs_lookup_tabel( int *** S, int *** A, char * filename ){
	/* le tabela de ficheiro com formato */
	FILE *f;
	char line[130];
	char c[2], a;
	int i, j, val;
	len_alphabet = 0;

	f = fopen(filename, "r");

	while( fscanf(f,"#%[^\n]\n", line) ); /* ler comentarios */

	do
		len_alphabet++;
	while( fscanf(f,"%s%c", c, &a)>0 && a!='\n');

	(*A) = (int **) malloc( len_alphabet * sizeof(int *) );		/* alloc lookup table */
	for (i = 0; i < len_alphabet; ++i)
		(*A)[i] = (int *) malloc( 2 * sizeof(int) );

	(*S) = (int **) malloc( len_alphabet * sizeof(int *) );		/* alloc subs table */
	for (i = 0; i < len_alphabet; ++i)
		(*S)[i] = (int *) malloc( len_alphabet * sizeof(int) );
	
	
	for (i = 0; i < len_alphabet; ++i){
		fscanf(f,"%s", c);				
		(*A)[i][0] = c[0];		(*A)[i][1] = i;					/* lookup table */

		for (j = 0; j < len_alphabet; ++j){
			fscanf(f,"%d", &val);
			(*S)[i][j] = val;									/* subs table */
		}
		fscanf(f,"%[^\n]\n", line); 	/* saltar para proxima linha */
	}

	fclose(f);

	/* ordenar lookup para facilitar futuras traducoes */
	qsort( (*A), len_alphabet, sizeof(int[2]), compare );
}

void align( int ** M, int **S, int ** A, char **T, char *s1, char *s2, int n, int m){
	int i, j, m1, m2, m3;
	int indel_index = len_alphabet-1;

	for( i=1; i < n + 1; i++)
		for( j=1; j < m + 1; j++ ){

			m1 = M[i-1][j] + S[ indel_index ][ get_translation(A, s2[j-1]) ] ;
			m2 = M[i][j-1] + S[ get_translation(A, s1[i-1]) ][ indel_index ] ;
			m3 = M[i-1][j-1] + S[ get_translation(A, s1[i-1]) ]
								[ get_translation(A, s2[j-1]) ];

			if( m1 >= m2 && m1 >= m3 ){
				M[i][j] = m1;
				T[i][j] = 'u';
			}
			else if( m2 >= m3 ){
				M[i][j] = m2;
				T[i][j] = 'l';
			}
			else{
				M[i][j] = m3;
				T[i][j] = 'd';
			}		
		}
}

void traceback( int **M, char **T, char *s1, char *s2, char *s1_, char *s2_ ){
	
	int i = strlen(s1) , j = strlen(s2);
	int len1 = 0, len2 = 0;

	while( i!=0 && j!=0 ){
		
		switch( T[i][j] ){
			case 'u':			
				i--;
				s2_[len2++] = '-';
				s1_[len1++] = s1[i];				
			break;
			case 'l':			
				j--;
				s1_[len1++] = '-';
				s2_[len2++] = s2[j];		
			break;
			case 'd':				
				i--; j--;
				s1_[len1++] = s1[i];
				s2_[len2++] = s2[j];							
			break;						
		}
	}
	s1_[len1] = '\0';
	s2_[len2] = '\0';
}

int compare(const int ** a,const int ** b) {
/* ordena por ordem alfabetica na tabela lookup */

  if ( (*a)[0]==(*b)[0] )
    return 0;
  else
    if ( (*a)[0] < (*b)[0] )
        return -1;
     else
      return 1;
}

int get_translation(int **A, int c){
	int i;
	for( i=0; i<len_alphabet && A[i][0]<c ; i++ );
	return A[i][0] == c ? A[i][1] : -1;
}

void print_tabel(int ** M, int n, int m){
	int i, j;
	for(i=0; i<n; i++){
		for(j=0; j<m; j++)
			printf("%d ", M[i][j]);
		
		printf("\n");
	}
}

void print_tabel_(char ** M, int n, int m){
	int i, j;
	for(i=0; i<n; i++){
		for(j=0; j<m; j++)
			printf("%c ", M[i][j]);
		
		printf("\n");
	}
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

