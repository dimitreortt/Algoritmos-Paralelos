// Programa: multMatrizes.c
// Autor: Dimitre Ortt
// Data: 02/12/2018
//
// O Diálogo: Este programa recebe como entrada
// um inteiro representando a dimensão 'n'
// e duas matrizes da entrada padrão e realiza 
// o calculo da multiplicação das matrizes 
// utilizando a técnica vista em sala de 
// particionamento das matrizes, uma em faixas
// horizontaislinhas e a outra em faixas 
// verticais. 
//
// É utilizada a biblioteca MPI e suas 
// funcionalidades: MPI_Bcast, MPI_Scatter 
// (envio simultâneo de dados para todos os 
// outros processadores) e MPI_Gather para 
// coleta de dados de todos os outros 
// processadores.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

/*Obs1.(IMPORTANTE): ressaltamos que este aplicativo NÃO FUNCIONA corretamente caso as dimensões
	'n' não sejam múltiplas do número de processadores dados como parâmetros na execução*/

/*Obs2.: foram implementados 4 métodos auxiliares para formatação dos dados
	para tráfego e para os cálculos, todo tráfego é feito através de ponteiros
		de inteiros, e todos os cálculos são feitos com matrizes devidamente formatadas*/

/*Ressaltamos que tais métodos foram de grande valia e facilitou muito a resolução do problema*/

/*recebe uma matriz src e um array dst (com meméria alocada) de destino, e transfere os dados
	de src para dst **linha por linha** */
void MatrixToVecHorizontal(int num_rows, int num_cols, int ** src, int * dst){
	int i;

	for(i = 0; i < num_rows*num_cols; i++){
		dst[i] = src[i/num_cols][i%num_cols];
	}
}

/*recebe uma matriz src e um array (com meméria alocada) de destino, e transfere os dados
	de src para dst **coluna por coluna** */
void MatrixToVecVertical(int num_rows, int num_cols, int ** src, int * dst){
	int i;

	for(i = 0; i < num_rows*num_cols; i++){
		dst[i] = src[i%num_rows][i/num_rows];
	}
}

/*realiza o inverso de MatrixToVecHorizontal*/
void VecToMatrixHorizontal(int * src, int size, int num_rows, int num_cols, int ** dst){
	int i;
	for(i = 0; i < size; i++){
		dst[i/num_cols][i%num_cols] = src[i];
	}
}

/*realiza o inverso de MatrixToVecVertical*/
void VecToMatrixVertical(int * src, int size, int num_rows, int num_cols, int ** dst){
	int i;
	for(i = 0; i < size; i++){
		dst[i%num_rows][i/num_rows] = src[i];
	}
}

/*prodOfLanes calcula o produto A'B' e armazena em C' com o devido offset assim
	como no algoritmo do slide*/
/*é interessante observar que é fundamental que o offset recebido esteja correto*/
void prodOfLanes(int lane_sz, int n, int ** a, int ** b, int offset, int ** c){
	int i, j, k, sum;
	int off_set = offset*lane_sz;

	for(i = 0; i < lane_sz; i++){
		for(j = 0; j < lane_sz; j++){
			sum = 0;
			for(k = 0; k < n; k++){
				sum += a[i][k]*b[k][j];
			}
			//esta aritmética insere os dados corretamente em C'
			c[i][j + off_set] = sum;
		}
	}	
}

/*função para alocação de memória para uma matriz*/
void aloc(int *** dst, int nrows, int ncols){
	int i;
	*dst = (int **)malloc(nrows*sizeof(int*));
	for(i = 0; i < nrows; i++){
		(*dst)[i] = (int*)malloc(ncols*sizeof(int));
	}
}

/*imprime uma matriz*/
void printMatrix(int nrows, int ncols, int ** M){
	int i, j;

	printf("\n");
	for(i = 0; i < nrows; i++){
		for(j = 0; j < ncols; j++){
			printf("%3d ", M[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

int main(void){
	int comm_sz;
	int my_rank;
	int n; // n corresponde à dimensão das matrizes
	int raiz_p;
	int i;
	int j;
	int * abuf; // abuf armazena A em formato de array
	int * bbuf; // bbuf armazena B em formato de array
	int ** A;  
	int ** B;
	MPI_Status status;

	/*devida inicialização*/
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if(my_rank == 0){

		printf("Entre com a dimensão 'n' das matrizes A e B\n\n");
		scanf("%d", &n);
		printf("As dimensoes sao: %d\n\n", n);

		aloc(&A, n, n);
		aloc(&B, n, n);
		raiz_p = (int)sqrt(comm_sz);

		printf("Entre com a matriz A:\n");
		for(i = 0; i < n; i++){
			for(j = 0; j < n; j++){
				scanf("%d", &A[i][j]);
			}
		}
		printMatrix(n, n, A);

		printf("entre com a matriz B:\n");
		for(i = 0; i < n; i++){
			for(j = 0; j < n; j++){
				scanf("%d", &B[i][j]);
			}
		}
		printMatrix(n, n, B);

		abuf = (int*)malloc(n*n*sizeof(int));
		bbuf = (int*)malloc(n*n*sizeof(int));

		/*traduz A e B para formato de array, para movimentação de dados*/
		MatrixToVecHorizontal(n, n, A, abuf);
		MatrixToVecVertical(n, n, B, bbuf);
	}
	
	/*sincronização*/
	MPI_Barrier(MPI_COMM_WORLD);

	/*comunicação a respeito das dimensoes n*/
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/*lane_sz, tamanho da 'faixa', corresponde a dimensão n dividida entre os processadores*/
	int lane_sz = n/comm_sz;

	int ** a; 				 // a minúsculo corresponde à A' no algoritmo do slide
	aloc(&a, lane_sz, n);	 //(int**)malloc(lane_sz*sizeof(int*));
	int ** b; 				 // B'
	aloc(&b, n, lane_sz);	 //= (int**)malloc(n*sizeof(int*));
	int ** c;	 			 // C'
	aloc(&c, lane_sz, n);	 //= (int**)malloc(lane_sz*sizeof(int*));

	/*A' e B' em formato de vetor, para o tráfego de dados*/
	int * vec_a = (int *)calloc(lane_sz*n, sizeof(int));
	int * vec_b = (int *)calloc(lane_sz*n, sizeof(int));

	/*sincronização*/
	MPI_Barrier(MPI_COMM_WORLD);
	
	/*método para distribuição de dados das matrizes A e B pelo processo 0*/
	MPI_Scatter(abuf, lane_sz*n, MPI_INT, vec_a, lane_sz*n, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(bbuf, lane_sz*n, MPI_INT, vec_b, lane_sz*n, MPI_INT, 0, MPI_COMM_WORLD);

	/*sincronização*/
	MPI_Barrier(MPI_COMM_WORLD);	
	
	/*tradução dos dados recebidos em formato de vetor para formato de matrizes*/
		/*'Horizontal' e 'Vertical' dizem à respeito do protocolo de envio de dados */
	VecToMatrixHorizontal(vec_a, lane_sz*n, lane_sz, n, a);
	VecToMatrixVertical(vec_b, lane_sz*n, n, lane_sz, b);

//PASSO 1.
	/*a função prodOfLanes calcula o produto de duas matrizes compatíveis a e b e armazena
		o resultado na matriz c (com memória alocada anteriormente)*/
	prodOfLanes(lane_sz, n, a, b, my_rank, c);
	//void prodOfLanes(int lane_sz, int n, int ** a, int ** b, int offset, int ** c){

//PASSO 2.
	for(j = 1; j < comm_sz; j++){
//PASSO 2.1.
		/*intercâmbio de dados*/
		MPI_Send(vec_b, lane_sz*n, MPI_INT, (my_rank+1)%comm_sz, 0, MPI_COMM_WORLD);
		MPI_Recv(vec_b, lane_sz*n, MPI_INT, (my_rank-1+comm_sz)%comm_sz, 0, MPI_COMM_WORLD, &status);

//PASSO 2.2.
		/*calcula C' e armazena em C' com o devido offset*/
		VecToMatrixVertical(vec_b, lane_sz*n, n, lane_sz, b);
		/*note a aritmética para cálculo do offset correto*/
		prodOfLanes(lane_sz, n, a, b, (my_rank-j+comm_sz)%comm_sz, c);
	}

	/*reunião dos dados*/
	MatrixToVecHorizontal(lane_sz, n, c, vec_b);
	MPI_Gather(vec_b, lane_sz*n, MPI_INT, bbuf, lane_sz*n, MPI_INT, 0, MPI_COMM_WORLD);

	/*processador 0 é o encarregado de imprimir o resultado*/
	if(my_rank == 0){
		VecToMatrixHorizontal(bbuf, n*n, n, n, B);
		printf("\n<<<<<<<<  A matriz resultado C é: >>>>>>>>>\n");
		printMatrix(n, n, B);
	}

	/*conclusão apropriada da aplicação*/
	MPI_Finalize();
	return 0;
}