// Programa: multMatrizes.c
// Programadores: Dimitre Ortt e Matheus Delmondes
// Data: 02/12/2018

// O Diálogo: Este programa recebe como entrada
// um inteiro representando a dimensão 'n'
// e duas matrizes da entrada padrão e realiza 
// o calculo da multiplicação das matrizes 
// utilizando a técnica vista em sala de 
// particionamento das matrizes em blocos e 
// distribuindo-os entre os processadores.

// Utiliza-se a biblioteca MPI e suas 
// funcinalidades:
// MPI_Bcast, MPI_Scatter (envio simultâneo de 
// dados para todos os outros processadores) e
// MPI_Gather para coleta de dados de todos os
// outros processadores.

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

/*Ressaltamos que tais métodos foram de grande valia e facilitaram muito a resolução do problema*/

/*recebe a matriz B lida e a traduz para um array da maneira ideal:
	após separar B em blocos de dimensão n/raiz_p, escreve no array destino os blocos
	um a um, começando do mais acima à esquerda, escrevendo-o linha a linha,
	e então passando para o próximo bloco e se deslocando nos blocos em um uma 'faixa' <<Horizontal>>*/
void MatrixBToVecBlocks(int n, int raiz_p, int ** src, int * dst){
	int pos = 0, i, j, k, l;

	for(i = 0; i < raiz_p; i++){
		for(j = 0; j < raiz_p; j++){
			for(k = 0; k < n/raiz_p; k++){
				for(l = 0; l < n/raiz_p; l++){
					dst[pos] = src[(n/raiz_p)*j + k][(n/raiz_p*i) + l];
					pos++;
				}
			}
		}
	}
}

/*recebe a matriz B lida e a traduz para um array da maneira ideal:
	após separar B em blocos de dimensão n/raiz_p, escreve no array destino os blocos
	um a um, começando do mais acima à esquerda, escrevendo-o linha a linha,
	e então passando para o próximo bloco se deslocando nos blocos em um uma 'faixa' <<Vertical>>*/
void MatrixAToVecBlocks(int n, int raiz_p, int ** src, int * dst){
	int pos = 0, i, j, k, l;

	for(i = 0; i < raiz_p; i++){
		for(j = 0; j < raiz_p; j++){
			for(k = 0; k < n/raiz_p; k++){
				for(l = 0; l < n/raiz_p; l++){
					dst[pos] = src[(n/raiz_p)*i + k][(n/raiz_p*j) + l];
					pos++;
				}
			}
		}
	}
}

/*recebe um bloco em seu formato de array e o trancreve para uma matriz correspondente, que
	apelidamos de bloco*/
void VecToBlock(int * src, int size, int ** dst){
	int pos = 0, i, j, k, l;

	for(i = 0; i < size; i++){
		for(j = 0; j < size; j++){
			dst[i][j] = src[pos];
			pos++;
			//bbuf[(i*n/raiz_p) * (j*n/raiz_p) + k*n/raiz_p] = A[n/raiz_p*i + k][raiz_p*j + l];
		}
	}
}

/*recebe um bloco e o transcreve ao seu formato de array*/
void BlockToVec(int ** src, int size, int * dst){
	int pos = 0, i, j, k, l;

	for(i = 0; i < size; i++){
		for(j = 0; j < size; j++){
			dst[pos] = src[i][j];
			pos++;
			//bbuf[(i*n/raiz_p) * (j*n/raiz_p) + k*n/raiz_p] = A[n/raiz_p*i + k][raiz_p*j + l];
		}
	}
}


/*prodOfMatrices calcula o produto das matrizes A'B' e armazena na matriz C', todas com
	memória alocada antecipadamente*/
	/*note que a minúsculo representa A' (resp. b, c)*/
void prodOfMatrices(int size, int ** a, int ** b, int ** c){
	int i, j, k, sum;

	for(i = 0; i < size; i++){
		for(j = 0; j < size; j++){
			sum = 0;
			for(k = 0; k < size; k++){
				sum += a[i][k]*b[k][j];
			}
			//esta aritmética insere os dados corretamente em C'
			c[i][j] = sum;
		}
	}	
}

/*calcula a soma das matrizes a e b e armazena em c*/
void sumOfMatrices(int size, int ** a, int ** b, int ** c){
	int i, j, k, sum;

	for(i = 0; i < size; i++){
		for(j = 0; j < size; j++){
			c[i][j] = a[i][j] + b[i][j];
		}
	}	
}

/*escreve parcela dos resultados calculados na porção de resultados da Matriz final C,
	cada com id múltiplo de raiz_p é responsável por essa escrita, esta função é invocada
	raiz_p vezes, e em cada vez com um offset adequado*/
void WriteC(int size, int ** src, int offset, int ** c){
	int i, j;
	int off_set = size*offset;

	for(i = 0; i < size; i++){
		for(j = 0; j < size; j++){
			c[i][j + off_set] = src[i][j];
		}
	}
}

/*imprime a saída com os resultados no meio da saída definido*/
void printC(int block_size, int raiz_p, int ** c){
	int i, j;
	for(i = 0; i < block_size; i++){
		for(j = 0; j < raiz_p*block_size; j++){
			printf("%4d ", c[i][j]);
		}
		printf("\n");
	}
	//printf("\n");
}

/*função para alocação de memória para uma matriz*/
void aloc(int *** dst, int nrows, int ncols){
	int i;
	*dst = (int **)calloc(nrows, sizeof(int*));
	for(i = 0; i < nrows; i++){
		(*dst)[i] = (int*)calloc(ncols, sizeof(int));
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
		MatrixAToVecBlocks(n, raiz_p, A, abuf);
		MatrixBToVecBlocks(n, raiz_p, B, bbuf);
		//MatrixToVecVertical(n, n, B, bbuf);

	}

	/*comunicação a respeito das dimensoes n*/
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&raiz_p, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/*sincronização*/
	MPI_Barrier(MPI_COMM_WORLD);

	/*block_size corresponde as dimensões dos blocos distribuídos 
		entre os processadores*/
	int block_size = n/raiz_p;
	
	/*a seguir são declaradas matrizes para realização dos cálculos
		A', B' correspondem à porção dos dados recebidos de A e B*/
	int ** a;
	int ** b;
	/*C' é a matriz que comportará os resultados*/
	int ** c;
	/*D' é uma matriz para o cálculo intermediário de A'B'*/
	int ** d;
	/*E' é uma matriz auxiliar para cálculos que se fizeram necessários*/
	int ** e;

	aloc(&a, block_size, block_size);	 
	aloc(&b, block_size, block_size);
	aloc(&c, block_size, raiz_p*block_size);	
	aloc(&e, block_size, block_size);
	aloc(&d, block_size, block_size);

	/*sincronização*/
	MPI_Barrier(MPI_COMM_WORLD);

	
	/*A', B', C' e D' em formato de vetor, para o tráfego de dados*/
	int * vec_a = (int *)calloc((n*n)/comm_sz, sizeof(int));
	int * vec_b = (int *)calloc((n*n)/comm_sz, sizeof(int));
	int * vec_c = (int *)calloc((n*n*raiz_p)/comm_sz, sizeof(int));
	int * vec_d = (int *)calloc((n*n)/comm_sz, sizeof(int));

	/*sincronização*/
	MPI_Barrier(MPI_COMM_WORLD);	

	/*variável que armazena o número de dados em cada bloco*/
	int lengh = (n*n)/comm_sz;

	/*distribuição dos dados das matrizes A e B já formatados, entre os processadores*/
	MPI_Scatter(abuf, lengh, MPI_INT, vec_a, (n*n)/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(bbuf, lengh, MPI_INT, vec_b, (n*n)/comm_sz, MPI_INT, 0, MPI_COMM_WORLD);

	/*sincronização*/
	MPI_Barrier(MPI_COMM_WORLD);
	
	/*resgate dos dados recebidos, passando-os para sua forma em blocos*/
	VecToBlock(vec_a, block_size, a);
	VecToBlock(vec_b, block_size, b);

	/*Marker: Neste ponto a distribuição dos dados foi devidamente realizada
		agora começamos os passos do algoritmo*/

//PASSO 1.
	/*calcula o produto A'B' de dimensoes block_size e armazena em D*/
	prodOfMatrices(block_size, a, b, d);

	/*tradução para formato de array para tráfego dos dados*/
	BlockToVec(d, block_size, vec_d);

	/*definição dos i e j locais do Pij na representação dos processadores cobrindo as matrizes*/
	int my_i = my_rank/raiz_p;
	int my_j = my_rank%raiz_p;

	if(my_j != 0){
	//PASSO 2.
		MPI_Send(vec_d, lengh, MPI_INT, my_i*raiz_p, 0, MPI_COMM_WORLD);
	}
	else{
	//PASSO 3.
		for(j = 1; j < raiz_p; j++){
			MPI_Recv(vec_d, lengh, MPI_INT, my_i*raiz_p + j, 0, MPI_COMM_WORLD, &status);
			VecToBlock(vec_d, block_size, e);

			//void sum (int src_dimension, int ** src1, int ** src2, int ** dst);
			sumOfMatrices(block_size, e, d, d);

		}
		/*if(my_rank == 20)
			printMatrix(block_size, block_size, d);*/
		//void WriteC(int src_dimension, int root_p, int ** src, int offset, int ** dst);
		WriteC(block_size, d, my_i%raiz_p,c);
	}

	//PASSO 4.
	for(i = 1; i < raiz_p; i++){

		/*adequa bloco B' para ser enviada*/
		BlockToVec(b, block_size, vec_b);

		//PASSO 4.1.
		MPI_Send(vec_b, lengh, MPI_INT, (my_rank-raiz_p+comm_sz)%comm_sz, 0, MPI_COMM_WORLD);
		MPI_Recv(vec_b, lengh, MPI_INT, (my_rank+raiz_p)%comm_sz, 0, MPI_COMM_WORLD, &status);
		
		/*resgata valor recebido para o bloco B'*/
		VecToBlock(vec_b, block_size, b);
		
		//PASSO 4.2.
		prodOfMatrices(block_size, a, b, d);
		
		BlockToVec(d, block_size, vec_d);

		//PASSO 4.3.
		if(my_j != 0){
			MPI_Send(vec_d, lengh, MPI_INT, my_i*raiz_p, 0, MPI_COMM_WORLD);
		}
		else{
			//PASSO 4.4.1
			for(j = 1; j < raiz_p; j++){
				memset(vec_d, 0, sizeof(vec_d));
				MPI_Recv(vec_d, lengh, MPI_INT, my_i*raiz_p + j, 0, MPI_COMM_WORLD, &status);
				VecToBlock(vec_d, block_size, e);

				sumOfMatrices(block_size, e, d, d);
			}
			//PASSO 4.4.2
			WriteC(block_size, d, (my_i+i)%raiz_p,c);					
		}
	}

	/*sincronização*/
	MPI_Barrier(MPI_COMM_WORLD);

	/*a partir daqui é feito um protocolo para a devida escrita da matriz resultado C
		o processador 0 escreve sua parcela da resposta em C e em seguida ordena as impressões
		da cada processador também encarregado de tal função*/
	if(my_rank == 0){
		/*imprime resultado armazenado em C na saída definida (globalmente)*/
		printC(block_size, raiz_p, c);
		/*cada processador espera sua vez de imprimir, e envia a resposta de que terminou tal processo*/
		for(i = 1; i < raiz_p; i++){
			/*estes Send e Recv possuem a exclusiva função de sincronização e não trafegam
				dados relevantes*/
			MPI_Send(vec_b, 1, MPI_INT, i*raiz_p, 0, MPI_COMM_WORLD);
			MPI_Recv(vec_b, 1, MPI_INT, i*raiz_p, 0, MPI_COMM_WORLD, &status);
		}
	}
	else if(my_j == 0){
		/*processadores cujo my_j é 0 são encarregados de direcionar a saída*/

		/*espera a sua vez*/
		MPI_Recv(vec_b, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		/*imprime*/
		printC(block_size, raiz_p, c);
		/*alega que terminou para o próximo receber autorização*/
		MPI_Send(vec_b, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	/*conclusão apropriada da aplicação*/
	MPI_Finalize();
	return 0;
}

