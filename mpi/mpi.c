#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <assert.h>

#define DEFAULT 0
#define ORDEM 2048      /* Ordem da matriz */
#define QTD_GEN 2000    /* Quantidade de gerações */
#define VIVA 1          /* Define a vida da celula de maneira booleana */
#define MORTA 0         /* Define a morte da celula de maneira booleana */


int** matriz, ** matrizAtu; /* Definição das matrizes */
int bufferEnvio[ORDEM], bufferRecebimento[ORDEM]; /* Buffers para troca de informação */

/* Verifica quantos vizinhos vivos uma célula tem */
int verificaVizinhosVivos(int i, int j) {
    int quantidadeVivos = 0;

    quantidadeVivos += matriz[i][((j + 1) % ORDEM)];
    quantidadeVivos += matriz[((i + 1) % ORDEM)][((j + 1) % ORDEM)];
    quantidadeVivos += matriz[((i + 1) % ORDEM)][j];
    quantidadeVivos += matriz[((i + 1) % ORDEM)][(ORDEM + (j - 1)) % ORDEM];
    quantidadeVivos += matriz[i][(ORDEM + (j - 1)) % ORDEM];
    quantidadeVivos += matriz[(ORDEM + (i - 1)) % ORDEM][(ORDEM + (j - 1)) % ORDEM];
    quantidadeVivos += matriz[(ORDEM + (i - 1)) % ORDEM][j];
    quantidadeVivos += matriz[(ORDEM + (i - 1)) % ORDEM][((j + 1) % ORDEM)];

    return quantidadeVivos;
}

/* Configura e inicia a próxima geração */
void proximaGeracao(int processo, int rank) {
    int i, j;
    int limite = processo * rank;

    for (i = limite; i < processo * (rank + 1); i++) {
        for (j = 0; j < ORDEM; j++) {
            /* Aplicação das regras */
            if (matriz[i][j] == 1) { /* Verifica as celulas vivas */
                if (verificaVizinhosVivos(i, j) < 2 || verificaVizinhosVivos(i, j) > 3)
                    matrizAtu[i][j] = MORTA;
                else
                    matrizAtu[i][j] = VIVA;
            }
            else { /* Verifica as celulas mortas */
                if (verificaVizinhosVivos(i, j) == 3)
                    matrizAtu[i][j] = VIVA;
                else
                    matrizAtu[i][j] = MORTA;
            }
        }
    }

    for (i = 0; i < ORDEM; i++) {
        for (j = 0; j < ORDEM; j++) {
            matriz[i][j] = matrizAtu[i][j];
        }
    }
}

/* Verifica o total de células vivas no momento atual */
int verificaTotalVivas() {
    int i, j, qtdVivas = 0;
    for (i = 0; i < ORDEM; i++) {
        for (j = 0; j < ORDEM; j++) {
            if (matriz[i][j])
                qtdVivas++;
        }
    }
    return qtdVivas;
}

/* Execução do processo primario */
void processoPrimario(int qtdProc) {
    int i, j, k, tag = 0, div = 0, rank, atual;
    int processoAtu = ORDEM / qtdProc, linha, coluna;

    clock_t begin = clock();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Separa a memoria para as matrizes */
    matriz = malloc(sizeof(int*) * ORDEM);
    matrizAtu = malloc(sizeof(int*) * ORDEM);

    for (i = 0; i < ORDEM; i++) {
        matriz[i] = malloc(sizeof(int) * ORDEM);
        matrizAtu[i] = malloc(sizeof(int) * ORDEM);
    }

    /* Marca os pontos iniciais */
    for (i = 0; i < (processoAtu); i++) {
        for (j = 0; j < ORDEM; j++) {
            matriz[i][j] = 0;
        }
    }

    /* GLIDER */
    linha = 1; coluna = 1;
    matriz[linha][coluna + 1] = 1;
    matriz[linha + 1][coluna + 2] = 1;
    matriz[linha + 2][coluna] = 1;
    matriz[linha + 2][coluna + 1] = 1;
    matriz[linha + 2][coluna + 2] = 1;

    /* R-pentomino */
    linha = 10; coluna = 30;
    matriz[linha][coluna + 1] = 1;
    matriz[linha][coluna + 2] = 1;
    matriz[linha + 1][coluna] = 1;
    matriz[linha + 1][coluna + 1] = 1;
    matriz[linha + 2][coluna + 1] = 1;

    div = i;

    for (k = 1; k < qtdProc; k++) {
        for (i = div; i < (k + 1) * processoAtu; i++) {
            MPI_Recv(bufferRecebimento, ORDEM, MPI_INT, k, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (j = 0; j < ORDEM; j++) {
                matriz[i][j] = bufferRecebimento[j];
            }
        }

        div = i;
    }

    printf("Condicao Inicial: %d\n", verificaTotalVivas());

    /* Crias as proximas geracoes a partir da primeira  */
    for (atual = 0; atual < QTD_GEN; atual++) {
        for (i = 0; i < ORDEM; i++) { // Envia tabela
            for (j = 0; j < ORDEM; j++) bufferEnvio[j] = matriz[i][j];
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(bufferEnvio, ORDEM, MPI_INT, DEFAULT, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        proximaGeracao(processoAtu, rank);

        for (k = 1; k < qtdProc; k++) {
            for (i = (processoAtu * k); i < (processoAtu * (k + 1)); i++) {
                MPI_Recv(bufferRecebimento, ORDEM, MPI_INT, k, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (j = 0; j < ORDEM; j++) {
                    matriz[i][j] = bufferRecebimento[j];
                }
            }
        }
    }

    printf("Ultima Geracao: %d\n", verificaTotalVivas());
    clock_t end = clock();

    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("%lf\n", time_spent);

}

/* Execução do processo secundario */
void processoSecundario(int qtdProc) {
    int tag = 0, dest = 0, i, j, atual, rank, processoAtu = ORDEM / qtdProc;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Separa a memoria para as matrizes */
    matriz = malloc(sizeof(int*) * ORDEM);
    matrizAtu = malloc(sizeof(int*) * ORDEM);
    for (i = 0; i < ORDEM; i++) {
        matriz[i] = malloc(sizeof(int) * ORDEM);
        matrizAtu[i] = malloc(sizeof(int) * ORDEM);
    }

    /* Marca a segunda parte na primeira geração */
    for (i = 0; i < ORDEM; i++) {
        for (j = 0; j < ORDEM; j++) {
            if (i >= rank * processoAtu && i < (rank + 1) * processoAtu) {
                matriz[i][j] = 0;
            }
        }
    }

    for (i = rank * processoAtu; i < (rank + 1) * processoAtu; i++) {
        for (j = 0; j < ORDEM; j++) bufferEnvio[j] = matriz[i][j];
        MPI_Send(bufferEnvio, ORDEM, MPI_INT, dest, tag, MPI_COMM_WORLD);
    }

    for (atual = 0; atual < QTD_GEN; atual++) {
        for (i = 0; i < ORDEM; i++) {
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(bufferRecebimento, ORDEM, MPI_INT, DEFAULT, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            for (j = 0; j < ORDEM; j++) matriz[i][j] = bufferRecebimento[j];
        }

        proximaGeracao(processoAtu, rank);

        for (i = (processoAtu * rank); i < processoAtu * (rank + 1); i++) {
            for (j = 0; j < ORDEM; j++) bufferEnvio[j] = matriz[i][j];
            MPI_Send(bufferEnvio, ORDEM, MPI_INT, dest, tag, MPI_COMM_WORLD);
        }
    }
}

int main(int argc, char* argv[]) {
    int rank;           /* ID do processo */
    int qtdProcessos;   /* Quantidade de processos */


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &qtdProcessos);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        processoPrimario(qtdProcessos);
    }
    else {
        processoSecundario(qtdProcessos);
    }

    MPI_Finalize();

    return 0;
}
