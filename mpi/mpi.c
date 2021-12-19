#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include<sys/time.h>

#define NUM_GEN 2000 // Numero de geracoes
#define TAM 2048 // Tamanho N da matriz NxN
#define c_viva 1
#define c_morta 0
#define PRINCIPAL 0

int **grid, **newgrid; // matrizes de geracao
int bufSnd[TAM], bufRcv[TAM]; // Buffers de comunicacao

/// Calculo de tempo
typedef struct {
    int secs;
    int usecs;
}Duracao;

Duracao * tempo_decorrido (struct timeval *start, struct timeval *end){
    Duracao *total = (Duracao*)malloc(sizeof(Duracao));

    if(start->tv_sec == end->tv_sec){
        total->secs = 0;
        total->usecs = end->tv_usec - start->tv_usec;
    }
    else{
        total->usecs = 1000000 - start->tv_usec;
        total->secs = end->tv_sec - (start->tv_sec + 1);
        total->usecs += end->tv_usec;
        if(total->usecs >= 1000000){
            total->usecs -= 1000000;
            total->secs += 1;
        }
    }
    return total;
}

// Retorna a quantidade de vizinhos vivos de cada celula
int getNeighbors(int i, int j) {
    int num_vizinhos=0;

    num_vizinhos += grid[i][((j+1)%TAM)];
    num_vizinhos += grid[((i+1)%TAM)][((j+1)%TAM)];
    num_vizinhos += grid[((i+1)%TAM)][j];
    num_vizinhos += grid[((i+1)%TAM)][(TAM+(j-1))%TAM];
    num_vizinhos += grid[i][(TAM+(j-1))%TAM];
    num_vizinhos += grid[(TAM+(i-1))%TAM][(TAM+(j-1))%TAM];
    num_vizinhos += grid[(TAM+(i-1))%TAM][j];
    num_vizinhos += grid[(TAM+(i-1))%TAM][((j+1)%TAM)];

    return num_vizinhos;
}

// Cria uma nova geracao de acordo com as regras estabelecidas
void prox_rodada(int proc_atual, int rank){
    int i, j;
    int limite = proc_atual*rank;

    for(i=limite;i<proc_atual*(rank+1); i++){
        for(j = 0; j<TAM; j++){
            if (grid[i][j] == 1){ //le celulas vivas
                if (getNeighbors(i,j) < 2 || getNeighbors(i,j) > 3) // Regra 1
                    newgrid[i][j] = c_morta;
                else // Regra 3
                    newgrid[i][j] = c_viva;
            }
            else{ // le células mortas
                if(getNeighbors(i,j) == 3) // Regra 2
                    newgrid[i][j] = c_viva;
                else // Regra 3
                    newgrid[i][j] = c_morta;
            }
        }
    }

    for(i=0;i<TAM; i++){
        for(j = 0; j<TAM; j++){
            grid[i][j] = newgrid[i][j];
        }
    }
}

// Conta quantas celulas estao vivas na geracao
int total_vivos(){
    int i,j,numv = 0;
        for(i=0;i<TAM; i++){
            for(j = 0; j<TAM; j++){
                if (grid[i][j])
                    numv++;
            }
        }
    return numv;
}

// Rotina do processo principal
void prinProc(int nproc){
    int i,j,inicio,tag=0,div=0,rank,atual;
    int proc_atual = TAM/nproc, lin, col;
    Duracao *valor;
    struct timeval start, end;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    gettimeofday (&start, NULL);

    // Alocacao das matrizes
    grid = malloc(sizeof(int*)*TAM);
    newgrid = malloc(sizeof(int*)*TAM);
    for(i=0;i<TAM;i++){
        grid[i] = malloc(sizeof(int)*TAM);
        newgrid[i] = malloc(sizeof(int)*TAM);
    }
    /// Definição dos pontos iniciais
    for(i=0;i<(proc_atual); i++){
        for(j = 0; j<TAM; j++){
            grid[i][j] = 0;
        }
    }
  //GLIDER
    lin = 1; col = 1;
    grid[lin  ][col+1] = 1;
    grid[lin+1][col+2] = 1;
    grid[lin+2][col  ] = 1;
    grid[lin+2][col+1] = 1;
    grid[lin+2][col+2] = 1;

    //R-pentomino
    lin = 10; col = 30;
    grid[lin  ][col+1] = 1;
    grid[lin  ][col+2] = 1;
    grid[lin+1][col  ] = 1;
    grid[lin+1][col+1] = 1;
    grid[lin+2][col+1] = 1;

    div = i;
    for(inicio=1;inicio<nproc;inicio++){
        for(i=div;i<(inicio+1)*proc_atual;i++){
            MPI_Recv(bufRcv,TAM,MPI_INT,inicio,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            for(j=0;j<TAM;j++){
                grid[i][j] = bufRcv[j];
            }
        }
        div=i;
    }

    printf("Condicao Inicial: %d\n", total_vivos());

    // Gera NUM_GEN geracoes a partir da primeira
    for(atual=0;atual<NUM_GEN;atual++){
        for(i=0;i<TAM;i++){ // Envia tabela
            for(j=0;j<TAM;j++) bufSnd[j] = grid[i][j];
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(bufSnd,TAM,MPI_INT,PRINCIPAL,MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        prox_rodada(proc_atual,rank);

        for(inicio=1;inicio<nproc;inicio++){
            for(i=(proc_atual*inicio);i<(proc_atual*(inicio+1));i++){
                MPI_Recv(bufRcv,TAM,MPI_INT,inicio,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                for(j=0;j<TAM;j++){
                    grid[i][j] = bufRcv[j];
                }
            }
        }
    }

    printf("Ultima Geracao: %d\n", total_vivos());

    gettimeofday (&end, NULL);
    valor = tempo_decorrido(&start, &end);
    printf("Tempo: %d,%d s\n",valor->secs,valor->usecs);
}

// Rotina do processo secundario
void secProc(int nproc){
    int tag=0,dest=0,i,j,atual,rank,proc_atual=TAM/nproc;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // Alocacao das matrizes
    grid = malloc(sizeof(int*)*TAM);
    newgrid = malloc(sizeof(int*)*TAM);
    for(i=0;i<TAM;i++){
        grid[i] = malloc(sizeof(int)*TAM);
        newgrid[i] = malloc(sizeof(int)*TAM);
    }

    // Gera a segunda metade da primeira geracao
    for(i=0;i<TAM; i++){
        for(j = 0; j<TAM; j++){
            if(i >= rank*proc_atual && i < (rank+1)*proc_atual){
                grid[i][j] = 0;
            }
        }
    }

    for(i=rank*proc_atual;i<(rank+1)*proc_atual;i++){
        for(j=0;j<TAM;j++) bufSnd[j] = grid[i][j];
        MPI_Send(bufSnd,TAM,MPI_INT,dest,tag,MPI_COMM_WORLD);
    }

    for(atual=0;atual<NUM_GEN;atual++){
        for(i=0;i<TAM;i++){ // recebe a tabela
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(bufRcv,TAM,MPI_INT,PRINCIPAL,MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            for(j=0;j<TAM;j++) grid[i][j] = bufRcv[j];
        }

        prox_rodada(proc_atual,rank);

        for(i=(proc_atual*rank);i<proc_atual*(rank+1);i++){ // envia tabela
            for(j=0;j<TAM;j++) bufSnd[j] = grid[i][j];
            MPI_Send(bufSnd,TAM,MPI_INT,dest,tag,MPI_COMM_WORLD);
        }
    }
}

int main(int argc, char *argv[]){
    int rank; // ID do processo
    int nproc; // numero de processos

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0){ // caso seja o processo principal
        prinProc(nproc);
    }else{
        secProc(nproc);
    }

    MPI_Finalize();
    return 0;
}
