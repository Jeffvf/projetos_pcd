#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

#define tamanho 700
double ideal = 4000;
#define numeropaises 5
#define taxacrescimento 100
#define MAX_THREADS 2

int tempo, cont_grau, decrementador = 0;
int *familias, **relacoes;
int p, desastre, tt, nascidos;
double porcentagens_paises, probabilidade_de_desastre_natural;
double decrementador_desastre, cont_esqucimento = 0, numero_pessoas;
double *paises;
double cont_mudanca, incrementador, porcentagem_real;

// Calculo de tempo
typedef struct
{
    int secs;
    int usecs;
} Duracao;
// função utilizada a partir da leitura do site:
// https://qastack.com.br/programming/5248915/execution-time-of-c-program
Duracao *tempo_decorrido(struct timeval *start, struct timeval *end)
{
    Duracao *total = (Duracao *)malloc(sizeof(Duracao));

    if (start->tv_sec == end->tv_sec)
    {
        total->secs = 0;
        total->usecs = end->tv_usec - start->tv_usec;
    }
    else
    {
        total->usecs = 1000000 - start->tv_usec;
        total->secs = end->tv_sec - (start->tv_sec + 1);
        total->usecs += end->tv_usec;
        if (total->usecs >= 1000000)
        {
            total->usecs -= 1000000;
            total->secs += 1;
        }
    }
    return total;
}

int aleatorio(int n)
{
    return rand() % n;
}

void inicializador()
{
    int i, j;
    relacoes = malloc(sizeof(int *) * tamanho);
    familias = malloc(sizeof(int *) * tamanho);
    paises = malloc(sizeof(int *) * numeropaises);
    for (i = 0; i < tamanho; i++)
    {
        relacoes[i] = malloc(sizeof(int *) * tamanho);
    }

    for (i = 0; i < tamanho; i++)
    {
        for (j = 0; j < tamanho; j++)
        { // esta matriz define as relações entre diferentes familias do pais
            relacoes[i][j] = 0;
        }
        familias[i] = 5; // este e um contador do numero de membros que cada familia possui
        numero_pessoas = numero_pessoas + familias[i];
    }
    for (i = 0; i < numeropaises; i++)
    { // este vetor define as influencias dos paises vizinhos
        paises[i] = 0;
    }
}

// aqui definimos as informações iniciais --> define os valores da matriz e os tipos de relações entre os paises
void relacoes_paises()
{
    int i, j;
    porcentagens_paises = 0;
    // abaixo são definidos os tipos de relações entre paises, fica a escolha de quem esta rodando o codigo
    printf("\n- Decida os tipos de relacoes que cada pais tera com o pais principal\n\n");
    printf("    0 - neutro: nao influencia\n");
    printf("    1 - aliado: influencia positivamente pouco\n");
    printf("    2 - aliado: influencia positivamente muito\n");
    printf("    3 - inimigo: influencia negativamente pouco\n");
    printf("    4 - inimigo: influencia negativamente muito\n");
    printf("\nCaso digite qualquer outro numero, sera considerado como 0\n\n");

    for (i = 0; i < numeropaises; i++)
    { //  influencia aleatoria do pais
        printf("- Pais %d: ", i);
        scanf("%d", &j);
        if (i == 0)
        {
            if (j == 2)
                paises[i] = 15;
            if (j == 1)
                paises[i] = 10;
            if (j == 0)
                paises[i] = 5;
            if (j == 3)
                paises[i] = -15;
            if (j == 4)
                paises[i] = -20;
        }
        else if (i == 1)
        {
            if (j == 2)
                paises[i] = 10;
            if (j == 1)
                paises[i] = 7;
            if (j == 0)
                paises[i] = 2;
            if (j == 3)
                paises[i] = -3;
            if (j == 4)
                paises[i] = -5;
        }
        else if (i == 2)
        {
            if (j == 2)
                paises[i] = 1;
            else
                paises[i] = 0;
        }
        else if (i == 3)
        {
            if (j == 2)
                paises[i] = 5;
            if (j == 1)
                paises[i] = 2.5;
            if (j == 0)
                paises[i] = 1;
            if (j == 3)
                paises[i] = -3;
            if (j == 4)
                paises[i] = -5;
        }
        else
        {
            if (j == 2)
                paises[i] = 10;
            if (j == 1)
                paises[i] = 5;
            if (j == 0)
                paises[i] = 0;
            if (j == 3)
                paises[i] = -5;
            if (j == 4)
                paises[i] = -10;
        }
        // esta variavel representa a porcentagem total de relações que os paises externos influenciam no crescimento da população
        porcentagens_paises = porcentagens_paises + ((paises[i]) / 100);
    }
}

// conta o nome de ligaçõpes que uma celula da matriz pode ter e já monta a matriz
void graus()
{
    int i, j, k, l, grau;
    printf("\n- Digite qual o grau do grafo desejado, sendo ele \n\n    1: grau 2\n    2: grau 4\n    3: grau 8\n\nCaso digite algum numero diferente sera considerado grau 2\n");
    scanf("%d", &grau);
    // aqui pode escolher o grau
    if (grau == 2)
    {
        cont_grau = 3;
    }
    else if (grau == 3)
    {
        cont_grau = 7;
    }
    else
    {
        cont_grau = 1;
    }

    // faz uma insercao de grau cont_grau
    for (i = 0; i < tamanho; i++)
    { // insere uma relação em i com uma familia aleatoria
        for (j = 0; j < cont_grau; j++)
        {
            k = aleatorio(tamanho);
            // escolhe o tipo de relação que familia i tem com familia k
            l = aleatorio(3); // 0 definido como neutro, 1 define como positivo; 2 define como negativo
            if (l == 1)
            {
                // define o nivel da relação sendo positibva
                l = aleatorio(2);
            }
            else if (l == 2)
            {
                // define o nivel da relação sendo negativa
                l = 0 - aleatorio(2);
            }
            if (i == k)
            { // relação da familia com ela mesma, ocorrera um aumento drastico dentro de nivel
                relacoes[i][k] = relacoes[k][i] = 10 * l;
            }
            else
            { // esta eh a porcentagem da relação-> mais pra frente sera transformada em porcentagem(float menor que 1)
                relacoes[i][k] = relacoes[k][i] = l;
            }
        }
    }
    numero_pessoas = 0;
    for (i = 0; i < tamanho; i++)
    {
        numero_pessoas = familias[i] + numero_pessoas;
    }
}

void vida()
{
    int i, j, k, l;
    #pragma omp parallel for num_threads(MAX_THREADS)
    for (tempo = 0; tempo < 100000; tempo++)
    {
        porcentagem_real = (double)(ideal / numero_pessoas);
        nascidos = taxacrescimento;
        // confere se a população ainda lembra do desastre --> se esqueceu cont_esqucimento = 0 ->decrementador_desastre = 0
        // se nao esqueceu -->cont_esqucimento --
        // se a população esqueceu do desastre-> não ocorre mais decremento na população devido ao desastre
        if (desastre != 1)
        {
            decrementador_desastre = 0;
        }
        else
        {

            if (cont_esqucimento == 0)
            {
                decrementador_desastre = 0;
            }
            else
            {
                cont_esqucimento--;
            }
            // teste para ver se ocorre desastre natural
            probabilidade_de_desastre_natural = aleatorio(10000); // possibilidade muito minuscula para ocorrencia

            if (probabilidade_de_desastre_natural == 99)
            {   // ocorreu o desastre
                // o desastre vai de 1 a 9, sendo 1 um nível não preocupante e 9 um nível muito grave
                // rand para decidir qual o desastre
                probabilidade_de_desastre_natural = aleatorio(45001); // mesmo que ocorra o desastre --> ainda tem uma chance de ele ser tao fraco que se pode ignorar ele
                if (probabilidade_de_desastre_natural < 10000)
                {
                    decrementador_desastre = 1;
                }
                else if (probabilidade_de_desastre_natural < 18000 && probabilidade_de_desastre_natural > 9999)
                {
                    decrementador_desastre = 2;
                }
                else if (probabilidade_de_desastre_natural < 25000 && probabilidade_de_desastre_natural > 17999)
                {
                    decrementador_desastre = 3;
                }
                else if (probabilidade_de_desastre_natural < 31000 && probabilidade_de_desastre_natural > 24999)
                {
                    decrementador_desastre = 4;
                }
                else if (probabilidade_de_desastre_natural < 36000 && probabilidade_de_desastre_natural > 30999)
                {
                    decrementador_desastre = 5;
                }
                else if (probabilidade_de_desastre_natural < 40000 && probabilidade_de_desastre_natural > 35999)
                {
                    decrementador_desastre = 6;
                }
                else if (probabilidade_de_desastre_natural < 43000 && probabilidade_de_desastre_natural > 39999)
                {
                    decrementador_desastre = 7;
                }
                else if (probabilidade_de_desastre_natural < 45000 && probabilidade_de_desastre_natural > 42999)
                {
                    decrementador_desastre = 8;
                }
                else if (probabilidade_de_desastre_natural == 45000)
                {
                    decrementador_desastre = 9;
                }
                else
                    decrementador = 0;
                for (l = 0; l < tamanho; l++)
                {
                    if (decrementador_desastre == 9)
                        familias[l] = familias[l] * (0.999);
                    else
                        familias[l] = familias[l] * (1 - decrementador_desastre / 10); // população diminui de tamanho--> pessoas morrem
                }

                cont_esqucimento = (decrementador_desastre)*15; // define o tempo que a população se lembrara deste desastre
                printf("\n- Ocorreu um desastre de nivel %f no ano %d\n", decrementador_desastre, tempo);

                if (decrementador_desastre == 9)
                    printf("\n-- A humanidade chegou perto do apocalipse no ano %d\n", tempo);
            }
        }

        decrementador = (int)numero_pessoas * 0.002; // por ciclo 0.6% da população morre

        if (decrementador < 1){
          decrementador = 1;
        }

        // teste de incremento e decremento
        // se o rand for aceito --> ocorre o incremento ou decremento
        // taxacrescimento eh o numero maximo de pessoas que a população pode aumentar

        for (i = 0; i < tamanho; i++)
        {
            p = 1;
            cont_mudanca = 0;
            #pragma omp parallel for reduction(+: cont_mudanca) num_threads(MAX_THREADS)
            for (l = 0; l < tamanho; l++)
            {
                cont_mudanca = cont_mudanca + ((double)relacoes[i][l]) / 100; // transforma a matriz em porcentagem
            }
            // define a porcentagem de incremento
            incrementador = (porcentagem_real * porcentagem_real * (1 + porcentagens_paises) * (1 - decrementador_desastre / 10) * (1 + cont_mudanca)) * 0.01;

            // esta parte foi colocada pois o programa não estava aceitando double para numeros menores que 0.000001
            if (incrementador == 0)
                p = 10000000;
            else
            {
                while (incrementador < 1)
                {
                    incrementador = incrementador * 10;
                    p = p * 10;
                }
            }

            j = aleatorio(p);
            if (j < (int)incrementador && nascidos > 0)
            { // ocorre o incremento de 1 pessoa na população
                familias[i] = familias[i] + 1;
                nascidos--;
            }
            if (decrementador > 1)
            {
                tt = aleatorio(decrementador);
                if (tt == 1)
                {
                    familias[i]--;
                    decrementador--;
                }
            }
        }

        i = 0;
        while (decrementador > 0)
        {
            familias[i]--;
            i++;
            decrementador--;
        }

        // faz a atualização do numero total de pessoas na sociedade
        numero_pessoas = 0;

        #pragma omp parallel for private(i) reduction(+: numero_pessoas) num_threads(MAX_THREADS)
        for (i = 0; i < tamanho; i++)
        {
            numero_pessoas = familias[i] + numero_pessoas;
        }

        // teste para mudar o tipo de relação atual---> grafo Erdös-Rényi
        for (i = 0; i < tamanho; i++)
        {
            j = aleatorio(100000000); // se der certo, pode ocorrer a mudança de relação entre familias.

            if (j == 1)
            {
                for (j = 0; j < cont_grau; j++)
                {
                    k = aleatorio(tamanho);
                    // agora escolhe o tipo de relação que familia i tem com familia k
                    l = aleatorio(3); // 0 definido como neutro, 1 define como positivo; 2 define como negativo
                    if (l == 1)
                    {
                        // define o nivel da relação sendo 9 muito bom
                        l = aleatorio(2);
                    }
                    else if (l == 2)
                    {
                        // define o nivel da relação sendo -9 muito ruim
                        l = 0 - aleatorio(2);
                    }
                    if (i == k)
                    { // relação da familia com ela mesma, ocorrera um aumento drastico dentro de nivel
                        relacoes[i][k] = relacoes[k][i] = 10 * l;
                    }
                    else
                    {
                        relacoes[i][k] = relacoes[k][i] = l;
                    }
                }
            }
        }

        // teste para ver se a população passou do limite aceitavel
        // teste final: mesmo com os decrementadores, ocorre super população
        if ((1 / porcentagem_real) > 3)
        {
            //printf("\n\n\n\n\n atualmente existem %d neste pais\n", (int)numero_pessoas);
            // população gigante leva a falta de comida --> lutas por comidas -> leva a guerra --> diminuição da população
            numero_pessoas = 0;

            #pragma omp parallel for private(i) reduction(+: numero_pessoas) num_threads(MAX_THREADS)
            for (i = 0; i < tamanho; i++)
            {
                familias[i] = familias[i] * 0.1;
                if (familias[i] < 1)
                    familias[i] = 1;
                numero_pessoas = familias[i] + numero_pessoas;
            }

            printf("\n- Atualmente existem %d pessoas neste pais\n", (int)numero_pessoas);
        }
    }
}

int main()
{
    Duracao *valor;
    struct timeval start, end;

    inicializador();
    valor = tempo_decorrido(&start, &end);


    relacoes_paises();
    valor = tempo_decorrido(&start, &end);

    graus();
    valor = tempo_decorrido(&start, &end);

    printf("\nSe voce deseja que o programa funcione com desastres naturais ativos, digite 1, senao digite outro valor qualquer:\n");
    scanf("%d", &desastre);

    gettimeofday(&start, NULL);
    vida();
    gettimeofday(&end, NULL);
    valor = tempo_decorrido(&start, &end);
    printf("\nTempo para calcular a populacao final: %d,%d s\n", valor->secs, valor->usecs);

    return 0;
}
