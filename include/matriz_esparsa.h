#ifndef MATRIZ_ESPARSA_H
#define MATRIZ_ESPARSA_H

/*
 * Matriz esparsa de double com listas encadeadas ortogonais.
 *
 * Apenas os elementos nao nulos sao armazenados. Cada no pertence a duas
 * listas ao mesmo tempo:
 *   - a lista da sua LINHA, ordenada por coluna crescente (ponteiro 'direita');
 *   - a lista da sua COLUNA, ordenada por linha crescente (ponteiro 'abaixo').
 * Indices sao base 0. Definir um elemento como 0.0 remove o no.
 *
 * Notacao de custo: l = linhas, c = colunas, nnz = nao zeros,
 * nnz_i = nao zeros da linha i, nnz_j = nao zeros da coluna j.
 */

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

/* Tipo opaco: nos e cabecalhos ficam escondidos em src/matriz_esparsa.c. */
typedef struct MatrizEsparsa MatrizEsparsa;

/* Cria uma matriz l x c toda nula. Retorna NULL se l ou c forem 0 ou se
 * faltar memoria. O(l + c) */
MatrizEsparsa *matriz_criar(size_t linhas, size_t colunas);

/* Libera todos os nos e a matriz. Aceita NULL. O(nnz + l + c) */
void matriz_destruir(MatrizEsparsa *m);

/* Dimensoes e quantidade de nao zeros. O(1) */
size_t matriz_linhas(const MatrizEsparsa *m);
size_t matriz_colunas(const MatrizEsparsa *m);
size_t matriz_nao_zeros(const MatrizEsparsa *m);

/* Define m[i][j] = valor. Valor 0.0 remove o no (se existir).
 * Retorna false se (i, j) estiver fora da matriz ou faltar memoria.
 * O(nnz_i + nnz_j) */
bool matriz_definir(MatrizEsparsa *m, size_t i, size_t j, double valor);

/* Le m[i][j] (0.0 se o elemento nao estiver armazenado ou (i, j) for
 * invalido). O(nnz_i) */
double matriz_obter(const MatrizEsparsa *m, size_t i, size_t j);

/* Retorna uma NOVA matriz a + b, ou NULL se as dimensoes forem diferentes
 * ou faltar memoria. Somas que se anulam nao geram no. O(nnz_a + nnz_b + l + c) */
MatrizEsparsa *matriz_somar(const MatrizEsparsa *a, const MatrizEsparsa *b);

/* Retorna uma NOVA matriz a * b (a: l x n, b: n x c), ou NULL se
 * colunas(a) != linhas(b) ou faltar memoria. Algoritmo de Gustavson,
 * linha a linha com acumulador esparso.
 * O(l + c + f + soma_i r_i log r_i), onde f = multiplicacoes escalares
 * efetivamente feitas e r_i = nao zeros da linha i do resultado. */
MatrizEsparsa *matriz_multiplicar(const MatrizEsparsa *a, const MatrizEsparsa *b);

/* Retorna uma NOVA matriz transposta, ou NULL se faltar memoria.
 * O(nnz + l + c) */
MatrizEsparsa *matriz_transpor(const MatrizEsparsa *m);

/* y = m * x, com x de tamanho colunas(m) e y de tamanho linhas(m).
 * x e y nao podem se sobrepor. Retorna false se algum ponteiro for NULL
 * ou se x e y se sobrepuserem. O(nnz + l) */
bool matriz_multiplicar_vetor(const MatrizEsparsa *m, const double *x, double *y);

/* Monta uma matriz a partir de um vetor denso em ordem de linhas
 * (dados[i * colunas + j]). Retorna NULL em erro. O(l * c) */
MatrizEsparsa *matriz_de_densa(const double *dados, size_t linhas, size_t colunas);

/* Escreve a matriz em 'dados' (l * c posicoes, ordem de linhas).
 * Retorna false se algum ponteiro for NULL. O(l * c) */
bool matriz_para_densa(const MatrizEsparsa *m, double *dados);

/* Imprime a matriz completa (inclusive os zeros), uma linha por linha.
 * O(l * c) */
void matriz_imprimir_densa(const MatrizEsparsa *m, FILE *saida);

/* Imprime apenas os nao zeros como "(i, j) = valor", em ordem de linhas.
 * O(nnz + l) */
void matriz_imprimir_esparsa(const MatrizEsparsa *m, FILE *saida);

/* Confere todas as invariantes da estrutura (listas ordenadas, cada no na
 * linha e na coluna certas, contagem de nao zeros, nenhum zero armazenado).
 * Usada nos testes. O(nnz + l + c) */
bool matriz_verificar(const MatrizEsparsa *m);

#endif /* MATRIZ_ESPARSA_H */
