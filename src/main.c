#include <stdio.h>
#include <stdlib.h>

#include "matriz_esparsa.h"

static void mostrar(const char *titulo, const MatrizEsparsa *m)
{
    size_t total = matriz_linhas(m) * matriz_colunas(m);
    printf("%s (%zu x %zu, %zu nao zeros, densidade %.1f%%):\n", titulo,
           matriz_linhas(m), matriz_colunas(m), matriz_nao_zeros(m),
           100.0 * (double)matriz_nao_zeros(m) / (double)total);
    matriz_imprimir_densa(m, stdout);
}

int main(void)
{
    int status = EXIT_FAILURE;
    MatrizEsparsa *a = NULL, *b = NULL, *t = NULL, *s = NULL, *p = NULL, *g = NULL;

    /* Matriz 4 x 5 montada elemento a elemento. */
    a = matriz_criar(4, 5);
    if (a == NULL) {
        goto fim;
    }
    struct { size_t i, j; double v; } elementos[] = {
        {0, 0, 5.0}, {0, 3, -1.5}, {1, 2, 2.0}, {2, 4, 7.25}, {3, 1, 3.0}, {3, 3, 4.0},
    };
    for (size_t k = 0; k < sizeof elementos / sizeof elementos[0]; k++) {
        if (!matriz_definir(a, elementos[k].i, elementos[k].j, elementos[k].v)) {
            goto fim;
        }
    }

    puts("== Matriz esparsa com listas ortogonais ==\n");
    mostrar("A", a);
    puts("A em triplas (so os nao zeros):");
    matriz_imprimir_esparsa(a, stdout);
    printf("obter A[2][4] = %g, obter A[2][3] = %g\n\n",
           matriz_obter(a, 2, 4), matriz_obter(a, 2, 3));

    t = matriz_transpor(a);
    if (t == NULL) {
        goto fim;
    }
    mostrar("A^T (transposta)", t);
    putchar('\n');

    const double dados_b[] = {
        1, 0, 0, 1.5, 0,
        0, 0, -2, 0, 0,
        0, 0, 0, 0, -7.25,
        0, 1, 0, 0, 0,
    };
    b = matriz_de_densa(dados_b, 4, 5);
    if (b == NULL) {
        goto fim;
    }
    mostrar("B", b);
    putchar('\n');
    s = matriz_somar(a, b);
    if (s == NULL) {
        goto fim;
    }
    mostrar("A + B", s);
    puts("  (as posicoes (0,3), (1,2) e (2,4) se anulam: nenhum no e criado)");
    putchar('\n');

    p = matriz_multiplicar(a, t);
    if (p == NULL) {
        goto fim;
    }
    mostrar("A * A^T", p);
    putchar('\n');

    const double x[] = {1, 2, 3, 4, 5};
    double y[4];
    if (!matriz_multiplicar_vetor(a, x, y)) {
        goto fim;
    }
    printf("A * x, com x = [1 2 3 4 5]: y = [");
    for (size_t i = 0; i < 4; i++) {
        printf(" %g", y[i]);
    }
    puts(" ]\n");

    matriz_definir(a, 0, 3, 0.0);
    matriz_definir(a, 2, 4, 0.0);
    mostrar("A apos definir A[0][3] = 0 e A[2][4] = 0", a);
    putchar('\n');

    /* Matriz grande e muito esparsa: tridiagonal 100000 x 100000. */
    const size_t n = 100000;
    g = matriz_criar(n, n);
    if (g == NULL) {
        goto fim;
    }
    for (size_t i = 0; i < n; i++) {
        if (!matriz_definir(g, i, i, 2.0) ||
            (i > 0 && !matriz_definir(g, i, i - 1, -1.0)) ||
            (i + 1 < n && !matriz_definir(g, i, i + 1, -1.0))) {
            goto fim;
        }
    }
    size_t nnz = matriz_nao_zeros(g);
    /* No = linha + coluna + valor + 2 ponteiros; mais os 2 vetores de cabecalhos. */
    double bytes_no = 2.0 * sizeof(size_t) + sizeof(double) + 2.0 * sizeof(void *);
    double mib = 1024.0 * 1024.0;
    printf("tridiagonal %zu x %zu: %zu nao zeros\n", n, n, nnz);
    printf("  armazenamento denso: %.0f MiB | listas ortogonais: ~%.1f MiB\n",
           (double)n * (double)n * sizeof(double) / mib,
           ((double)nnz * bytes_no + 2.0 * (double)n * sizeof(void *)) / mib);

    status = EXIT_SUCCESS;
fim:
    if (status != EXIT_SUCCESS) {
        fprintf(stderr, "erro: memoria insuficiente\n");
    }
    matriz_destruir(a);
    matriz_destruir(b);
    matriz_destruir(t);
    matriz_destruir(s);
    matriz_destruir(p);
    matriz_destruir(g);
    return status;
}
