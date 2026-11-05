#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "matriz_esparsa.h"

static int testes_passaram = 0;

#define RODAR(teste)                          \
    do {                                      \
        teste();                              \
        testes_passaram++;                    \
        printf("  [ok] %s\n", #teste);        \
    } while (0)

/* ---------------------------------------------------------------------- */
/* Utilitarios: gerador deterministico e implementacao densa de referencia */

static uint32_t semente = 123456789u;
static uint32_t aleatorio(void)
{
    semente ^= semente << 13;
    semente ^= semente >> 17;
    semente ^= semente << 5;
    return semente;
}

/* Valor nao nulo: inteiro em [-9, 9] (contas exatas em double) ou real. */
static double valor_aleatorio(bool inteiro)
{
    if (inteiro) {
        int v = (int)(aleatorio() % 18) - 9;
        return v >= 0 ? v + 1 : v; /* nunca zero */
    }
    return ((double)(aleatorio() % 2000001) - 1000000.0) / 1000.0 + 0.0005;
}

/* Matriz densa l x c com cada posicao nao nula com probabilidade densidade%. */
static double *densa_aleatoria(size_t l, size_t c, unsigned densidade, bool inteiro)
{
    double *d = malloc(l * c * sizeof *d);
    assert(d != NULL);
    for (size_t k = 0; k < l * c; k++) {
        d[k] = (aleatorio() % 100 < densidade) ? valor_aleatorio(inteiro) : 0.0;
    }
    return d;
}

static size_t contar_densa(const double *d, size_t n)
{
    size_t nnz = 0;
    for (size_t k = 0; k < n; k++) {
        nnz += d[k] != 0.0;
    }
    return nnz;
}

static void multiplicar_densa(const double *a, const double *b, double *c,
                              size_t l, size_t n, size_t m)
{
    for (size_t i = 0; i < l; i++) {
        for (size_t j = 0; j < m; j++) {
            double s = 0.0;
            for (size_t k = 0; k < n; k++) {
                s += a[i * n + k] * b[k * m + j];
            }
            c[i * m + j] = s;
        }
    }
}

static bool quase_igual(double x, double y)
{
    return fabs(x - y) <= 1e-9 * (1.0 + fabs(x) + fabs(y));
}

/* Confere a matriz esparsa contra a densa: elemento a elemento (via obter e
 * via para_densa), contagem de nao zeros e invariantes internas. */
static void conferir(const MatrizEsparsa *m, const double *d, size_t l, size_t c,
                     bool exato)
{
    assert(m != NULL);
    assert(matriz_verificar(m));
    assert(matriz_linhas(m) == l && matriz_colunas(m) == c);
    double *copia = malloc(l * c * sizeof *copia);
    assert(copia != NULL);
    assert(matriz_para_densa(m, copia));
    for (size_t i = 0; i < l; i++) {
        for (size_t j = 0; j < c; j++) {
            double esperado = d[i * c + j];
            double v = matriz_obter(m, i, j);
            assert(v == copia[i * c + j]);
            assert(exato ? v == esperado : quase_igual(v, esperado));
        }
    }
    if (exato) {
        assert(matriz_nao_zeros(m) == contar_densa(d, l * c));
    }
    free(copia);
}

/* ---------------------------------------------------------------------- */

static void teste_criar(void)
{
    assert(matriz_criar(0, 3) == NULL);
    assert(matriz_criar(3, 0) == NULL);
    MatrizEsparsa *m = matriz_criar(3, 4);
    assert(m != NULL);
    assert(matriz_linhas(m) == 3 && matriz_colunas(m) == 4);
    assert(matriz_nao_zeros(m) == 0);
    assert(matriz_obter(m, 2, 3) == 0.0);
    assert(matriz_obter(m, 3, 0) == 0.0); /* fora da matriz */
    assert(matriz_verificar(m));
    matriz_destruir(m);
}

static void teste_definir_obter_remover(void)
{
    MatrizEsparsa *m = matriz_criar(4, 5);
    assert(m != NULL);
    assert(matriz_definir(m, 1, 3, 2.5));
    assert(matriz_definir(m, 1, 0, -1.0)); /* antes na mesma linha */
    assert(matriz_definir(m, 0, 3, 7.0));  /* acima na mesma coluna */
    assert(matriz_definir(m, 3, 3, 4.0));  /* abaixo na mesma coluna */
    assert(matriz_nao_zeros(m) == 4);
    assert(matriz_obter(m, 1, 3) == 2.5);
    assert(matriz_obter(m, 1, 0) == -1.0);
    assert(matriz_verificar(m));

    assert(matriz_definir(m, 1, 3, 9.0)); /* atualiza sem criar no */
    assert(matriz_nao_zeros(m) == 4);
    assert(matriz_obter(m, 1, 3) == 9.0);

    assert(matriz_definir(m, 1, 3, 0.0)); /* zero remove */
    assert(matriz_nao_zeros(m) == 3);
    assert(matriz_obter(m, 1, 3) == 0.0);
    assert(matriz_definir(m, 2, 2, 0.0)); /* zero onde ja era zero */
    assert(matriz_nao_zeros(m) == 3);
    assert(matriz_verificar(m));

    assert(!matriz_definir(m, 4, 0, 1.0));
    assert(!matriz_definir(m, 0, 5, 1.0));
    matriz_destruir(m);
}

static void teste_definir_aleatorio_contra_densa(void)
{
    enum { L = 13, C = 17 };
    double d[L * C] = {0};
    MatrizEsparsa *m = matriz_criar(L, C);
    assert(m != NULL);
    for (int op = 0; op < 5000; op++) {
        size_t i = aleatorio() % L, j = aleatorio() % C;
        double v = (aleatorio() % 3 == 0) ? 0.0 : valor_aleatorio(true);
        assert(matriz_definir(m, i, j, v));
        d[i * C + j] = v;
        if (op % 250 == 0) {
            conferir(m, d, L, C, true);
        }
    }
    conferir(m, d, L, C, true);
    /* Zerar tudo deve liberar todos os nos. */
    for (size_t i = 0; i < L; i++) {
        for (size_t j = 0; j < C; j++) {
            assert(matriz_definir(m, i, j, 0.0));
        }
    }
    assert(matriz_nao_zeros(m) == 0);
    assert(matriz_verificar(m));
    matriz_destruir(m);
}

static void teste_densa_ida_e_volta(void)
{
    const double d[] = {
        0, 0, 3,
        4, 0, 0,
        0, 0, 0,
        0, 5, 6,
    };
    MatrizEsparsa *m = matriz_de_densa(d, 4, 3);
    conferir(m, d, 4, 3, true);
    assert(matriz_nao_zeros(m) == 4);
    assert(matriz_de_densa(NULL, 2, 2) == NULL);
    matriz_destruir(m);
}

static void teste_somar(void)
{
    for (int rodada = 0; rodada < 40; rodada++) {
        size_t l = 1 + aleatorio() % 25, c = 1 + aleatorio() % 25;
        unsigned dens = 5 + aleatorio() % 40;
        bool inteiro = rodada % 2 == 0;
        double *da = densa_aleatoria(l, c, dens, inteiro);
        double *db = densa_aleatoria(l, c, dens, inteiro);
        double *ds = malloc(l * c * sizeof *ds);
        assert(ds != NULL);
        for (size_t k = 0; k < l * c; k++) {
            ds[k] = da[k] + db[k];
        }
        MatrizEsparsa *a = matriz_de_densa(da, l, c);
        MatrizEsparsa *b = matriz_de_densa(db, l, c);
        MatrizEsparsa *s = matriz_somar(a, b);
        conferir(s, ds, l, c, true); /* soma de dois doubles: exata nos dois casos */
        matriz_destruir(a);
        matriz_destruir(b);
        matriz_destruir(s);
        free(da);
        free(db);
        free(ds);
    }
}

static void teste_somar_cancelamento_e_dimensoes(void)
{
    enum { L = 9, C = 11 };
    double *d = densa_aleatoria(L, C, 40, false);
    double neg[L * C];
    for (size_t k = 0; k < L * C; k++) {
        neg[k] = -d[k];
    }
    MatrizEsparsa *a = matriz_de_densa(d, L, C);
    MatrizEsparsa *b = matriz_de_densa(neg, L, C);
    MatrizEsparsa *s = matriz_somar(a, b);
    assert(s != NULL && matriz_nao_zeros(s) == 0); /* A + (-A) nao guarda nos */
    assert(matriz_verificar(s));

    MatrizEsparsa *outra = matriz_criar(C, L);
    assert(matriz_somar(a, outra) == NULL);
    assert(matriz_somar(a, NULL) == NULL);
    matriz_destruir(outra);
    matriz_destruir(s);
    matriz_destruir(a);
    matriz_destruir(b);
    free(d);
}

static void teste_transpor(void)
{
    for (int rodada = 0; rodada < 30; rodada++) {
        size_t l = 1 + aleatorio() % 30, c = 1 + aleatorio() % 30;
        double *d = densa_aleatoria(l, c, 1 + aleatorio() % 50, false);
        double *dt = malloc(l * c * sizeof *dt);
        assert(dt != NULL);
        for (size_t i = 0; i < l; i++) {
            for (size_t j = 0; j < c; j++) {
                dt[j * l + i] = d[i * c + j];
            }
        }
        MatrizEsparsa *m = matriz_de_densa(d, l, c);
        MatrizEsparsa *t = matriz_transpor(m);
        conferir(t, dt, c, l, true);
        MatrizEsparsa *tt = matriz_transpor(t); /* (A^T)^T = A */
        conferir(tt, d, l, c, true);
        matriz_destruir(m);
        matriz_destruir(t);
        matriz_destruir(tt);
        free(d);
        free(dt);
    }
}

static void teste_multiplicar(void)
{
    for (int rodada = 0; rodada < 40; rodada++) {
        size_t l = 1 + aleatorio() % 20, n = 1 + aleatorio() % 20, m = 1 + aleatorio() % 20;
        unsigned dens = 3 + aleatorio() % 45;
        bool inteiro = rodada % 2 == 0;
        double *da = densa_aleatoria(l, n, dens, inteiro);
        double *db = densa_aleatoria(n, m, dens, inteiro);
        double *dc = malloc(l * m * sizeof *dc);
        assert(dc != NULL);
        multiplicar_densa(da, db, dc, l, n, m);

        MatrizEsparsa *a = matriz_de_densa(da, l, n);
        MatrizEsparsa *b = matriz_de_densa(db, n, m);
        MatrizEsparsa *c = matriz_multiplicar(a, b);
        /* Inteiros pequenos: produto exato. Reais: ordem de soma pode mudar
         * o ultimo bit, entao compara com tolerancia relativa. */
        conferir(c, dc, l, m, inteiro);
        matriz_destruir(a);
        matriz_destruir(b);
        matriz_destruir(c);
        free(da);
        free(db);
        free(dc);
    }
}

static void teste_multiplicar_casos_especiais(void)
{
    /* Identidade: I * A = A. */
    enum { N = 6 };
    MatrizEsparsa *id = matriz_criar(N, N);
    assert(id != NULL);
    for (size_t i = 0; i < N; i++) {
        assert(matriz_definir(id, i, i, 1.0));
    }
    double *d = densa_aleatoria(N, N, 30, false);
    MatrizEsparsa *a = matriz_de_densa(d, N, N);
    MatrizEsparsa *r = matriz_multiplicar(id, a);
    conferir(r, d, N, N, true);
    matriz_destruir(r);

    /* Cancelamento: [1 1] * [1; -1] = [0] nao gera no. */
    const double lin[] = {1, 1}, col[] = {1, -1};
    MatrizEsparsa *x = matriz_de_densa(lin, 1, 2);
    MatrizEsparsa *y = matriz_de_densa(col, 2, 1);
    r = matriz_multiplicar(x, y);
    assert(r != NULL && matriz_nao_zeros(r) == 0 && matriz_verificar(r));
    matriz_destruir(r);

    /* Dimensoes incompativeis. */
    assert(matriz_multiplicar(x, x) == NULL);
    assert(matriz_multiplicar(NULL, x) == NULL);

    matriz_destruir(x);
    matriz_destruir(y);
    matriz_destruir(a);
    matriz_destruir(id);
    free(d);
}

static void teste_multiplicar_vetor(void)
{
    for (int rodada = 0; rodada < 30; rodada++) {
        size_t l = 1 + aleatorio() % 40, c = 1 + aleatorio() % 40;
        double *d = densa_aleatoria(l, c, 20, false);
        double *x = malloc(c * sizeof *x);
        double *y = malloc(l * sizeof *y);
        assert(x != NULL && y != NULL);
        for (size_t j = 0; j < c; j++) {
            x[j] = valor_aleatorio(false);
        }
        MatrizEsparsa *m = matriz_de_densa(d, l, c);
        assert(matriz_multiplicar_vetor(m, x, y));
        for (size_t i = 0; i < l; i++) {
            double s = 0.0;
            for (size_t j = 0; j < c; j++) {
                s += d[i * c + j] * x[j];
            }
            assert(quase_igual(y[i], s));
        }
        assert(!matriz_multiplicar_vetor(m, NULL, y));
        assert(!matriz_multiplicar_vetor(m, x, NULL));
        matriz_destruir(m);
        free(d);
        free(x);
        free(y);
    }
    /* x e y sobrepostos sao recusados. */
    MatrizEsparsa *q = matriz_criar(3, 3);
    double v[3] = {1, 2, 3};
    assert(q != NULL && !matriz_multiplicar_vetor(q, v, v));
    matriz_destruir(q);
}

static void teste_matriz_grande_e_esparsa(void)
{
    /* 2000 x 2000 com ~0,1% de nao zeros: a versao densa teria 4 milhoes de
     * posicoes; aqui so ~4000 nos. A^T * A e conferido em alguns pontos. */
    enum { N = 2000 };
    MatrizEsparsa *a = matriz_criar(N, N);
    assert(a != NULL);
    for (size_t k = 0; k < 4000; k++) {
        assert(matriz_definir(a, aleatorio() % N, aleatorio() % N, valor_aleatorio(true)));
    }
    assert(matriz_verificar(a));
    MatrizEsparsa *t = matriz_transpor(a);
    MatrizEsparsa *p = matriz_multiplicar(t, a);
    assert(p != NULL && matriz_verificar(p));
    for (int k = 0; k < 200; k++) {
        size_t i = aleatorio() % N, j = aleatorio() % N;
        double s = 0.0;
        for (size_t r = 0; r < N; r++) {
            s += matriz_obter(a, r, i) * matriz_obter(a, r, j);
        }
        assert(matriz_obter(p, i, j) == s);
    }
    /* Diagonal de A^T A = soma dos quadrados da coluna. */
    for (size_t j = 0; j < N; j += 97) {
        double s = 0.0;
        for (size_t r = 0; r < N; r++) {
            double v = matriz_obter(a, r, j);
            s += v * v;
        }
        assert(matriz_obter(p, j, j) == s);
    }
    matriz_destruir(a);
    matriz_destruir(t);
    matriz_destruir(p);
}

static void teste_ponteiros_nulos(void)
{
    double v = 0;
    assert(matriz_linhas(NULL) == 0);
    assert(matriz_colunas(NULL) == 0);
    assert(matriz_nao_zeros(NULL) == 0);
    assert(!matriz_definir(NULL, 0, 0, 1.0));
    assert(matriz_obter(NULL, 0, 0) == 0.0);
    assert(matriz_transpor(NULL) == NULL);
    assert(!matriz_multiplicar_vetor(NULL, &v, &v));
    assert(!matriz_para_densa(NULL, &v));
    assert(!matriz_verificar(NULL));
    matriz_imprimir_densa(NULL, stdout);
    matriz_imprimir_esparsa(NULL, stdout);
    matriz_destruir(NULL);
}

int main(void)
{
    puts("Testes da matriz esparsa (comparando com implementacao densa):");
    RODAR(teste_criar);
    RODAR(teste_definir_obter_remover);
    RODAR(teste_definir_aleatorio_contra_densa);
    RODAR(teste_densa_ida_e_volta);
    RODAR(teste_somar);
    RODAR(teste_somar_cancelamento_e_dimensoes);
    RODAR(teste_transpor);
    RODAR(teste_multiplicar);
    RODAR(teste_multiplicar_casos_especiais);
    RODAR(teste_multiplicar_vetor);
    RODAR(teste_matriz_grande_e_esparsa);
    RODAR(teste_ponteiros_nulos);
    printf("%d testes passaram\n", testes_passaram);
    return 0;
}
