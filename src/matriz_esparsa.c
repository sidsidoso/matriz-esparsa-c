#include "matriz_esparsa.h"

#include <stdint.h>
#include <stdlib.h>

typedef struct No {
    size_t linha;
    size_t coluna;
    double valor;
    struct No *direita; /* proximo nao zero da mesma linha  */
    struct No *abaixo;  /* proximo nao zero da mesma coluna */
} No;

struct MatrizEsparsa {
    size_t linhas;
    size_t colunas;
    size_t nao_zeros;
    No **inicio_linha;  /* inicio_linha[i]: primeiro no da linha i  */
    No **inicio_coluna; /* inicio_coluna[j]: primeiro no da coluna j */
};

/* ---------------------------------------------------------------------- */
/* Criacao e destruicao */

MatrizEsparsa *matriz_criar(size_t linhas, size_t colunas)
{
    if (linhas == 0 || colunas == 0) {
        return NULL;
    }
    MatrizEsparsa *m = malloc(sizeof *m);
    if (m == NULL) {
        return NULL;
    }
    m->linhas = linhas;
    m->colunas = colunas;
    m->nao_zeros = 0;
    /* calloc zera os ponteiros e confere estouro de linhas * sizeof. */
    m->inicio_linha = calloc(linhas, sizeof *m->inicio_linha);
    m->inicio_coluna = calloc(colunas, sizeof *m->inicio_coluna);
    if (m->inicio_linha == NULL || m->inicio_coluna == NULL) {
        free(m->inicio_linha);
        free(m->inicio_coluna);
        free(m);
        return NULL;
    }
    return m;
}

void matriz_destruir(MatrizEsparsa *m)
{
    if (m == NULL) {
        return;
    }
    /* Cada no esta em exatamente uma lista de linha: basta percorre-las. */
    for (size_t i = 0; i < m->linhas; i++) {
        No *no = m->inicio_linha[i];
        while (no != NULL) {
            No *prox = no->direita;
            free(no);
            no = prox;
        }
    }
    free(m->inicio_linha);
    free(m->inicio_coluna);
    free(m);
}

size_t matriz_linhas(const MatrizEsparsa *m)
{
    return m == NULL ? 0 : m->linhas;
}

size_t matriz_colunas(const MatrizEsparsa *m)
{
    return m == NULL ? 0 : m->colunas;
}

size_t matriz_nao_zeros(const MatrizEsparsa *m)
{
    return m == NULL ? 0 : m->nao_zeros;
}

/* ---------------------------------------------------------------------- */
/* Acesso individual */

bool matriz_definir(MatrizEsparsa *m, size_t i, size_t j, double valor)
{
    if (m == NULL || i >= m->linhas || j >= m->colunas) {
        return false;
    }
    /* Ponteiro para ponteiro: 'elo_l' e o campo que aponta (ou apontaria)
     * para o no (i, j) na lista da linha, sem caso especial para o inicio. */
    No **elo_l = &m->inicio_linha[i];
    while (*elo_l != NULL && (*elo_l)->coluna < j) {
        elo_l = &(*elo_l)->direita;
    }
    No *existente = (*elo_l != NULL && (*elo_l)->coluna == j) ? *elo_l : NULL;

    if (existente != NULL && valor != 0.0) {
        existente->valor = valor; /* so atualiza: nao mexe na coluna */
        return true;
    }
    if (existente == NULL && valor == 0.0) {
        return true; /* ja era zero */
    }

    No **elo_c = &m->inicio_coluna[j];
    while (*elo_c != NULL && (*elo_c)->linha < i) {
        elo_c = &(*elo_c)->abaixo;
    }

    if (existente != NULL) {
        /* valor == 0: desliga o no das duas listas. */
        *elo_l = existente->direita;
        *elo_c = existente->abaixo;
        free(existente);
        m->nao_zeros--;
        return true;
    }

    No *novo = malloc(sizeof *novo);
    if (novo == NULL) {
        return false;
    }
    novo->linha = i;
    novo->coluna = j;
    novo->valor = valor;
    novo->direita = *elo_l;
    novo->abaixo = *elo_c;
    *elo_l = novo;
    *elo_c = novo;
    m->nao_zeros++;
    return true;
}

double matriz_obter(const MatrizEsparsa *m, size_t i, size_t j)
{
    if (m == NULL || i >= m->linhas || j >= m->colunas) {
        return 0.0;
    }
    for (const No *no = m->inicio_linha[i]; no != NULL && no->coluna <= j;
         no = no->direita) {
        if (no->coluna == j) {
            return no->valor;
        }
    }
    return 0.0;
}

/* ---------------------------------------------------------------------- */
/* Construcao sequencial: as operacoes geram o resultado em ordem de linhas
 * (e, dentro da linha, de colunas crescentes). Assim cada no novo vai para o
 * FIM da sua linha e da sua coluna, em O(1), sem nenhuma busca. */

typedef struct {
    MatrizEsparsa *m;
    No **ultimo_coluna; /* ultimo no ja anexado em cada coluna */
    No **elo_linha;     /* onde o proximo no da linha atual sera ligado */
} Construtor;

static bool construtor_iniciar(Construtor *c, size_t linhas, size_t colunas)
{
    c->m = matriz_criar(linhas, colunas);
    c->ultimo_coluna = calloc(colunas, sizeof *c->ultimo_coluna);
    c->elo_linha = NULL;
    if (c->m == NULL || c->ultimo_coluna == NULL) {
        matriz_destruir(c->m);
        free(c->ultimo_coluna);
        return false;
    }
    return true;
}

static void construtor_nova_linha(Construtor *c, size_t i)
{
    c->elo_linha = &c->m->inicio_linha[i];
}

static bool construtor_anexar(Construtor *c, size_t i, size_t j, double valor)
{
    No *no = malloc(sizeof *no);
    if (no == NULL) {
        return false;
    }
    no->linha = i;
    no->coluna = j;
    no->valor = valor;
    no->direita = NULL;
    no->abaixo = NULL;

    *c->elo_linha = no;
    c->elo_linha = &no->direita;

    if (c->ultimo_coluna[j] != NULL) {
        c->ultimo_coluna[j]->abaixo = no;
    } else {
        c->m->inicio_coluna[j] = no;
    }
    c->ultimo_coluna[j] = no;
    c->m->nao_zeros++;
    return true;
}

/* Finaliza: devolve a matriz pronta (ou NULL e libera tudo se 'ok' for false). */
static MatrizEsparsa *construtor_finalizar(Construtor *c, bool ok)
{
    free(c->ultimo_coluna);
    if (!ok) {
        matriz_destruir(c->m);
        return NULL;
    }
    return c->m;
}

/* ---------------------------------------------------------------------- */
/* Operacoes */

MatrizEsparsa *matriz_somar(const MatrizEsparsa *a, const MatrizEsparsa *b)
{
    if (a == NULL || b == NULL || a->linhas != b->linhas || a->colunas != b->colunas) {
        return NULL;
    }
    Construtor c;
    if (!construtor_iniciar(&c, a->linhas, a->colunas)) {
        return NULL;
    }
    bool ok = true;
    for (size_t i = 0; i < a->linhas && ok; i++) {
        construtor_nova_linha(&c, i);
        /* Intercalacao de duas listas ordenadas por coluna (como no merge). */
        const No *p = a->inicio_linha[i];
        const No *q = b->inicio_linha[i];
        while ((p != NULL || q != NULL) && ok) {
            size_t j;
            double soma;
            if (q == NULL || (p != NULL && p->coluna < q->coluna)) {
                j = p->coluna; soma = p->valor; p = p->direita;
            } else if (p == NULL || q->coluna < p->coluna) {
                j = q->coluna; soma = q->valor; q = q->direita;
            } else {
                j = p->coluna; soma = p->valor + q->valor;
                p = p->direita; q = q->direita;
            }
            if (soma != 0.0) {
                ok = construtor_anexar(&c, i, j, soma);
            }
        }
    }
    return construtor_finalizar(&c, ok);
}

MatrizEsparsa *matriz_transpor(const MatrizEsparsa *m)
{
    if (m == NULL) {
        return NULL;
    }
    Construtor c;
    if (!construtor_iniciar(&c, m->colunas, m->linhas)) {
        return NULL;
    }
    bool ok = true;
    /* A coluna j de m, lida de cima para baixo, e exatamente a linha j da
     * transposta em ordem crescente de coluna: nao precisa ordenar nada. */
    for (size_t j = 0; j < m->colunas && ok; j++) {
        construtor_nova_linha(&c, j);
        for (const No *no = m->inicio_coluna[j]; no != NULL && ok; no = no->abaixo) {
            ok = construtor_anexar(&c, j, no->linha, no->valor);
        }
    }
    return construtor_finalizar(&c, ok);
}

static int comparar_size_t(const void *x, const void *y)
{
    size_t a = *(const size_t *)x, b = *(const size_t *)y;
    return (a > b) - (a < b);
}

MatrizEsparsa *matriz_multiplicar(const MatrizEsparsa *a, const MatrizEsparsa *b)
{
    if (a == NULL || b == NULL || a->colunas != b->linhas) {
        return NULL;
    }
    size_t n_col = b->colunas;
    /* Acumulador esparso (Gustavson):
     *  acumulado[j] - soma parcial de c[i][j];
     *  marca[j]     - i + 1 se a coluna j ja foi tocada na linha i;
     *  tocadas[]    - colunas tocadas na linha atual (para nao varrer c). */
    double *acumulado = malloc(n_col * sizeof *acumulado);
    size_t *marca = calloc(n_col, sizeof *marca);
    size_t *tocadas = malloc(n_col * sizeof *tocadas);
    Construtor c;
    bool ok = acumulado != NULL && marca != NULL && tocadas != NULL &&
              construtor_iniciar(&c, a->linhas, n_col);
    if (!ok) {
        free(acumulado);
        free(marca);
        free(tocadas);
        return NULL;
    }

    for (size_t i = 0; i < a->linhas && ok; i++) {
        size_t qtd = 0;
        for (const No *pa = a->inicio_linha[i]; pa != NULL; pa = pa->direita) {
            /* a[i][k] * (linha k de b) contribui para a linha i de c. */
            for (const No *pb = b->inicio_linha[pa->coluna]; pb != NULL; pb = pb->direita) {
                size_t j = pb->coluna;
                if (marca[j] != i + 1) {
                    marca[j] = i + 1;
                    acumulado[j] = 0.0;
                    tocadas[qtd++] = j;
                }
                acumulado[j] += pa->valor * pb->valor;
            }
        }
        qsort(tocadas, qtd, sizeof *tocadas, comparar_size_t);
        construtor_nova_linha(&c, i);
        for (size_t t = 0; t < qtd && ok; t++) {
            size_t j = tocadas[t];
            if (acumulado[j] != 0.0) { /* cancelamentos nao viram no */
                ok = construtor_anexar(&c, i, j, acumulado[j]);
            }
        }
    }
    free(acumulado);
    free(marca);
    free(tocadas);
    return construtor_finalizar(&c, ok);
}

bool matriz_multiplicar_vetor(const MatrizEsparsa *m, const double *x, double *y)
{
    if (m == NULL || x == NULL || y == NULL) {
        return false;
    }
    /* Sobreposicao comparada por endereco inteiro (comparar ponteiros de
     * objetos diferentes com < seria comportamento indefinido). */
    uintptr_t xi = (uintptr_t)x, yi = (uintptr_t)y;
    uintptr_t x_fim = xi + m->colunas * sizeof *x;
    uintptr_t y_fim = yi + m->linhas * sizeof *y;
    if (xi < y_fim && yi < x_fim) {
        return false;
    }
    for (size_t i = 0; i < m->linhas; i++) {
        double soma = 0.0;
        for (const No *no = m->inicio_linha[i]; no != NULL; no = no->direita) {
            soma += no->valor * x[no->coluna];
        }
        y[i] = soma;
    }
    return true;
}

/* ---------------------------------------------------------------------- */
/* Conversao e impressao */

MatrizEsparsa *matriz_de_densa(const double *dados, size_t linhas, size_t colunas)
{
    if (dados == NULL) {
        return NULL;
    }
    Construtor c;
    if (!construtor_iniciar(&c, linhas, colunas)) {
        return NULL;
    }
    bool ok = true;
    for (size_t i = 0; i < linhas && ok; i++) {
        construtor_nova_linha(&c, i);
        for (size_t j = 0; j < colunas && ok; j++) {
            double v = dados[i * colunas + j];
            if (v != 0.0) {
                ok = construtor_anexar(&c, i, j, v);
            }
        }
    }
    return construtor_finalizar(&c, ok);
}

bool matriz_para_densa(const MatrizEsparsa *m, double *dados)
{
    if (m == NULL || dados == NULL) {
        return false;
    }
    for (size_t i = 0; i < m->linhas; i++) {
        double *linha = dados + i * m->colunas;
        for (size_t j = 0; j < m->colunas; j++) {
            linha[j] = 0.0;
        }
        for (const No *no = m->inicio_linha[i]; no != NULL; no = no->direita) {
            linha[no->coluna] = no->valor;
        }
    }
    return true;
}

void matriz_imprimir_densa(const MatrizEsparsa *m, FILE *saida)
{
    if (m == NULL || saida == NULL) {
        return;
    }
    for (size_t i = 0; i < m->linhas; i++) {
        /* Um cursor anda pela lista da linha junto com j: O(c) por linha. */
        const No *no = m->inicio_linha[i];
        fputs("  [", saida);
        for (size_t j = 0; j < m->colunas; j++) {
            double v = 0.0;
            if (no != NULL && no->coluna == j) {
                v = no->valor;
                no = no->direita;
            }
            fprintf(saida, " %7.2f", v);
        }
        fputs(" ]\n", saida);
    }
}

void matriz_imprimir_esparsa(const MatrizEsparsa *m, FILE *saida)
{
    if (m == NULL || saida == NULL) {
        return;
    }
    for (size_t i = 0; i < m->linhas; i++) {
        for (const No *no = m->inicio_linha[i]; no != NULL; no = no->direita) {
            fprintf(saida, "  (%zu, %zu) = %g\n", no->linha, no->coluna, no->valor);
        }
    }
}

bool matriz_verificar(const MatrizEsparsa *m)
{
    if (m == NULL) {
        return false;
    }
    /* esperado[j] = proximo no que a lista da coluna j deve conter. Varrendo
     * as linhas em ordem, os nos da coluna j aparecem de cima para baixo;
     * entao cada no visto pela linha tem de ser exatamente esperado[j]. Isso
     * prova que linhas e colunas encadeiam os MESMOS nos, na ordem certa. */
    const No **esperado = malloc(m->colunas * sizeof *esperado);
    if (esperado == NULL) {
        return false;
    }
    for (size_t j = 0; j < m->colunas; j++) {
        esperado[j] = m->inicio_coluna[j];
    }
    bool ok = true;
    size_t contados = 0;
    for (size_t i = 0; i < m->linhas && ok; i++) {
        const No *ant = NULL;
        for (const No *no = m->inicio_linha[i]; no != NULL && ok; no = no->direita) {
            ok = no->linha == i && no->coluna < m->colunas && no->valor != 0.0 &&
                 (ant == NULL || ant->coluna < no->coluna) &&
                 esperado[no->coluna] == no;
            if (ok) {
                esperado[no->coluna] = no->abaixo;
                ant = no;
                contados++;
            }
        }
    }
    for (size_t j = 0; j < m->colunas && ok; j++) {
        ok = esperado[j] == NULL; /* nenhuma coluna tem nos "sobrando" */
    }
    free(esperado);
    return ok && contados == m->nao_zeros;
}
