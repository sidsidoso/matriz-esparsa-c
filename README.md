# Matriz Esparsa com Listas Ortogonais em C

[![CI](https://github.com/sidsidoso/matriz-esparsa-c/actions/workflows/ci.yml/badge.svg)](https://github.com/sidsidoso/matriz-esparsa-c/actions/workflows/ci.yml)
![Linguagem](https://img.shields.io/badge/linguagem-C11-blue.svg)
[![Licença: MIT](https://img.shields.io/badge/licen%C3%A7a-MIT-green.svg)](LICENSE)

Implementação de uma **matriz esparsa** de `double` em C11 usando **listas
encadeadas ortogonais**: só os elementos não nulos ocupam memória, e cada um
deles está ao mesmo tempo na lista da sua linha e na lista da sua coluna.
Suporta definir/obter elementos (definir 0 remove o nó), soma, multiplicação
(algoritmo de Gustavson), transposição, produto matriz × vetor, impressão densa
e contagem de não zeros. Os testes comparam tudo com uma implementação densa
em matrizes aleatórias.
Projeto desenvolvido para a disciplina de **Estrutura de Dados**.

---

## Sumário

- [Conceito](#conceito)
- [Por que listas ortogonais e não CSR?](#por-que-listas-ortogonais-e-não-csr)
- [Operações e complexidade](#operações-e-complexidade)
- [Decisões de implementação e otimizações](#decisões-de-implementação-e-otimizações)
- [Estrutura de pastas](#estrutura-de-pastas)
- [Como compilar, executar e testar](#como-compilar-executar-e-testar)
- [Exemplo de saída](#exemplo-de-saída)
- [Referências](#referências)

---

## Conceito

Uma matriz é **esparsa** quando a maioria dos elementos é zero — por exemplo,
matrizes de adjacência de grafos, sistemas lineares de elementos finitos ou a
matriz tridiagonal do exemplo, em que 100 000 × 100 000 posições (≈ 76 GiB
densa) têm só ≈ 300 mil não zeros (≈ 13 MiB aqui).

Cada não zero vira um **nó** com suas coordenadas e dois ponteiros:

```c
typedef struct No {
    size_t linha, coluna;
    double valor;
    struct No *direita;  /* próximo não zero da mesma linha  */
    struct No *abaixo;   /* próximo não zero da mesma coluna */
} No;
```

A matriz guarda dois vetores de cabeçalhos: `inicio_linha[i]` e
`inicio_coluna[j]`. Cada lista de linha é ordenada por coluna e cada lista de
coluna é ordenada por linha:

```
 inicio_coluna:  [0]          [1]          [2]          [3]
                  |            |            |            |
 inicio_linha     v            |            |            v
   [0] -------> (0,0) 5 ---------------------------> (0,3) -1.5 --> NULL
                  |            |            v            |
   [1] ------------------------------> (1,2) 2 ----------------> NULL
                  |            |            |            |
   [2] --> NULL   |            |            |            |
                  |            v            |            v
   [3] ----------------> (3,1) 3 ------------------> (3,3) 4 -----> NULL
                  v            v            v            v
                 NULL         NULL         NULL         NULL

   setas horizontais = 'direita' (mesma linha)
   setas verticais   = 'abaixo'  (mesma coluna)

            [ 5  0  0  -1.5 ]
  densa:    [ 0  0  2   0   ]
            [ 0  0  0   0   ]
            [ 0  3  0   4   ]
```

**Definir `m[i][j]`**: percorre a linha *i* até a coluna *j* e a coluna *j*
até a linha *i* usando ponteiro para ponteiro (`No **elo`), então inserir ou
remover no início de uma lista não é caso especial:

```
 definir(1, 0, 9.0):  linha 1:  [1,2] ...        -> [1,0| 9.0] -> [1,2] ...
                      coluna 0: [0,0] -> NULL    -> [0,0] -> [1,0| 9.0] -> NULL
 definir(1, 0, 0.0):  desliga o nó das duas listas e libera (nnz diminui)
```

## Por que listas ortogonais e não CSR?

O formato **CSR** (*Compressed Sparse Row*) guarda três vetores (`valores`,
`colunas`, `inicio_linha`) e é imbatível para percorrer linhas em sequência
(cache e memória). Mas ele é **estático**: inserir ou remover um não zero
exige deslocar O(nnz) posições, e acessar uma **coluna** exige varrer a matriz
inteira (ou montar o CSC à parte).

Este projeto adota **listas ortogonais** porque:

| Critério                          | Listas ortogonais           | CSR                         |
|-----------------------------------|-----------------------------|-----------------------------|
| `definir` / remover um elemento   | O(nnz_i + nnz_j), sem deslocar nada | O(nnz) (desloca vetores) |
| Acesso a uma coluna inteira       | O(nnz_j) direto             | O(nnz) ou CSC extra         |
| Transposição                      | O(nnz + l + c), já ordenada | O(nnz + l + c) (contagem)   |
| Memória por não zero              | 40 bytes (64 bits)          | 16 bytes                    |
| Localidade de cache               | pior (nós espalhados)       | ótima                       |

Como o trabalho pede uma estrutura **mutável** (definir e zerar elementos a
qualquer momento) e operações que usam linhas **e** colunas (transposição),
as listas ortogonais são a escolha natural. Para recuperar parte da eficiência
do CSR nas operações em lote, os resultados de soma, multiplicação,
transposição e conversão de densa são **construídos em ordem** (ver abaixo),
sem nenhuma busca.

## Operações e complexidade

Notação: *l* linhas, *c* colunas, *nnz* não zeros, *nnz_i* / *nnz_j* não zeros
da linha *i* / coluna *j*.

| Operação                     | Função                      | Tempo                              | Espaço extra |
|------------------------------|-----------------------------|------------------------------------|--------------|
| Criar                        | `matriz_criar`              | O(l + c)                           | —            |
| Destruir                     | `matriz_destruir`           | O(nnz + l + c)                     | O(1)         |
| Definir (0 remove)           | `matriz_definir`            | O(nnz_i + nnz_j)                   | O(1)         |
| Obter                        | `matriz_obter`              | O(nnz_i)                           | O(1)         |
| Somar A + B                  | `matriz_somar`              | O(nnz_A + nnz_B + l + c)           | O(c)         |
| Multiplicar A · B            | `matriz_multiplicar`        | O(l + c + f + Σ rᵢ log rᵢ) ¹       | O(c)         |
| Transpor                     | `matriz_transpor`           | O(nnz + l + c)                     | O(l)         |
| Matriz × vetor               | `matriz_multiplicar_vetor`  | O(nnz + l)                         | O(1)         |
| Contar não zeros / dimensões | `matriz_nao_zeros`, ...     | O(1)                               | O(1)         |
| Imprimir densa               | `matriz_imprimir_densa`     | O(l · c)                           | O(1)         |
| Imprimir triplas             | `matriz_imprimir_esparsa`   | O(nnz + l)                         | O(1)         |
| Converter de/para densa      | `matriz_de_densa`, `matriz_para_densa` | O(l · c)                | O(c)         |
| Verificar invariantes        | `matriz_verificar`          | O(nnz + l + c)                     | O(c)         |

¹ *f* = número de multiplicações escalares realmente feitas (pares
`a[i][k] · b[k][j]` com ambos não nulos) e *rᵢ* = não zeros da linha *i* do
resultado. Para matrizes esparsas isso é muito menor que o O(l · n · c) da
multiplicação densa.

Espaço total: **O(nnz + l + c)**.

## Decisões de implementação e otimizações

- **Construção sequencial O(1) por nó** (`Construtor`): soma, multiplicação,
  transposição e `matriz_de_densa` geram o resultado em ordem de linhas.
  Guardando o último nó de cada coluna e o "elo" da linha atual, cada nó novo
  é anexado ao fim das duas listas sem percorrer nada.
- **Soma por intercalação**: cada par de linhas é fundido como no merge sort
  (listas já ordenadas); somas que se anulam **não** geram nó.
- **Multiplicação de Gustavson**: para cada linha *i* de A, cada `a[i][k]`
  é multiplicado pela linha *k* de B e acumulado num vetor denso de *c*
  posições. Um vetor de **marcas** (`marca[j] = i + 1`) evita zerar o
  acumulador a cada linha, e só as colunas tocadas são ordenadas e anexadas.
  Cancelamentos exatos não viram nó.
- **Transposição sem ordenação**: a coluna *j* de A, lida de cima para baixo,
  já é a linha *j* de Aᵀ em ordem crescente — basta anexar.
- **Ponteiro para ponteiro** em `matriz_definir`: inserir, atualizar e remover
  usam o mesmo código para início, meio e fim das listas; atualizar um valor
  existente nem toca na lista da coluna.
- **Invariantes verificáveis** (`matriz_verificar`): confere ordenação,
  coordenadas, ausência de zeros armazenados, contagem de não zeros e que as
  listas de linha e de coluna encadeiam **os mesmos nós** — usada em todos os
  testes.
- **Robustez**: tipo opaco, todos os `malloc`/`calloc` verificados (em falha,
  nada vaza e a função devolve `NULL`/`false`), índices validados,
  `matriz_multiplicar_vetor` recusa `x` e `y` sobrepostos.
- **Testes contra implementação densa** em dezenas de matrizes aleatórias de
  formatos e densidades variados: valores inteiros (comparação **exata**) e
  reais (tolerância relativa de 1e-9, pois a ordem das somas pode mudar o
  último bit), além de A + (−A), identidade, cancelamentos e uma matriz
  2000 × 2000 com Aᵀ·A.

## Estrutura de pastas

```
matriz-esparsa-c/
├── .github/workflows/ci.yml       # integração contínua (testes com sanitizers + demo)
├── include/matriz_esparsa.h       # API pública (struct opaca + documentação)
├── src/matriz_esparsa.c           # implementação
├── src/main.c                     # programa de demonstração
├── tests/test_matriz_esparsa.c    # testes comparando com a versão densa
├── Makefile
├── LICENSE
└── README.md
```

## Como compilar, executar e testar

### Linux / macOS (com `make`)

```bash
make            # compila demo e testes em build/
make run        # executa a demonstração
make test       # executa os testes
make test SANITIZE=1   # testes com AddressSanitizer + UBSan
make clean      # remove build/
```

### Windows (MinGW-w64 / GCC, sem `make`)

```bat
mkdir build
gcc -std=c11 -Wall -Wextra -Wpedantic -O2 -Iinclude src/matriz_esparsa.c src/main.c -o build/demo.exe
gcc -std=c11 -Wall -Wextra -Wpedantic -O2 -Iinclude src/matriz_esparsa.c tests/test_matriz_esparsa.c -o build/test_matriz_esparsa.exe
build\demo.exe
build\test_matriz_esparsa.exe
```

(Com o MinGW, `mingw32-make` também funciona no lugar de `make`.)

Saída esperada dos testes:

```
Testes da matriz esparsa (comparando com implementacao densa):
  [ok] teste_criar
  [ok] teste_definir_obter_remover
  [ok] teste_definir_aleatorio_contra_densa
  [ok] teste_densa_ida_e_volta
  [ok] teste_somar
  [ok] teste_somar_cancelamento_e_dimensoes
  [ok] teste_transpor
  [ok] teste_multiplicar
  [ok] teste_multiplicar_casos_especiais
  [ok] teste_multiplicar_vetor
  [ok] teste_matriz_grande_e_esparsa
  [ok] teste_ponteiros_nulos
12 testes passaram
```

## Exemplo de saída

Saída real de `make run`:

```
== Matriz esparsa com listas ortogonais ==

A (4 x 5, 6 nao zeros, densidade 30.0%):
  [    5.00    0.00    0.00   -1.50    0.00 ]
  [    0.00    0.00    2.00    0.00    0.00 ]
  [    0.00    0.00    0.00    0.00    7.25 ]
  [    0.00    3.00    0.00    4.00    0.00 ]
A em triplas (so os nao zeros):
  (0, 0) = 5
  (0, 3) = -1.5
  (1, 2) = 2
  (2, 4) = 7.25
  (3, 1) = 3
  (3, 3) = 4
obter A[2][4] = 7.25, obter A[2][3] = 0

A^T (transposta) (5 x 4, 6 nao zeros, densidade 30.0%):
  [    5.00    0.00    0.00    0.00 ]
  [    0.00    0.00    0.00    3.00 ]
  [    0.00    2.00    0.00    0.00 ]
  [   -1.50    0.00    0.00    4.00 ]
  [    0.00    0.00    7.25    0.00 ]

B (4 x 5, 5 nao zeros, densidade 25.0%):
  [    1.00    0.00    0.00    1.50    0.00 ]
  [    0.00    0.00   -2.00    0.00    0.00 ]
  [    0.00    0.00    0.00    0.00   -7.25 ]
  [    0.00    1.00    0.00    0.00    0.00 ]

A + B (4 x 5, 3 nao zeros, densidade 15.0%):
  [    6.00    0.00    0.00    0.00    0.00 ]
  [    0.00    0.00    0.00    0.00    0.00 ]
  [    0.00    0.00    0.00    0.00    0.00 ]
  [    0.00    4.00    0.00    4.00    0.00 ]
  (as posicoes (0,3), (1,2) e (2,4) se anulam: nenhum no e criado)

A * A^T (4 x 4, 6 nao zeros, densidade 37.5%):
  [   27.25    0.00    0.00   -6.00 ]
  [    0.00    4.00    0.00    0.00 ]
  [    0.00    0.00   52.56    0.00 ]
  [   -6.00    0.00    0.00   25.00 ]

A * x, com x = [1 2 3 4 5]: y = [ -1 6 36.25 22 ]

A apos definir A[0][3] = 0 e A[2][4] = 0 (4 x 5, 4 nao zeros, densidade 20.0%):
  [    5.00    0.00    0.00    0.00    0.00 ]
  [    0.00    0.00    2.00    0.00    0.00 ]
  [    0.00    0.00    0.00    0.00    0.00 ]
  [    0.00    3.00    0.00    4.00    0.00 ]

tridiagonal 100000 x 100000: 299998 nao zeros
  armazenamento denso: 76294 MiB | listas ortogonais: ~13.0 MiB
```

## Referências

- KNUTH, D. E. *The Art of Computer Programming*, v. 1: *Fundamental Algorithms*. 3. ed. Addison-Wesley, 1997. (Seção 2.2.6 — Arrays and orthogonal lists.)
- GUSTAVSON, F. G. Two fast algorithms for sparse matrices: multiplication and permuted transposition. *ACM Transactions on Mathematical Software*, v. 4, n. 3, p. 250–269, 1978.
- SAAD, Y. *Iterative Methods for Sparse Linear Systems*. 2. ed. SIAM, 2003. (Cap. 3 — Sparse matrices, formatos COO e CSR.)
- TENENBAUM, A. M.; LANGSAM, Y.; AUGENSTEIN, M. J. *Estruturas de Dados Usando C*. São Paulo: Makron Books, 1995.
- CELES, W.; CERQUEIRA, R.; RANGEL, J. L. *Introdução a Estruturas de Dados: com técnicas de programação em C*. 2. ed. Elsevier, 2016.

## Licença

Distribuído sob a licença MIT. Veja [LICENSE](LICENSE).
