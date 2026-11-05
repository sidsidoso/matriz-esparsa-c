# Makefile - matriz-esparsa-c
# Uso: make | make run | make test | make test SANITIZE=1 | make clean

CC      ?= cc
CFLAGS  = -std=c11 -Wall -Wextra -Wpedantic -O2
CPPFLAGS += -Iinclude
LDFLAGS ?=
LDLIBS  = -lm

# Com SANITIZE=1 os binarios sao instrumentados com AddressSanitizer e
# UndefinedBehaviorSanitizer (vazamentos, acessos invalidos e UB abortam).
ifeq ($(SANITIZE),1)
CFLAGS  += -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer -g
LDFLAGS += -fsanitize=address,undefined
endif

# No Windows (MinGW) o gcc gera .exe; o sufixo evita recompilar a toa.
ifeq ($(OS),Windows_NT)
EXE     := .exe
endif

BUILD   := build
LIB_SRC := src/matriz_esparsa.c
HEADER  := include/matriz_esparsa.h
DEMO    := $(BUILD)/demo$(EXE)
TESTES  := $(BUILD)/test_matriz_esparsa$(EXE)

.PHONY: all run test clean

all: $(DEMO) $(TESTES)

$(BUILD):
	mkdir -p $@

$(DEMO): $(LIB_SRC) src/main.c $(HEADER) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LIB_SRC) src/main.c -o $@ $(LDFLAGS) $(LDLIBS)

$(TESTES): $(LIB_SRC) tests/test_matriz_esparsa.c $(HEADER) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LIB_SRC) tests/test_matriz_esparsa.c -o $@ $(LDFLAGS) $(LDLIBS)

run: $(DEMO)
	./$(DEMO)

test: $(TESTES)
	./$(TESTES)

clean:
	rm -rf $(BUILD)
