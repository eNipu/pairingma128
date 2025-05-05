CC = gcc
CFLAGS = -Wall -Wextra -O2
LDFLAGS = -lgmp

SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

BLS12_SRC = $(wildcard BLS12/*.c)
BLS12_OBJ = $(BLS12_SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
BLS12_BIN = $(BIN_DIR)/bls12

BN_SRC = $(wildcard BN/*.c)
BN_OBJ = $(BN_SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
BN_BIN = $(BIN_DIR)/bn

KSS16_SRC = $(wildcard KSS16/*.c)
KSS16_OBJ = $(KSS16_SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
KSS16_BIN = $(BIN_DIR)/kss16

all: $(BLS12_BIN) $(BN_BIN) $(KSS16_BIN)

$(BLS12_BIN): $(BLS12_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BN_BIN): $(BN_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(KSS16_BIN): $(KSS16_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean
