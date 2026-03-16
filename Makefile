export TMP := C:\Users\Casper\AppData\Local\Temp
export TEMP := C:\Users\Casper\AppData\Local\Temp

CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -Iinclude
LDFLAGS = -lm

SRC_DIR = src
INC_DIR = include
TEST_DIR = test
BUILD_DIR = build

SRCS = $(SRC_DIR)/centroid.c $(SRC_DIR)/catalog.c $(SRC_DIR)/star_tracker.c $(SRC_DIR)/attitude.c
OBJS = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRCS))

TARGET = $(BUILD_DIR)/star_tracker

.PHONY: all clean

all: $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/main.o: $(TEST_DIR)/main.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJS) $(BUILD_DIR)/main.o
	$(CC) $(OBJS) $(BUILD_DIR)/main.o -o $@ $(LDFLAGS)

clean:
	rm -rf $(BUILD_DIR)
