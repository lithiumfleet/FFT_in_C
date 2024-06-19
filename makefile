# Define the compiler
CC = gcc

# Define compiler flags
CFLAGS = -Wall -g -O

.PHONY: run
run:
	@echo "start compile"
%.exe : %.c
	fn=$(basename $(notdir $@))
	$(CC) $(CFLAGS) ./src/$(basename $(notdir $@)).c -o ./bin/$(basename $(notdir $@)).exe
	./bin/$(basename $(notdir $@)).exe

# Rule to clean the generated files
.PHONY: clean
clean:
	rm -f ./bin/*.exe
