
CC = clang
CFLAGS += -Wall -Werror -g

CITYCOUNT ?= 128
THREADCOUNT ?= 1
  
#0 -> Synchronous, 1 -> Asynchronous
MODE ?= 0

all:
	$(CC) $(CFLAGS) -DN=$(CITYCOUNT) -std=c99 -o main main.c -lm -lpthread
	$(CC) $(CFLAGS) -o datagen datagen.c
	$(CC) $(CFLAGS) -o serialtester serialtester.c

test: all
	./datagen $(CITYCOUNT)
	./main $(THREADCOUNT) $(MODE)
	./serialtester

maintest: all
	./datagen $(CITYCOUNT)
	./main 1 0
	./main 2 0
	./main 4 0
	./main 8 0
	./main 1 1
	./main 2 1
	./main 4 1
	./main 8 1
	./serialtester

fixtest: all
	./datagen $(CITYCOUNT)
	./main 1 1
	./serialtester
#	./main 2 1
#	./serialtester
#	./main 4 1
#	./serialtester

clean:
	rm -rf main
	rm -rf datagen
	rm -rf serialtest