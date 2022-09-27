CC = gcc # cc, gcc, cl

CFLAGS = -g
#CFLAGS = -fast

LIBS = -lm # -lM


HIVLateTree : tools.c treesub.c mcmctree.c
	$(CC) $(CFLAGS) -o $@ tools.c mcmctree.c $(LIBS) -O3

clean:
	rm -f HIVLateTree
