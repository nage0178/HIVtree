CC = gcc # cc, gcc, cl

CFLAGS = -g -w
#CFLAGS = -fast

LIBS = -lm # -lM


HIVtree : tools.c treesub.c mcmctree.c
	$(CC) $(CFLAGS) -o $@ tools.c mcmctree.c $(LIBS) -O3

clean:
	rm -f HIVtree
