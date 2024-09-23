
make:
	clang src/plot.c -Iinclude -o obj/plot.o -c
	clang src/matrix.c -Iinclude -o obj/matrix.o -c
