
make:
	clang src/plot.c -Iinclude -o obj/plot.o -c -Wall -Werror
	clang src/matrix.c -Iinclude -o obj/matrix.o -c -Wall -Werror
	clang src/nmds.c -Iinclude -o obj/nmds.o -c -Wall -Werror
