
SOURCES=dllmain.c functions.c hermite.c
NAME=HermiteLib
FLAGS=-shared -static -g -Wall
DEBUG_DIR=../Debug
RELEASE_DIR=../Release

all: gcc32 move

gcc32: gcc32_q

gcc32_d:
	i686-w64-mingw32-gcc -m32 $(FLAGS) $(SOURCES) -o $(NAME)_$@.dll
	cp $(NAME)_$@.dll $(NAME).dll	

gcc32_ld:
	i686-w64-mingw32-gcc -m32 $(FLAGS) -D__LONG_DOUBLE $(SOURCES) -o $(NAME)_$@.dll
	cp $(NAME)_$@.dll $(NAME).dll

gcc32_q:
	i686-w64-mingw32-gcc -m32 $(FLAGS) -D__FLOAT_128 $(SOURCES) -o $(NAME)_$@.dll -lquadmath
	cp $(NAME)_$@.dll $(NAME).dll

move:
	mv $(NAME).dll $(DEBUG_DIR)

clean:
	rm -f *.dll
