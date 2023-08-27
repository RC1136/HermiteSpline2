
SOURCES=dllmain.c functions.c hermite.c
NAME=HermiteLib
DEBUG_FLAGS=-shared -static -Wall -g -O0
RELEASE_FLAGS=-shared -static -w -O2
DEBUG_DIR=../Debug
RELEASE_DIR=../Release

all: gcc32 move

gcc32: gcc32_q

gcc32_d:
	i686-w64-mingw32-gcc -m32 $(DEBUG_FLAGS) $(SOURCES) -o $(NAME)_$@_D.dll
	i686-w64-mingw32-gcc -m32 $(RELEASE_FLAGS) $(SOURCES) -o $(NAME)_$@_R.dll
	cp $(NAME)_$@_D.dll $(NAME)_D.dll
	cp $(NAME)_$@_R.dll $(NAME)_R.dll

gcc32_ld:
	i686-w64-mingw32-gcc -m32 $(DEBUG_FLAGS) -D__LONG_DOUBLE $(SOURCES) -o $(NAME)_$@_D.dll
	i686-w64-mingw32-gcc -m32 $(RELEASE_FLAGS) -D__LONG_DOUBLE $(SOURCES) -o $(NAME)_$@_R.dll
	cp $(NAME)_$@_D.dll $(NAME)_D.dll
	cp $(NAME)_$@_R.dll $(NAME)_R.dll

gcc32_q:
	i686-w64-mingw32-gcc -m32 $(DEBUG_FLAGS) -D__FLOAT_128 $(SOURCES) -o $(NAME)_$@_D.dll -lquadmath
	i686-w64-mingw32-gcc -m32 $(RELEASE_FLAGS) -D__FLOAT_128 $(SOURCES) -o $(NAME)_$@_R.dll -lquadmath
	cp $(NAME)_$@_D.dll $(NAME)_D.dll
	cp $(NAME)_$@_R.dll $(NAME)_R.dll

move:
	mv $(NAME)_R.dll $(RELEASE_DIR)/$(NAME).dll
	mv $(NAME)_D.dll $(DEBUG_DIR)/$(NAME).dll

clean:
	rm -f *.dll

