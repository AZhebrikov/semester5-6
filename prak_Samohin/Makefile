run: prog
				prog
prog: prog.o Makefile
				g++ -W -Wall -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wcast-align -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wno-long-long -Wcast-align -Winline -Werror -pedantic -pedantic-errors -Wunused -Wuninitialized --param inline-unit-growth=1000000 --param max-inline-insns-single=10000000 --param large-function-growth=10000000 -fPIC -O2 -std=c++11 -fno-elide-constructors -o prog prog.o
prog.o: prog.cpp Makefile
				g++ -W -Wall -Wfloat-equal -Wpointer-arith -Wwrite-strings -Wcast-align -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wno-long-long -Wcast-align -Winline -Werror -pedantic -pedantic-errors -Wunused -Wuninitialized --param inline-unit-growth=1000000 --param max-inline-insns-single=10000000 --param large-function-growth=10000000 -fPIC  -O2 -c prog.cpp
