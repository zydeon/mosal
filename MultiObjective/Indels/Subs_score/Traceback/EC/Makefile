CC   = gcc -Wall -O3
OBJS = main.o
PROG = prog

# GENERIC

all:	 ${PROG}

clean:
		 rm ${OBJS} ${PROG}

${PROG}: ${OBJS}
		 ${CC} ${OBJS} -o $@
.c.o:
		 ${CC} $< -c -o $@

#############################

main.o: main.c
