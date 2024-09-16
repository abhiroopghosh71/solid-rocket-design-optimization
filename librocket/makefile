CC = gcc

default: libapprocket.a

libapprocket.a: approcket.o
	ar rcs $@ $^
    
approcket.o: approcket.c
	$(CC) -g -fPIC -c $<
	$(CC) -g test_rocket.c -o test_rocket.out -lm

clean:
	rm *.o *.a
