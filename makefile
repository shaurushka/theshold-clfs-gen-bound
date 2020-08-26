OBJS = main.o gen_bounds_calculation.o base_functions.o  
CC = g++ -std=gnu++11
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

gen_bounds_calc : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o gen_bounds_calc

main.o : main.cpp gen_bounds_calculation.h base_functions.h 
	$(CC) $(CFLAGS) main.cpp

gen_bounds_calculation.o : gen_bounds_calculation.cpp gen_bounds_calculation.h base_functions.h 
	$(CC) $(CFLAGS) gen_bounds_calculation.cpp 
	
base_functions.o : base_functions.cpp base_functions.h
	$(CC) $(CFLAGS) base_functions.cpp

clean:
	\rm *.o gen_bounds_calc



