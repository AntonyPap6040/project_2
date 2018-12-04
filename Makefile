OBJS1	= cube.o
OBJS2	= cube.cpp
OBJS3	= item.cpp
OBJS4	= item.h
OBJS5	= euclidean.cpp
OBJS6	= euclidean.h
OBJS7	= cosine.cpp
OBJS8	= cosine.h
OBJS9	= hyper.cpp
OBJS10	= hyper.h
OBJS11	= item.o
OBJS12	= euclidean.o
OBJS13	= cosine.o
OBJS14	= hyper.o
OBJS21	= initialization.h
OBJS20	= initialization.cpp
OBJS22	= initialization.o
OBJS18	= assignment.cpp
OBJS19	= assignment.h
OBJS23	= assignment.o
OBJS15	= argClass.cpp
OBJS16	= argClass.h
OBJS17	= argClass.o
OBJS25 = cluster.cpp
out	= cluster

CC	= g++
FLAGS1	= -c -w -Wall -g3
FLAGS2	= -w -o


# main
all:
	$(CC) $(FLAGS1) $(OBJS3) $(OBJS4) 
	$(CC) $(FLAGS1) $(OBJS5) $(OBJS6) 
	$(CC) $(FLAGS1) $(OBJS7) $(OBJS8) 
	$(CC) $(FLAGS1) $(OBJS9) $(OBJS10)
	$(CC) $(FLAGS1) $(OBJS15) $(OBJS16) 
	$(CC) $(FLAGS1) $(OBJS20) $(OBJS21)  
	$(CC) $(FLAGS1) $(OBJS18) $(OBJS19)
	$(CC) $(FLAGS1) $(OBJS2) $(OBJS23)  
	$(CC) $(FLAGS2) $(out) $(OBJS25) $(OBJS11) $(OBJS12) $(OBJS13) $(OBJS14) $(OBJS17) $(OBJS22) $(OBJS23) $(OBJS1)
	
# clean up
clean:
	rm  $(OUT1)
	rm  $(OUT2)
	rm  $(OBJS11)
	rm  $(OBJS12)
	rm  $(OBJS13)
	rm  $(OBJS14)
	
