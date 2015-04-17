
include makefile.in

DIRS	=	./ 
	 
OBJ	=	./main.o ./riemi.o ./spzce.o		

all: $(PROGRAM_NAME)

$(PROGRAM_NAME) :$(OBJ) 
	$(LINKER) $(L_FLAGS) $(OBJ) -o $@

%.o: %.c
	$(CC) $(C_FLAGS) -c $< -o $@

clean: 
	@for dir in $(DIRS) ;\
		do \
		echo cleaning $$dir ;\
		(cd $$dir;  rm -rf *.a *.o *.mod *.ex* *core* out* ./extras/core* \#*  *~) ;\
		done
