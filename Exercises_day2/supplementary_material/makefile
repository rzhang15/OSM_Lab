#########################################
all: normalize_vec.exec dot_prod.exec serial_pi.exec
########################################
normalize_vec.exec : normalize_vec.cpp
	g++ normalize_vec.cpp -fopenmp -o normalize_vec.exec

dot_prod.exec : dot_prod.cpp
	g++ dot_prod.cpp -fopenmp -o dot_prod.exec

serial_pi.exec : serial_pi.cpp
	g++ serial_pi.cpp -fopenmp -o serial_pi.exec
########################################
clean :
	rm -f *.exec
