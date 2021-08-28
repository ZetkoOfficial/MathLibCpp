test: compile_test
	./test

compile_test: 
	g++ -O2 test.cpp -o test