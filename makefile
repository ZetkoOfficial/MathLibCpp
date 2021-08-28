test: compile_test
	./test 0
test_github: compile_test
	./test 1
compile_test: 
	g++ -O2 test.cpp -o test