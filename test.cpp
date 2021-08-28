#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <chrono>

#include "matrixvector.hpp"
#include "calculus.hpp"


using namespace std;

const string RESET = "\033[0m", RED = "\033[0;31m", GREEN = "\033[0;32m", RED_B = "\033[1;31m", GREEN_B = "\033[1;32m";

double time_f(function<bool(void)> f){
    auto start = chrono::system_clock::now();
    for(int i = 0; i < 100000; i++) f();
    auto stop = chrono::system_clock::now();

    return chrono::duration_cast<chrono::nanoseconds>(stop-start).count() / 100000.0;
}

bool test_operation(bool github, vector<function<bool(void)>> check_functions, string test_name, bool time = true){
    int index = 1; bool success = 1, tmp; double total_time = 0;
    if(github) cout << "::group::";
    cout << RESET << "Running " << check_functions.size() << " subtests (" << test_name << ")" << endl;

    for(auto f : check_functions) {
        cout << RESET << "Test case " << index++ << ": " << ((tmp = f()) ? (GREEN + "passed.") : (RED + "failed.")) << endl;
        if(time) { total_time += time_f(f); } success &= tmp;
    }
    cout << RESET <<(success ? (GREEN_B + "[ Test passed ]") : (RED_B + "[ Test failed ]")) << RESET;
    if(time) cout << " in "  << total_time << " ns" << endl;
    
    if(github) cout << "::endgroup::" << endl;
    else cout << endl;

    return success;
}

int main(int argc, char* argv[]) {
    bool github = stoi(argv[1]);
    cout << endl;

    bool s = 1;
    
    s &= test_operation(github, {
        []()->bool{ return matrixvector::vector2(3, 5) + matrixvector::vector2(6, -5) == matrixvector::vector2(9, 0); },
        []()->bool{ return matrixvector::vector3(3, 5, -1) + matrixvector::vector3(6, -5, 2) == matrixvector::vector3(9, 0, 1); },
        []()->bool{ return matrixvector::vector2(3, 5) - matrixvector::vector2(6, -5) == matrixvector::vector2(-3, 10); },
        []()->bool{ return matrixvector::vector3(3, 5, -1) - matrixvector::vector3(6, -5, 2) == matrixvector::vector3(-3, 10, -3); }
    }, "vector addition and subtraction");

    s &= test_operation(github, {
        []()->bool{ return (matrixvector::vector2(3, 5) * 2) == matrixvector::vector2(3*2, 5*2); },
        []()->bool{ return (matrixvector::vector3(4, 2, 3) / 3) == matrixvector::vector3(4/3.0, 2/3.0, 3/3.0); },
        []()->bool{ return (matrixvector::vector_f<4>({4, 2, 3, 5}) * -1.5) == matrixvector::vector_f<4>({-4*1.5, -2*1.5, -3*1.5, -5*1.5}); }
    }, "vector scaling");

    s &= test_operation(github, {
        []()->bool{ return (matrixvector::vector2(3, 5) * matrixvector::vector2(6, -5)) == (3*6 - 5*5); },
        []()->bool{ return (matrixvector::vector3(4, 2, 3) * matrixvector::vector3(9, 4, 1)) == (4*9 + 2*4 + 3*1); },
        []()->bool{ return (matrixvector::vector_f<4>({4, 2, 3, 5}) * matrixvector::vector_f<4>({9, 4, 1, -3})) == (4*9 + 2*4 + 3*1 - 5*3); }
    }, "vector dot product");

    s &= test_operation(github, {
        []()->bool{ return matrixvector::vector2(3, 5).size() == 2; },
        []()->bool{ return matrixvector::vector3(4, 2, 3).size() == 3; },
        []()->bool{ return matrixvector::vector3(4, 2, 3)[0] = 4 && matrixvector::vector3(4, 2, 3)[1] == 2 && matrixvector::vector3(4, 2, 3)[2] == 3; },
        []()->bool{ return abs(matrixvector::abs(matrixvector::vector3(4, 2, 3)) - 5.38516481) < 1e-4; }
    }, "vector properties");

    s &= test_operation(github, {
        []()->bool{ return matrixvector::apply_function(matrixvector::vector2(3, 5), [](double d) { return d-1; }) == matrixvector::vector2(2, 4); },
        []()->bool{ return matrixvector::apply_function(matrixvector::vector3(3, 5, 5), [](double d) { return d-1; }) == matrixvector::vector3(2, 4, 4); },
        []()->bool{ return matrixvector::apply_function(matrixvector::vector2(3, 5), matrixvector::vector2(1, 2), [](double d1, double d2) { return d1*d2; }) == matrixvector::vector2(3, 10); },
        []()->bool{ return matrixvector::apply_function(matrixvector::vector3(3, 5, 5), matrixvector::vector3(1, 2, 3), [](double d1, double d2) { return d1*d2; }) == matrixvector::vector3(3, 10, 15); }
    }, "vector functions");

    cout << endl << (s ? GREEN_B + "[ All tests passed ]" : RED_B + "[ Tests failed ]") << RESET << endl;

    return s == 0;
}