# MathLibCpp
A C++ library that adds matrices, vectors and more.

The library is seperated into 3 parts: `matrixvector`, `calculus` and `complex`.

## Complex:

Example usage:
```c++
import "complex.hpp"
using complex::I;
...

int main(){
  complex::c_double a = 2 + 3*I;
  complex::c_double b = 3 + 4*I;
  
  cout << a * b << endl;                       // -6 + 17i
  cout << a / b << endl;                       // 0.72 + 0.04i
  cout << complex::root_p(3 + 2*I, 4) << endl; // 1.36312 + 0.201835i
  cout << complex::sqrt(-10) << endl;          // 0 + 3.16228i
}

```
