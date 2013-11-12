#include <assert.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

class A {
public:
    A();
public:
    int a;
    float f;
    char str[35];
    long *lp;
};

A::A() : a(), f(), str(), lp() {}

/* A::A()
{
    memset(this, 0, sizeof(*this));
}*/


int main()
{
    A test = A();
    char* teststr = &(test.str[0]);
    int i = test.str[0];
    cerr << " i = " << i << endl;
    cerr << "a = " << test.a << endl;
    cerr << "f = " << test.f << endl;
    cerr << "lp = " << (int) test.lp << endl;
    exit(0);
}
