#include <assert.h>
#include <iostream>

using namespace std;

int main()
{
    const int n = 100;
    double bchoosek[n+1];
    for (int k=0;k<n+1;k++) { bchoosek[k] = 0.0; } 
    // compute binomial coefficients
    bchoosek[0] = 1.0;
    for (int k=0;k<n;k++) {
	cout << "n = " << k+1 << endl;
	for (int i=k+1;i>0;i--) {
	    bchoosek[i] += bchoosek[i-1];
	}
	for (int i=0;i<=k+1;i++) {
	    cout << bchoosek[i] << " ";
	}
	cout << endl;
    }
    
}
