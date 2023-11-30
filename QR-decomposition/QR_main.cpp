#include <iostream>
#include <cstdlib>
#include <chrono>
#include <omp.h>
#include "QR_header.h"

int main()
{
	omp_set_num_threads(24);
	int count = 0;
	size_t n;
	cout << "Enter matrix dimentions: " << endl;
	cin >> n;
	while (count < 5)
	{
		QR <double> t(0, n);
		auto start{ chrono::steady_clock::now() };
		t.HHolder_A();
		t.HHolder_R();
		auto end{ chrono::steady_clock::now() };
		chrono::duration<double> elapsed_seconds = end - start;
		cout << "Time spent: " << elapsed_seconds.count() << " sec" << endl;
		count++;
	}
	/*start = chrono::steady_clock::now();
	t.HHolder_Q();
	t.transpQ();
	end = chrono::steady_clock::now();
	elapsed_seconds = end - start;
	cout << "Time spent: " << elapsed_seconds.count() << " sec" << endl;
	if (t.check())
		cout << "QR-decomposition is correct";
	else cout << "QR-decomposition is incorrect";
	*/
	return 0;
}