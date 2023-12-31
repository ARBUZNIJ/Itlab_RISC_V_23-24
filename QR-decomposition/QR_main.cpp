#include <iostream>
#include <cstdlib>
#include <omp.h>
#include "QR_header.h"

int main()
{
	int count;
	size_t n, block_size;
	cout << "Enter matrix dimentions: " << endl;
	cin >> n;
	while (1)
	{
		cout << "Enter block size: ";
		cin >> block_size;
		count = 0;
		while (count < 5)
		{

			QR <double> t(0, n, block_size);
			auto start_A{ chrono::steady_clock::now() };
			t.HHolder_A();
			auto end_A{ chrono::steady_clock::now() };
			chrono::duration<double> elapsed_seconds_A = end_A - start_A;
			auto start_Q = chrono::steady_clock::now();
			t.HHolder_Q();
			auto end_Q = chrono::steady_clock::now();
			chrono::duration<double>  elapsed_seconds_Q = end_Q - start_Q;
			cout << "Time spent: " << elapsed_seconds_A.count() << " sec + " << elapsed_seconds_Q.count() << " sec = ";
			elapsed_seconds_A += elapsed_seconds_Q;
			cout << elapsed_seconds_A.count() << endl;
			count++;
			if (t.check())
				cout << "QR-decomposition is correct";
			else cout << "QR-decomposition is incorrect";
			cout << endl;
		}
	}
	/*t.HHolder_Q();
	end = chrono::steady_clock::now();
	elapsed_seconds = end - start;
	cout << "Time spent: " << elapsed_seconds.count() << " sec" << endl;
	if (t.check())
		cout << "QR-decomposition is correct";
	else cout << "QR-decomposition is incorrect";
	cout << endl;*/
	/*QR <double> t(0, n, block_size);
	t.out('R');*/

	return 0;
}