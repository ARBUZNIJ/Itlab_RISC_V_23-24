#include <iostream>
#include <cstdlib>
#include <chrono>
#include <omp.h>
#include "QR_header.h"

int main()
{
	//omp_set_num_threads(24);
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
			auto start{ chrono::steady_clock::now() };

			t.HHolder_A();
			auto end{ chrono::steady_clock::now() };
			chrono::duration<double> elapsed_seconds = end - start;
			cout << "Time spent: " << elapsed_seconds.count() << " sec" << endl;
			start = chrono::steady_clock::now();
			t.HHolder_Q();
			end = chrono::steady_clock::now();
			elapsed_seconds = end - start;
			//cout << "Time spent: " << elapsed_seconds.count() << " sec" << endl;
			count++;
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
	return 0;
}