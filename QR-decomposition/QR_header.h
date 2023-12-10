#pragma once
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <omp.h>

//const size_t bs = 32;
using namespace std;
template <typename T>
class QR
{
private:

	size_t i, j, k, n, block_size;
	T* Q, * R, * REF_A, * v;
	T eps = 1e-10, gamma;

public:

	QR(bool flag, size_t size, size_t block_size) //flag = 0 - ���� ��������� �����, ����� ���� � ����������
	{
		n = size;
		Q = new T[n * n]();
		R = new T[n * n]();
		REF_A = new T[n * n]();
		
		this->block_size = block_size;

		if (flag)
		{

			cout << "Enter matrix:" << endl;
			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
					cin >> R[i * n + j];

			for (i = 0; i < n; i++)

				copy(R + i * n, R + (i + 1) * n, REF_A + i * n);

		}
		else {


			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)

					R[i * n + j] = T(rand()) / RAND_MAX;

			for (i = 0; i < n; i++)
				R[i * n + i] += n;

			//#pragma omp parallel for schedule(static) private(i)
			for (i = 0; i < n; i++)
				copy(R + i * n, R + (i + 1) * n, REF_A + i * n);
		}
	}
	void form_v_gamma(size_t ind)
	{
		v = new T[n](); //��� ���� T ������ ������������ ����������� �� ���������;
		T* u = new T[n - ind];
		T scl = 0;
		//#pragma omp parallel for simd reduction(+: scl) 
		 
		for (size_t i = 0; i < n - ind; i++)
		{
			u[i] = R[(i + ind) * n + ind];
			scl += u[i] * u[i];
		}

		if (scl < eps)
		{
			v[ind] = 1;
			gamma = 0.5;
			return;
		}

		else
		{
			scl = 1 / sqrt(scl);
			u[0] *= scl;
			gamma = (1 + abs(u[0]));
			v[ind] = sgn(u[0]) * gamma;

			//#pragma omp parallel for simd
			for (size_t i = ind + 1; i < n; i++)
				v[i] = u[i - ind] * scl;

			return;
		}
	}
	void form_v(size_t k)
	{
		v = new T[n]();

		//#pragma omp parallel for simd
		for (size_t m = k; m < n; m++)
			v[m] = R[(m + 1) * n + k];

		gamma = abs(v[k]);

		return;
	}
	T sgn(T val)                //������������ sign(x), ������� ���������� 1 ��� ���������������, ����� -1 68719476736
	{
		if (val >= 0)
			return 1;
		else return -1;
	}

	T scal(size_t ind)  //������� ���������� ������������ ��� ind ������a ������� A � ������� v
	{
		T res = 0;
		for (size_t i = 0; i < n; i++)
			res += R[i * n + ind] * v[i];
		return res;
	}
	void HHolder_A()
	{
		size_t m = 0;

		for (; m < (n / block_size); m++)
		
			HHolder_Block(m * block_size, block_size);

		if (n % block_size != 0)
			HHolder_Block(m * block_size, n - m * block_size);
	}
	void HHolder_Block(size_t i_start, size_t b_size)
	{

		T* factor;
		for (j = i_start; j < i_start + b_size; j++)
		{
			form_v_gamma(j);

			factor = new T[n - j]();

//#pragma omp parallel for private(k)
			for (k = j; k < n; k++)
				factor[k - j] = scal(k) / gamma;

			for (i = j; i < n; i++)
			{

//#pragma omp parallel for simd private(k) 
				for (k = j; k < n; k++)
					R[i * n + k] -= v[i] * factor[k - j];
			}
		}
	}
	void HHolder_Q()
	{
		T* factor;

		for (i = 0; i < n; i++)
			copy(REF_A + i * n, REF_A + (i + 1) * n, Q + i * n);

		for (i = 0; i < n; i++)
		{

			for (j = 0; j < n; j++)
				Q[j * n + i] /= R[i * n + i];

			factor = new T[n - i - 1];

			copy(R + i * (n + 1) + 1, R + (i + 1) * n, factor);

			for (k = 0; k < n; k++)
				for (j = i + 1; j < n; j++)

					Q[k * n + j] -= factor[j-(i+1)] * Q[k * n + i];

		}
	}


	bool check()
	{
		T* tmp = new T[n * n]();

		for (i = 0; i < n; i++)
			for (k = 0; k < n; k++)
				for (j = 0; j < n; j++)
					tmp[i * n + j] += Q[i * n + k] * R[k * n + j];

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				if (abs(REF_A[i * n + j] - tmp[i * n + j]) >= eps)
					return false;

		return true;
	}

	void out(char s)
	{
		cout << endl;
		if (s == 'R')

			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
					cout << R[i * n + j] << ' ';
				cout << endl;
			}

		else if (s == 'Q')

			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
					cout << Q[i * n + j] << ' ';
				cout << endl;
			}

		else if (s == 'A')

			for (i = 0; i < n + 1; i++)
			{
				for (j = 0; j < n; j++)
					cout << A[i * n + j] << ' ';
				cout << endl;
			}
		cout << endl;
	}

	void transpQ()
	{
		for (i = 0; i < n; i++)

			//#pragma omp parallel for private(j)
			for (j = i + 1; j < n; j++)
				swap(Q[i * n + j], Q[j * n + i]);

		return;
	}

	~QR()
	{
		delete[]Q;
		delete[]R;
		delete[]v;
		delete[]REF_A;

		Q = R = REF_A = v = nullptr;
		i = j = k = 0;
	}
};