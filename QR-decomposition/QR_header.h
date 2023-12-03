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
	T* A, * Q, * R, * TMP, * REF_A, * v;
	T eps = 1e-10, gamma;

public:

	QR(bool flag, size_t size, size_t block_size) //flag = 0 - ввод рандомных чисел, иначе ввод с клавиатуры
	{
		n = size;
		A = new T[(n + 1) * n]();
		Q = new T[n * n]();
		R = new T[n * n]();
		TMP = new T[n * n]();
		REF_A = new T[n * n]();
		
		this->block_size = block_size;

		if (flag)
		{

			cout << "Enter matrix:" << endl;
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
				{
					cin >> A[i * n + j];
					REF_A[i * n + j] = A[i * n + j];
				}
			}

		}
		else {


			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
					A[i * n + j] = T(rand()) / RAND_MAX;

			for (i = 0; i < n; i++)
				A[i * n + i] += n;

			//#pragma omp parallel for schedule(static) private(i)
			for (i = 0; i < n; i++)
				copy(A + i * n, A + (i + 1) * n, REF_A + i * n);
		}
	}
	void form_v_gamma(size_t ind)
	{
		v = new T[n](); //дл€ типа T должен существовать конструктор по умолчанию;
		T* u = new T[n - ind];
		T scl = 0;
		//#pragma omp parallel for simd reduction(+: scl) 
		 
		for (size_t i = 0; i < n - ind; i++)
		{
			u[i] = A[(i + ind) * n + ind];
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
			v[m] = A[(m + 1) * n + k];
		gamma = abs(v[k]);
		return;

	}
	T sgn(T val)                //неканоничный sign(x), который возвращает 1 дл€ неотрицательных, иначе -1 68719476736
	{
		if (val >= 0)
			return 1;
		else return -1;
	}

	T scal(size_t ind)  //функци€ скал€рного произведени€ дл€ ind столбцa матрицы A и вектора v
	{
		T res = 0;
		for (size_t i = 0; i < n; i++)
			res += A[i * n + ind] * v[i];
		return res;
	}
	void HHolder_A()
	{
		size_t m = 0;

		for (; m < (n / block_size); m++)
		{
			HHolder_Block(m * block_size, block_size);
		}

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
					A[i * n + k] -= v[i] * factor[k - j];
			}

			//#pragma omp parallel for simd schedule(static) private(k)
			for (k = j + 1; k <= n; k++)
				A[k * n + j] = v[k - 1];
		}
	}
	void HHolder_R()
	{
		//#pragma omp parallel for private(i)
		for (i = 0; i < n; i++)
			copy(A + i * (n + 1), A + (i + 1) * n, R + i * (n + 1));
	}
	void HHolder_Q()
	{
		//#pragma omp parallel for private(i)
		for (i = 0; i < n; i++)
			Q[i * n + i] = TMP[i * n + i] = 1.0;

		k = n - 1;
		form_v(k);

		gamma = -1 / gamma;
		Q[k * n + k] = v[k] * v[k] * gamma + 1.0;

		for (k = n - 2; k >= 0; k--)
		{
			form_v(k);
			gamma = -1 / gamma;

			for (i = k; i < n; i++)
				//#pragma omp parallel for private(j)
				for (j = k; j < n; j++)
				{
					TMP[i * n + j] = v[i] * v[j] * gamma;
					if (i == j) TMP[i * n + j]++;
				}

			Q_mult_TMP_put_Q(k);

			if (k == 0) break;
		}

	}
	void Q_mult_TMP_put_Q(size_t ind)
	{
		T* RES = new T[n * n]();

		#pragma omp parallel for private(i)
		for (size_t i = 0; i < ind; i++)
			RES[i * n + i] = 1.0;

		for (size_t i = ind; i < n; i++)
			for (size_t k = ind; k < n; k++)
				//#pragma omp parallel for simd
				for (size_t j = ind; j < n; j++)
					RES[i * n + j] += Q[i * n + k] * TMP[k * n + j];

		swap(RES, Q);

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
		delete[]A;
		delete[]TMP;
		delete[]REF_A;
		delete[]v;

		Q = R = A =	TMP = REF_A = v = nullptr;
		i = j = k = 0;
	}
};