#include<iostream>
#include<cstdlib>
#include<cmath>
#include <chrono>

using namespace std;
template <typename T>
class QR
{
private:
	size_t n, i, j, k;
	T** A; T** Q; T** R; T** TMP; T** REF_A;
	T* v;
	T eps = 1e-7, gamma;
public:
	QR(bool flag, size_t n) //flag = 0 - ввод рандомных чисел, иначе ввод с клавиатуры
	{
		this->n = n;
		this->n = n;
		A = new T * [n+1];
		REF_A = new T * [n];
		if (flag)
		{
			cout << "Enter matrix:" << endl;
			for (i = 0; i < n; i++)
			{
				A[i] = new T[n];
				REF_A[i] = new T[n];
				for (j = 0; j < n; j++)
				{
					cin >> A[i][j];
					REF_A[i][j] = A[i][j];
				}
			}
		}
		else {
			for (i = 0; i < n; i++)
			{
				A[i] = new T[n];
				REF_A[i] = new T[n];
				for (j = 0; j < n; j++)
				{
					REF_A[i][j] = A[i][j] = T(rand()) / RAND_MAX * 1000000 - 500000;
				}
			}
		}
	}
	void form_v_gamma(size_t ind)
	{
		v = new T[n](); //дл€ типа T должен существовать конструктор по умолчанию; предполагаетс€, что при его вызове у каждой компоненты вектора вектор обнулитс€
		T* u = new T[n-ind];
		T scl=0;
		for (size_t i = 0; i < n - ind; i++)
		{
			u[i] = A[i + ind][ind];						//получаем s и (s,s)
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
			for (size_t i = 0; i < n - ind; i++)
			{
				u[i] *= scl;
			}
			gamma = (1 + abs(u[0]));
			v[ind] = sgn(u[0]) * gamma;
			for (size_t i = ind + 1; i < n; i++)
			{
				v[i] = u[i - ind];
			}
			return;
		}
	}
	void form_v(size_t k)
	{
		v = new T[n]();
		for (size_t m = k; m < n; m++)
			v[m] = A[m+1][k];
		gamma = abs(v[k]);
		return;
	}
	T sgn(T val)                //неканоничный sign(x), который возвращает 1 дл€ неотрицательных, иначе -1
	{
		if (val >= 0)
			return 1;
		else return -1;
	}

	T scal(size_t ind)  //функци€ скал€рного произведени€ дл€ ind столбцa матрицы A и вектора v
	{
		T res = 0;
		for (size_t i = 0; i < n; i++)
		{
			res += A[i][ind] * v[i];
		}
		return res;
	}
	void HHolder_A(size_t r)
	{
		A[n] = new T[n];
		for (j = 0; j < r; j++)
		{
			form_v_gamma(j);
			for (k = j; k < n; k++)
			{
				T factor = scal(k) / gamma;
				for (i = 0; i < n; i++)
					A[i][k] = A[i][k] - factor * v[i];
			}
			for (k = j+1 ; k <= n; k++)
			{
				A[k][j] = v[k-1];
			}
		}			
	}
	void HHolder_R(size_t r)
	{
		R = new T * [n];
		for (i = 0; i < n; i++)
		{
			R[i] = new T[n]();
			for (j = i; j < r; j++)
				R[i][j] = A[i][j];
		}
	}
	void HHolder_Q()
	{
		Q = new T * [n];
		TMP = new T * [n];
		for (i = 0; i < n; i++)
		{
			Q[i] = new T[n]();
			Q[i][i] = 1.0;
			TMP[i] = new T[n]();
			TMP[i][i] = 1.0;
		}
		k = n - 1;
		form_v(k); //+
		gamma = -1 / gamma;
		Q[k][k] = v[k] * v[k] * gamma + 1.0;
		for (k = n - 2; k >= 0; k--)
		{
			form_v(k);
			gamma = -1 / gamma;
			for (i = k; i < n; i++)
				for (j = k; j < n; j++)
				{
					TMP[i][j] = v[i] * v[j] * gamma;
					if (i == j) TMP[i][j]++;
				}

			Q_mult_TMP_put_Q(k);

			if (k == 0) break;
		}
		return;
	}
	void Q_mult_TMP_put_Q(size_t ind)
	{
		T** RES = new T * [n];
		for (size_t i = 0; i < n; i++)
			RES[i] = new T[n]();
		for (size_t i = 0; i < ind; i++)
			RES[i][i] = 1.0;
		for (size_t i = ind; i < n; i++)
		{
			for (size_t k = ind; k < n; k++)
			{
				for (size_t j = ind; j < n; j++)
				{
					RES[i][j] += Q[i][k] * TMP[k][j];
				}
			}
		}
		swap(RES, Q);
		return;
	}
	bool check()
	{
		T** tmp = new T * [n];
		for (i = 0; i < n; i++)
		{
			tmp[i] = new T[n];
			for (j = 0; j < n; j++)
				tmp[i][j] = 0;
		}
		for (i = 0; i < n; i++)
			for (k = 0; k < n; k++)
				for (j = 0; j < n; j++)
					tmp[i][j] += Q[i][k] * R[k][j];
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				if (abs(REF_A[i][j] - tmp[i][j]) >= eps)
					return false;
		return true;
	}
	void out(char s)
	{
		if(s=='R')
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
				cout <<R[i][j] << ' ';
			cout << endl;
		}
		else if(s=='Q')
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
					cout << Q[i][j] << ' ';
				cout << endl;
			}
		else if (s == 'A')
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
					cout << A[i][j] << ' ';
				cout << endl;
			}
	}
	void transpQ()
	{
		for (i = 0; i < n; i++)
			for (j = i + 1; j < n; j++)
				swap(Q[i][j], Q[j][i]);
		return;
	}
	~QR()
	{
		delete[]Q;
		delete[]R;
		delete[]A;
		delete[]v;
		Q = R = A = nullptr;
		v = nullptr;
		i = j = k = n = 0;
	}
};

int main()
{
	size_t n, r;
	cout << "Enter matrix dimentions: " << endl;
	cin >> n;
	r = n;
	QR <double> t(0, n);
	auto start{ chrono::steady_clock::now() };
	t.HHolder_A(r);
	/*cout << endl;
	t.out('A');
	cout << endl;*/
	t.HHolder_R(r);
	/*cout << endl;
	t.out('R');
	cout << endl;*/
	t.HHolder_Q();
	t.transpQ();
	//cout << endl;
	//t.out('Q');
	//cout << endl;
	auto end{ chrono::steady_clock::now() };
	chrono::duration<double> elapsed_seconds = end - start;
	cout << "Time spent: " << elapsed_seconds.count() << " sec" << endl;

	cout << "QR-decomposition is finished, checking the result..." << endl;
	if (t.check())
		cout << "QR-decomposition is correct";
	else cout << "QR-decomposition is incorrect";

	return 0;
}