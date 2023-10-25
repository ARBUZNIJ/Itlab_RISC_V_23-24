#include<iostream>
#include<cstdlib>
#include<cmath>
#include <chrono>


using namespace std;
template <typename T>
class QR
{
private:
	size_t clmn, line, i, j, k;
	T** A; T** Q;T **R;
	T eps = 1e-7;
public:
	QR(bool flag, size_t line, size_t clmn) //flag = 0 - ввод рандомных чисел, иначе ввод с клавиатуры
	{
		this->line = line;
		this->clmn = clmn;
		A = new T * [line];
		if (flag)
		{
			cout << "Enter matrix:" << endl;
			for (i = 0; i < line; i++)
			{
				A[i] = new T[clmn];
				for (j = 0; j < clmn; j++)
				{
					cin >> A[i][j];
				}
			}
		}
		else {
			for (i = 0; i < line; i++)
			{
				A[i] = new T[clmn];
				for (j = 0; j < clmn; j++)
				{
					A[i][j] = T(rand()) / RAND_MAX * 1000000 - 500000;
				}
			}
		}
	}
	T scl_clmn(T** data_1, size_t ind_1_clmn, T** data_2, size_t ind_2_clmn)
	{
		T res=0;
		for (size_t i = 0; i < line; i++)
		{
			res += data_1[i][ind_1_clmn] * data_2[i][ind_2_clmn];
		}
		return res;
	}
	T scl_sqr_clmn(T** data, size_t ind)
	{
		return scl_clmn(data, ind, data, ind);
	}
	void ort()
	{
		Q = new T * [line];						//def Q=A
		T tmp;
		for (i = 0; i < line; i++)
		{
			Q[i] = new T[clmn];
			//memcpy(Q[i], A[i], clmn);
			for (j = 0; j < clmn; j++)	
				Q[i][j] = A[i][j];

		}
		for (j = 0; j < clmn; j++)
		{
			for (k = 0; k < j; k++)
			{
				tmp = scl_clmn(A, j, Q, k) / scl_sqr_clmn(Q, k);
				for (i = 0; i < line; i++)
					Q[i][j] -= tmp * Q[i][k];
			}
		}
	}
	void norm()
	{
		T tmp;
		for (j = 0; j < clmn; j++)
		{
			tmp = sqrt(scl_sqr_clmn(Q, j));       //каким образом можно распознать нули, не вводя формальной проверки?
			if (tmp <= eps)
				tmp = 0;
			else tmp = 1 / tmp;
			for (i = 0; i < line; i++)
				Q[i][j] *= tmp;
		}
	}
	void form_R()
	{
		R = new T * [line];						//def Q=A
		for (i = 0; i < line; i++)
		{
			R[i] = new T[clmn];
			for (j = 0; j < clmn; j++)
				R[i][j] = 0;
			//memset(Q[i], 0, clmn);            //почему неправильно?
		}
		for (i = 0; i < line; i++)
			for (k = 0; k < line; k++)
				for (j = 0; j < clmn; j++)
					R[i][j] += Q[k][i] * A[k][j];
	}
	bool check()
	{
		T** tmp = new T * [line];
		for (i = 0; i < line; i++)
		{
			tmp[i] = new T[clmn];
			for (j = 0; j < line; j++)
				tmp[i][j] = 0;
		}
		for (i = 0; i < line; i++)
			for (k = 0; k < line; k++)
				for (j = 0; j < clmn; j++)
					tmp[i][j] += Q[i][k] * R[k][j];
		for (i = 0; i < line; i++)
			for (j = 0; j < line; j++)
				if (abs(A[i][j] - tmp[i][j]) >= eps)
					return false;
		return true;
	}
	void out(char s)
	{
		if(s=='R')
		for (i = 0; i < line; i++)
		{
			for (j = 0; j < clmn; j++)
				cout <<R[i][j] << ' ';
			cout << endl;
		}
		else if(s=='Q')
			for (i = 0; i < line; i++)
			{
				for (j = 0; j < clmn; j++)
					cout << Q[i][j] << ' ';
				cout << endl;
			}
		else if (s == 'A')
			for (i = 0; i < line; i++)
			{
				for (j = 0; j < clmn; j++)
					cout << A[i][j] << ' ';
				cout << endl;
			}
	}
};

int main()
{
	size_t line, clmn;
	cout << "Enter matrix dimentions: " << endl;
	cin >> line >> clmn;
	QR <double> t(0, line, clmn);
	auto start{ chrono::steady_clock::now() };
	t.ort();
	t.norm();
	t.form_R();
	auto end{ chrono::steady_clock::now() };
	//t.out('A');
	cout <<"QR-decomposition is finished, checking the result..." << endl;
	//t.out('R');
	if (t.check())
		cout << "QR-decomposition is correct";
	else cout << "QR-decomposition is incorrect";
	chrono::duration<double> elapsed_seconds = end - start;
	cout << endl << "Time spent: " << elapsed_seconds.count() << " sec" << endl;
	return 0;
}