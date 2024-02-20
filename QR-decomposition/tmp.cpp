////РЕЗЕРВНАЯ КОПИЯ ПРОГРАММЫ + ЛИШНИЕ ФУНКЦИИ, ПО МЕРЕ МОЕГО ПРОСВЕЩЕНИЯ ОКАЗАВШИЕСЯ НЕНУЖНЫМИ
//#include<iostream>
//#include<cstdlib>
//#include<cmath>
//#include <chrono>
//
//using namespace std;
//template <typename T>
//class QR
//{
//private:
//	size_t n, i, j, k;
//	T** A; T** Q; T** R; T** TMP; T** REF_A;
//	T* v;
//	T eps = 1e-7, gamma;
//public:
//	QR(bool flag, size_t n) //flag = 0 - ввод рандомных чисел, иначе ввод с клавиатуры
//	{
//		this->n = n;
//		this->n = n;
//		A = new T * [n + 1];
//		REF_A = new T * [n];
//		if (flag)
//		{
//			cout << "Enter matrix:" << endl;
//			for (i = 0; i < n; i++)
//			{
//				A[i] = new T[n];
//				REF_A[i] = new T[n];
//				for (j = 0; j < n; j++)
//				{
//					cin >> A[i][j];
//					REF_A[i][j] = A[i][j];
//				}
//			}
//		}
//		else {
//			for (i = 0; i < n; i++)
//			{
//				A[i] = new T[n];
//				REF_A[i] = new T[n];
//				for (j = 0; j < n; j++)
//				{
//					REF_A[i][j] = A[i][j] = T(rand()) / RAND_MAX * 1000000 - 500000;
//				}
//			}
//		}
//	}
//	void form_v_gamma(size_t ind)
//	{
//		v = new T[n](); //для типа T должен существовать конструктор по умолчанию; предполагается, что при его вызове у каждой компоненты вектора вектор обнулится
//		T* u = new T[n - ind];
//		T scl = 0;
//		for (size_t i = 0; i < n - ind; i++)
//		{
//			u[i] = A[i + ind][ind];						//получаем s и (s,s)
//			scl += u[i] * u[i];
//		}
//		if (scl < eps)
//		{
//			v[ind] = 1;
//			gamma = 0.5;
//			return;
//		}
//		else
//		{
//			scl = 1 / sqrt(scl);
//			for (size_t i = 0; i < n - ind; i++)
//			{
//				u[i] *= scl;
//			}
//			gamma = (1 + abs(u[0]));
//			v[ind] = sgn(u[0]) * gamma;
//			for (size_t i = ind + 1; i < n; i++)
//			{
//				v[i] = u[i - ind];
//			}
//			return;
//		}
//	}
//	void form_v(size_t k)
//	{
//		v = new T[n]();
//		for (size_t m = k; m < n; m++)
//			v[m] = A[m + 1][k];
//		gamma = abs(v[k]);
//		return;
//	}
//	T sgn(T val)                //неканоничный sign(x), который возвращает 1 для неотрицательных, иначе -1
//	{
//		if (val >= 0)
//			return 1;
//		else return -1;
//	}
//	T scl_n(T** data_1, size_t ind_1_n, T** data_2, size_t ind_2_n)  //функция скалярного произведения для столбцов двух матриц одинакового размера
//	{
//		T res = 0;
//		for (size_t i = 0; i < n; i++)
//		{
//			res += data_1[i][ind_1_n] * data_2[i][ind_2_n];
//		}
//		return res;
//	}
//	T scl_sqr_n(T** data, size_t ind)
//	{
//		return scl_n(data, ind, data, ind);
//	}
//	T scal(size_t ind)  //функция скалярного произведения для ind столбцa матрицы A и вектора v
//	{
//		T res = 0;
//		for (size_t i = 0; i < n; i++)
//		{
//			res += A[i][ind] * v[i];
//		}
//		return res;
//	}
//	void ort()
//	{
//		Q = new T * [n];
//		T tmp;
//		for (i = 0; i < n; i++)
//		{
//			Q[i] = new T[n];
//			//memcpy(Q[i], A[i], n);
//			for (j = 0; j < n; j++)
//				Q[i][j] = A[i][j];
//
//		}
//		for (j = 0; j < n; j++)
//		{
//			for (k = 0; k < j; k++)
//			{
//				tmp = scl_n(A, j, Q, k) / scl_sqr_n(Q, k);
//				for (i = 0; i < n; i++)
//					Q[i][j] -= tmp * Q[i][k];
//			}
//		}
//	}
//	void norm()
//	{
//		T tmp;
//		for (size_t j = 0; j < n; j++)
//		{
//			tmp = sqrt(scl_sqr_n(Q, j));       //каким образом можно распознать нули, не вводя формальной проверки?
//			for (i = 0; i < n; i++)
//				Q[i][j] *= tmp;
//		}
//	}
//	void form_R()
//	{
//		R = new T * [n];						//def Q=A
//		for (i = 0; i < n; i++)
//		{
//			R[i] = new T[n];
//			for (j = 0; j < n; j++)
//				memset(Q[i], 0, n * sizeof(T));            //почему неправильно?
//		}
//		for (i = 0; i < n; i++)
//			for (k = 0; k < n; k++)
//				for (j = 0; j < n; j++)
//					R[i][j] += Q[k][i] * A[k][j];
//	}
//	void HHolder_R()
//	{
//		A[n] = new T[n];
//		for (j = 0; j < n; j++)
//		{
//			form_v_gamma(j);
//			for (k = j; k < n; k++)
//			{
//				T factor = scal(k) / gamma;
//				for (i = 0; i < n; i++)
//					A[i][k] = A[i][k] - factor * v[i];
//			}
//			for (k = j + 1; k <= n; k++)
//			{
//				A[k][j] = v[k - 1];
//			}
//		}
//		R = new T * [n];						//def Q=A
//		for (i = 0; i < n; i++)
//		{
//			R[i] = new T[n]();
//			for (j = i; j < n; j++)
//				R[i][j] = A[i][j];
//		}
//
//	}
//	void HHolder_Q()
//	{
//		Q = new T * [n];
//		TMP = new T * [n];
//		for (i = 0; i < n; i++)
//		{
//			Q[i] = new T[n];
//			TMP[i] = new T[n];
//		}
//		k = n - 1;
//		form_v(k);
//		gamma = -1 / gamma;
//		for (i = 0; i < n; i++)
//		{
//			for (j = 0; j < n; j++)
//			{
//				Q[i][j] = v[i] * v[j] * gamma;
//				if (i == j) Q[i][j] += 1;
//			}
//		}
//		for (k = n - 2; k >= 0; k--)
//		{
//			form_v(k);
//			gamma = -1 / gamma;
//			for (i = 0; i < n; i++)
//			{
//				for (j = 0; j < n; j++)
//				{
//					TMP[i][j] = v[i] * v[j] * gamma;
//					if (i == j) TMP[i][j] += 1;
//				}
//			}
//			Q_mult_TMP_put_Q();
//			if (k == 0) break;
//		}
//		return;
//	}
//	void Q_mult_TMP_put_Q()
//	{
//		T** RES = new T * [n];
//		for (size_t i = 0; i < n; i++)
//		{
//			RES[i] = new T[n]();
//		}
//		for (size_t i = 0; i < n; i++)
//		{
//			for (size_t k = 0; k < n; k++)
//			{
//				for (size_t j = 0; j < n; j++)
//				{
//					RES[i][j] += Q[i][k] * TMP[k][j];
//				}
//			}
//		}
//		delete[]Q;
//		Q = RES;
//		return;
//	}
//	bool check()
//	{
//		T** tmp = new T * [n];
//		for (i = 0; i < n; i++)
//		{
//			tmp[i] = new T[n];
//			for (j = 0; j < n; j++)
//				tmp[i][j] = 0;
//		}
//		for (i = 0; i < n; i++)
//			for (k = 0; k < n; k++)
//				for (j = 0; j < n; j++)
//					tmp[i][j] += Q[i][k] * R[k][j];
//		for (i = 0; i < n; i++)
//			for (j = 0; j < n; j++)
//				if (abs(REF_A[i][j] - tmp[i][j]) >= eps)
//					return false;
//		return true;
//	}
//	void out(char s)
//	{
//		if (s == 'R')
//			for (i = 0; i < n; i++)
//			{
//				for (j = 0; j < n; j++)
//					cout << R[i][j] << ' ';
//				cout << endl;
//			}
//		else if (s == 'Q')
//			for (i = 0; i < n; i++)
//			{
//				for (j = 0; j < n; j++)
//					cout << Q[i][j] << ' ';
//				cout << endl;
//			}
//		else if (s == 'A')
//			for (i = 0; i < n; i++)
//			{
//				for (j = 0; j < n; j++)
//					cout << A[i][j] << ' ';
//				cout << endl;
//			}
//	}
//	void transpQ()
//	{
//		for (i = 0; i < n; i++)
//			for (j = i + 1; j < n; j++)
//				swap(Q[i][j], Q[j][i]);
//		return;
//	}
//	~QR()
//	{
//		delete[]Q;
//		delete[]R;
//		delete[]A;
//		delete[]v;
//		Q = R = A = nullptr;
//		v = nullptr;
//		i = j = k = n = 0;
//	}
//};
//
//int main()
//{
//	size_t n;
//	cout << "Enter matrix dimentions: " << endl;
//	cin >> n;
//	QR <double> t(0, n);
//	//auto start{ chrono::steady_clock::now() };
//	//t.ort();
//	//t.norm();
//	//t.form_R();
//	//auto end{ chrono::steady_clock::now() };
//	////t.out('A');
//	//cout <<"QR-decomposition is finished, checking the result..." << endl;
//	////t.out('R');
//	//if (t.check())
//	//	cout << "QR-decomposition is correct";
//	//else cout << "QR-decomposition is incorrect";
//	//chrono::duration<double> elapsed_seconds = end - start;
//	//cout << endl << "Time spent: " << elapsed_seconds.count() << " sec" << endl;
//	//cout << endl << "Reference Q" << endl;
//	//t.out('Q');
//	//cout << endl << "Reference R" << endl;
//	//t.out('R');
//	//t.HHolder_R();
//	//cout << endl << "HHolder R" << endl;
//	//t.out('R');
//	//t.HHolder_Q();
//	//cout << endl << "HHolder Q" << endl;
//	//t.transpQ();
//	//t.out('Q');
//	auto start{ chrono::steady_clock::now() };
//	t.HHolder_R();
//	cout << endl;
//	t.out('R');
//	cout << endl;
//	t.HHolder_Q();
//	t.transpQ();
//	cout << endl;
//	t.out('Q');
//	cout << endl;
//	auto end{ chrono::steady_clock::now() };
//	chrono::duration<double> elapsed_seconds = end - start;
//	cout << "Time spent: " << elapsed_seconds.count() << " sec" << endl;
//	cout << "QR-decomposition is finished, checking the result..." << endl;
//	if (t.check())
//		cout << "QR-decomposition is correct";
//	else cout << "QR-decomposition is incorrect";
//	return 0;
//}
//
//
//+
//
//
//void ort()
//{
//	Q = new T * [n];
//	T tmp;
//	for (i = 0; i < n; i++)
//	{
//		Q[i] = new T[n];
//		//memcpy(Q[i], A[i], n);
//		for (j = 0; j < n; j++)
//			Q[i][j] = A[i][j];
//
//	}
//	for (j = 0; j < n; j++)
//	{
//		for (k = 0; k < j; k++)
//		{
//			tmp = scl_n(A, j, Q, k) / scl_sqr_n(Q, k);
//			for (i = 0; i < n; i++)
//				Q[i][j] -= tmp * Q[i][k];
//		}
//	}
//}
//void norm()
//{
//	T tmp;
//	for (size_t j = 0; j < n; j++)
//	{
//		tmp = sqrt(scl_sqr_n(Q, j));       //каким образом можно распознать нули, не вводя формальной проверки?
//		for (i = 0; i < n; i++)
//			Q[i][j] *= tmp;
//	}
//}
//void form_R()
//{
//	R = new T * [n];						//def Q=A
//	for (i = 0; i < n; i++)
//	{
//		R[i] = new T[n];
//		for (j = 0; j < n; j++)
//			memset(Q[i], 0, n * sizeof(T));            //почему неправильно?
//	}
//	for (i = 0; i < n; i++)
//		for (k = 0; k < n; k++)
//			for (j = 0; j < n; j++)
//				R[i][j] += Q[k][i] * A[k][j];
//}
//T scl_n(T** data_1, size_t ind_1_n, T** data_2, size_t ind_2_n)  //функция скалярного произведения для столбцов двух матриц одинакового размера
//{
//	T res = 0;
//	for (size_t i = 0; i < n; i++)
//	{
//		res += data_1[i][ind_1_n] * data_2[i][ind_2_n];
//	}
//	return res;
//}
//T scl_sqr_n(T** data, size_t ind)
//{
//	return scl_n(data, ind, data, ind);
//}
//
//void HHolder_Q()
//{
//	//#pragma omp parallel for private(i)
//	for (i = 0; i < n; i++)
//		Q[i * n + i] = TMP[i * n + i] = 1.0;
//
//	k = n - 1;
//	form_v(k);
//
//	gamma = -1 / gamma;
//	Q[k * n + k] = v[k] * v[k] * gamma + 1.0;
//
//	for (k = n - 2; k >= 0; k--)
//	{
//		form_v(k);
//		gamma = -1 / gamma;
//
//		for (i = k; i < n; i++)
//			//#pragma omp parallel for private(j)
//			for (j = k; j < n; j++)
//			{
//				TMP[i * n + j] = v[i] * v[j] * gamma;
//				if (i == j) TMP[i * n + j]++;
//			}
//
//		Q_mult_TMP_put_Q(k);
//
//		if (k == 0) break;
//	}
//
//}
//
//
//void Q_mult_TMP_put_Q(size_t ind)
//{
//	T* RES = new T[n * n]();
//
//#pragma omp parallel for private(i)
//	for (size_t i = 0; i < ind; i++)
//		RES[i * n + i] = 1.0;
//
//	for (size_t i = ind; i < n; i++)
//		for (size_t k = ind; k < n; k++)
//			//#pragma omp parallel for simd
//			for (size_t j = ind; j < n; j++)
//				RES[i * n + j] += Q[i * n + k] * TMP[k * n + j];
//
//	swap(RES, Q);
//
//}
//
//void HHolder_Block(size_t i_start, size_t b_size)
//{
//	for (j = i_start; j < i_start + b_size; j++)
//	{
//		form_v_gamma(j);
//
//#pragma omp parallel for simd private(k)       //??
//		for (k = j; k < n; k++)
//			factor[k - j] = scal(k, j) / gamma;
//
//#pragma omp parallel for simd private(i)       //??
//		for (i = j; i < n; i++)
//		{
//
//			//#pragma omp simd
//			for (k = j; k < n; k++)
//				R[i * n + k] -= v[i] * factor[k - j];
//		}
//	}
//}
//
//void form_v_gamma(size_t ind)
//{
//	//для типа T должен существовать конструктор по умолчанию;
//	T scl = 0;
//	//#pragma omp parallel for simd reduction(+: scl) 
//
//	for (size_t i = 0; i < n - ind; i++)
//	{
//		u[i] = R[(i + ind) * n + ind];
//		scl += u[i] * u[i];
//	}
//
//	if (scl < eps)
//	{
//		v[ind] = 1;
//		gamma = 0.5;
//		return;
//	}
//
//	else
//	{
//		scl = 1 / sqrt(scl);
//		u[0] *= scl;
//		gamma = (1 + abs(u[0]));
//		v[ind] = sgn(u[0]) * gamma;
//
//		//#pragma omp parallel for simd
//		for (size_t i = ind + 1; i < n; i++)
//			v[i] = u[i - ind] * scl;
//
//		return;
//	}
//}