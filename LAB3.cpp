#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;


//Функция транспонирования матрицы
template <typename T> void TransponMtx(T** matr, T** tMatr, int n) {
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			tMatr[j][i] = matr[i][j];
}


//Функция освобождения памяти
template <typename T> void FreeMem(T** matr, int n)
{
	for (int i = 0; i < n; i++)
		delete[] matr[i];
	delete[] matr;
}


//Функция заполнения матрицы
template <typename T> void SetMtx(T** matr, vector<vector<int>> M, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {

			matr[i][j] = M[i][j];
		}
}


//Функция печати матрицы
template <typename T> void PrintMtx(T** matr, int n)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			cout << setw(20) << setprecision(4) << matr[i][j] << " ";
		cout << endl;
	}
}


//Функция вычеркивания строки и столбца
void Get_matr(int** matr, int n, int** temp_matr, int indRow, int indCol)
{
	int ki = 0;
	for (int i = 0; i < n; i++) {
		if (i != indRow) {
			for (int j = 0, kj = 0; j < n; j++) {
				if (j != indCol) {
					temp_matr[ki][kj] = matr[i][j];
					kj++;
				}
			}
			ki++;
		}
	}
}


//Функция вычисления определителя матрицы
int Det(int** matr, int n)
{
	int temp = 0;   //Временная переменная для хранения определителя
	int k = 1;      //Степень
	if (n < 1) {
		cout << "Ошибка ввода: некорректный размер матрицы" << endl;
		return 0;
	}
	else if (n == 1)
		temp = matr[0][0];
	else if (n == 2)
		temp = matr[0][0] * matr[1][1] - matr[1][0] * matr[0][1];
	else {
		for (int i = 0; i < n; i++) {
			int m = n - 1;
			int** temp_matr = new int*[m];
			for (int j = 0; j < m; j++)
				temp_matr[j] = new int[m];
			Get_matr(matr, n, temp_matr, 0, i);
			temp = temp + k * matr[0][i] * Det(temp_matr, m);
			k = -k;
			FreeMem(temp_matr, m);
		}
	}
	return temp;
}


void matrFill(vector<vector<int>> &M, int n) { //Случайное заполнение матриц
	for (int i = 0; i < n; i++) {
		vector<int> temp;
		for (int j = 0; j < n; j++) {
			int znak = rand() % 2;
			if (znak == 0) temp.push_back(rand() % 50);
			if (znak == 1) temp.push_back((rand() % 50) * -1);
		}
		M.push_back(temp);
	}
}


void bimatrFill(vector<vector<int>> &TA, vector<vector<int>> &TB, int n) {
	for (int i = 0; i < 2; i++) {
		vector<int> temp;
		for (int j = 0; j < 2; j++) {
			temp.push_back(0);
		}
		TA.push_back(temp);
		TB.push_back(temp);
	}
	TA[0][0] = 6; 
	TB[0][0] = 7; 
	TA[0][1] = 8;
	TB[0][1] = 4;
	TA[1][0] = 2;
	TB[1][0] = 1; 
	TA[1][1] = 9; 
	TB[1][1] = 3;
}




void matrZeroFill(vector<vector<int>>& A, vector<vector<int>>& B, int n) { //Создание матриц для нахождения пересечения множеств 
	for (int i = 0; i < n; i++) {
		vector<int> temp;
		for (int j = 0; j < n; j++) {
			temp.push_back(0);
		}
		A.push_back(temp);
		B.push_back(temp);
	}
}


void matrPrint(vector<vector<int>> A, vector<vector<int>> B, int n) { //Вывод матрицы стратегий двух игроков

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << "(" << setw(3) << A[i][j] << ", " << setw(3) << B[i][j] << ") ";
		}
		cout << endl;
	}
}


int checkMaxStr(vector<vector<int>> M, int n, int x) { //Поиск максимума по строке
	int max = -100;
	int m;
	int unique = 0;


	for (int j = 0; j < n; j++) {


		if (M[x][j] > max) {
			max = M[x][j];
			m = j;
		}


	}


	for (int j = 0; j < n; j++) {


		if (M[x][j] == max) { //Проверка на уникальность максимума
			unique++;
		}


	}


	if (unique > 1) {
		//cout << "Для строки " << x + 1 << " уникальность максимума не достигнута (" << max << ")" << endl;
		return n;
	}
	return m;
}


int checkMaxStl(vector<vector<int>> M, int n, int x) { //Поиск максимума по столбцу
	int max = -100;
	int m;
	int unique = 0;


	for (int j = 0; j < n; j++) {


		if (M[j][x] > max) {
			max = M[j][x];
			m = j;
		}


	}

	for (int j = 0; j < n; j++) {


		if (M[j][x] == max) { //Проверка на уникальность максимума
			unique++;
		}


	}

	if (unique > 1) {
		//cout << "Для столбца " << x + 1 << " уникальность максимума не достигнута (" << max << ")" << endl;
		return n;
	}


	return m;
}


void nashCheck(vector<vector<int>> A, vector<vector<int>> B, vector<vector<int>> &tN, int n) {
	vector<vector<int>> nash;


	for (int i = 0; i < n; i++) {
		vector<int> temp;
		for (int j = 0; j < n; j++) {
			temp.push_back(0);
		}
		nash.push_back(temp);
	}


	for (int i = 0; i < n; i++) {
		if (checkMaxStl(A, n, i) != n) {
			//cout << "Макс. в столбце " << i + 1 << " для игрока А: " << A[checkMaxStl(A, n, i)][i] << endl;
			nash[checkMaxStl(A, n, i)][i]++;
		}

		if (checkMaxStr(B, n, i) != n) {
			//cout << "Макс. в строке " << i + 1 << " для игрока В: " << B[i][checkMaxStr(B, n, i)] << endl;
			nash[i][checkMaxStr(B, n, i)]++;
		}

	}


	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (nash[i][j] == 2) cout << setw(1) << " + ";
			if (nash[i][j] == 1) cout << setw(1) << " = ";
			if (nash[i][j] == 0) cout << setw(1) << " - ";
		}
		cout << endl;
	}


	cout << "Равновесие Нэша: " << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (nash[i][j] == 2) cout << "(" << setw(3) << A[i][j] << ", " << setw(3) << B[i][j] << ") " << endl;
		}

	}


	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			tN[i][j] = nash[i][j];
		}
	}


}


void matrExamples(vector<vector<int>> &Pa, vector<vector<int>> &Sa, vector<vector<int>> &Da,
	vector<vector<int>> &Pb, vector<vector<int>> &Sb, vector<vector<int>> &Db) { //Заполнение матриц-примеров
	for (int i = 0; i < 2; i++) {

		vector<int> temp;
		for (int j = 0; j < 2; j++) {
			temp.push_back(0);
		}


		Pa.push_back(temp);
		Pb.push_back(temp);
		Sa.push_back(temp);
		Sb.push_back(temp);
		Da.push_back(temp);
		Db.push_back(temp);


	}


	Pa[0][0] = -1; Pb[0][0] = -1; Pa[0][1] = -2; Pb[0][1] = 2;
	Pa[1][0] = 2; Pb[1][0] = -5; Pa[1][1] = -7; Pb[1][1] = -7;


	Sa[0][0] = 4; Sb[0][0] = 1; Sa[0][1] = 0; Sb[0][1] = 0;
	Sa[1][0] = 0; Sb[1][0] = 0; Sa[1][1] = 1; Sb[1][1] = 4;


	Da[0][0] = -5; Db[0][0] = -5; Da[0][1] = 0; Db[0][1] = -10;
	Da[1][0] = -10; Db[1][0] = 0; Da[1][1] = -1; Db[1][1] = -1;


}


void paretoCheck(vector<vector<int>> A, vector<vector<int>> B, vector<vector<int>> &tP, int n) {

	vector<vector<int>> pareto;
	for (int i = 0; i < n; i++) {
		vector<int> temp;
		for (int j = 0; j < n; j++) {
			temp.push_back(0);
		}
		pareto.push_back(temp);
	}


	//int c = 0;
	for (int i = 0; i < n; i++) {


		for (int j = 0; j < n; j++) {


			for (int i2 = 0; i2 < n; i2++) {


				for (int j2 = 0; j2 < n; j2++) {



					if (((A[i][j] <= A[i2][j2]) && (B[i][j] < B[i2][j2])) ||
						((A[i][j] < A[i2][j2]) && (B[i][j] <= B[i2][j2]))) {


						pareto[i][j]++;


					}


				}


			}


		}


	}


	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (pareto[i][j] > 0) cout << setw(1) << " - ";
			if (pareto[i][j] == 0) cout << setw(1) << " + ";
		}
		cout << endl;
	}


	cout << "Оптимальные по Парето: " << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (pareto[i][j] == 0) cout << "(" << setw(3) << A[i][j] << ", " << setw(3) << B[i][j] << ") " << endl;
		}
	}


	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			tP[i][j] = pareto[i][j];
		}
	}


}


void findCol(vector<vector<int>> A, vector<vector<int>> B, vector<vector<int>> tN, vector<vector<int>> tP, int n) {
	vector<vector<int>> col;
	int c = 0;


	for (int i = 0; i < n; i++) {
		vector<int> temp;
		for (int j = 0; j < n; j++) {
			temp.push_back(0);
		}
		col.push_back(temp);
	}


	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if ((tN[i][j] == 2) && (tP[i][j] == 0)) {
				col[i][j]++;
				c++;
			}
		}
	}


	if (c == 0) {
		cout << "Нет пересечения множеств" << endl;
	}
	else {
		cout << "Пересечение:" << endl;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (col[i][j] == 1) {
					cout << "(" << setw(3) << A[i][j] << ", " << setw(3) << B[i][j] << ") " << endl;
				}
			}
		}
	}
}


template <typename T> void FindMtx(T** tobr_matr, int n) {


	vector<float> res; 
	for (int i = 0; i < n; i++) {
		res.push_back(0);
	}


	vector<float> u;
	for (int i = 0; i < n; i++) {
		u.push_back(1);
	}


	for (int i = 0; i < n; i++)
	{
		float temp = 0;
		for (int j = 0; j < n; j++)
		{
			temp += tobr_matr[i][j] * u[i];
		}
		res[i] = temp;
	}


	cout << "uB^-1 = (";
	for (int i = 0; i < n; i++) {
		cout << res[i]; if (i != n - 1) cout << ", ";
	}
	cout << ") " << endl;


	float temp = 0; for (int i = 0; i < n; i++) temp += res[i];
	cout << "v1 = uB^-1u = " << 1 / temp << endl;


	cout << "x = [";
	for (int i = 0; i < n; i++) {
		cout << res[i] / temp; if (i != n - 1) cout << ", "; 
	}
	cout << "] " << endl;
}


template <typename T> void FindMty(T** tobr_matr, int n) {


	vector<float> res;
	for (int i = 0; i < n; i++) {
		res.push_back(0);
	}


	vector<float> u; 
	for (int i = 0; i < n; i++) {
		u.push_back(1);
	}


	for (int i = 0; i < n; i++)
	{
		float temp = 0;
		for (int j = 0; j < n; j++)
		{
			temp += tobr_matr[j][i] * u[i];
		}
		res[i] = temp;
	}


	cout << "A^-1u = (";
	for (int i = 0; i < n; i++) { 
		cout << res[i]; if (i != n - 1) cout << ", "; 
	}
	cout << ") " << endl;


	float temp = 0;
	for (int i = 0; i < n; i++) 
		temp += res[i];
	cout << "v2 = uA^-1u = " << 1 / temp << endl;


	cout << "y = [";
	for (int i = 0; i < n; i++) { 
		cout << res[i] / temp; if (i != n - 1) cout << ", "; 
	}
	cout << "] " << endl;
}




int main()
{
	srand(unsigned(time(0)));
	setlocale(0, "");
	cout << "================== Задание 1 ==================" << endl;
	int n = 10;
	vector<vector<int>> A; 
	vector<vector<int>> B; //Матрицы стратегий игроков А и В
	vector<vector<int>> PA, PB, SA, SB, DA, DB, tA, tB, tNas, tPar;


	matrFill(A, n);
	matrFill(B, n); //Случайное заполнение матриц стратегий игроков А и В
	matrPrint(A, B, n); //Вывод матрицы стратегий игроков А и В


	matrZeroFill(tNas, tPar, n); //Создание матриц для нахождения пересечения множеств


	cout << endl << "----------- Поиск равновесия Нэша -----------" << endl;
	nashCheck(A, B, tNas, n);
	cout << endl << "---------- Оптимальность по Парето ----------" << endl;
	paretoCheck(A, B, tPar, n);
	cout << endl << "---------- Нахождение пересечений ----------" << endl;
	findCol(A, B, tNas, tPar, n);


	cout << endl << "====== Проверка алгоритмов на примерах ======" << endl;
	matrExamples(PA, SA, DA, PB, SB, DB);

	cout << endl << "==== Перекресток (с разными отклонениями) ====" << endl;
	matrPrint(PA, PB, 2);
	cout << endl << "----------- Поиск равновесия Нэша -----------" << endl;
	nashCheck(PA, PB, tNas, 2);
	cout << endl << "---------- Оптимальность по Парето ----------" << endl;
	paretoCheck(PA, PB, tPar, 2);
	cout << endl << "---------- Нахождение пересечений ----------" << endl;
	findCol(PA, PB, tNas, tPar, 2);

	cout << endl << "=============== Семейный спор ===============" << endl;
	matrPrint(SA, SB, 2);
	cout << endl << "----------- Поиск равновесия Нэша -----------" << endl;
	nashCheck(SA, SB, tNas, 2);
	cout << endl << "---------- Оптимальность по Парето ----------" << endl;
	paretoCheck(SA, SB, tPar, 2);
	cout << endl << "---------- Нахождение пересечений ----------" << endl;
	findCol(SA, SB, tNas, tPar, 2);


	cout << endl << "=========== Дилемма заключенного ===========" << endl;
	matrPrint(DA, DB, 2);
	cout << endl << "----------- Поиск равновесия Нэша -----------" << endl;
	nashCheck(DA, DB, tNas, 2);
	cout << endl << "---------- Оптимальность по Парето ----------" << endl;
	paretoCheck(DA, DB, tPar, 2);
	cout << endl << "---------- Нахождение пересечений ----------" << endl;
	findCol(DA, DB, tNas, tPar, 2);


	cout << endl << "=============== Задание 2 ==================" << endl;
	bimatrFill(tA, tB, 2);
	matrPrint(tA, tB, 2);
	cout << endl << "----------- Поиск равновесия Нэша -----------" << endl;
	nashCheck(tA, tB, tNas, 2);
	cout << endl << "--------- Нахождение обратных матриц ---------" << endl;

	int det; n = 2; //Определитель и размер матрицы

	int** matr = new int*[n]; 
	int** matrb = new int*[n]; //Матрицы значений
	double** obr_matr = new double*[n];
	double** obr_matrb = new double*[n]; //Обратные матрицы
	double** tobr_matr = new double*[n]; 
	double** tobr_matrb = new double*[n]; //Транспонированные обратные матрицы
	for (int i = 0; i < n; i++) {
		matr[i] = new int[n]; 
		matrb[i] = new int[n];
		obr_matr[i] = new double[n]; 
		obr_matrb[i] = new double[n];
		tobr_matr[i] = new double[n]; 
		tobr_matrb[i] = new double[n];
	}
	SetMtx(matr, tA, n); SetMtx(matrb, tB, n);
	cout << "Матрица А: " << endl;
	PrintMtx(matr, n);
	det = Det(matr, n);
	cout << "Определитель матрицы = " << det << endl;
	if (det) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				int m = n - 1;
				int** temp_matr = new int*[m];
				for (int k = 0; k < m; k++)
					temp_matr[k] = new int[m];
				Get_matr(matr, n, temp_matr, i, j);
				obr_matr[i][j] = pow(-1.0, i + j + 2) * Det(temp_matr, m) / det;
				FreeMem(temp_matr, m);
			}
		}
	}
	else cout << "Определитель матрицы равен 0 => матрица является вырожденной и обратной нет" << endl;


	//Транспонирование матрицы
	TransponMtx(obr_matr, tobr_matr, n);
	//Печать обратной матрицы после транспонирования
	cout << "Обратная матрица: " << endl;
	PrintMtx(tobr_matr, n);

	cout << "Матрица В: " << endl;
	PrintMtx(matrb, n);
	det = Det(matrb, n);
	cout << "Определитель матрицы = " << det << endl;
	if (det) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				int m = n - 1;
				int** temp_matrb = new int*[m];
				for (int k = 0; k < m; k++)
					temp_matrb[k] = new int[m];
				Get_matr(matrb, n, temp_matrb, i, j);
				obr_matrb[i][j] = pow(-1.0, i + j + 2) * Det(temp_matrb, m) / det;
				FreeMem(temp_matrb, m);
			}
		}
	}
	else cout << "Определитель матрицы равен 0 => матрица является вырожденной и обратной нет" << endl;


	//Транспонирование матрицы
	TransponMtx(obr_matrb, tobr_matrb, n);
	//Печать обратной матрицы после транспонирования
	cout << "Обратная матрица: " << endl;
	PrintMtx(tobr_matrb, n);

	FindMtx(obr_matrb, n);
	FindMty(obr_matr, n);

	//Освобождение памяти
	FreeMem(tobr_matr, n);
	FreeMem(matr, n);
	FreeMem(obr_matr, n);
	FreeMem(tobr_matrb, n);
	FreeMem(matrb, n);
	FreeMem(obr_matrb, n);
	system("Pause");
}
