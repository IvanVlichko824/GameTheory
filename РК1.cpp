#include <iostream>

#include <iomanip>

#include <vector>

using TMatrix = std::vector<std::vector<double>>;

static const TMatrix G_MATRIX = {

{ 6, 15, 16 },

{ 15, 10, 0 },

{ 10, 7, 10}

};

TMatrix GetCoFactor(const TMatrix& matrix, size_t row, size_t col) {

	TMatrix minorMatrix(matrix.size() - 1, std::vector<double>(matrix[0].size() - 1));

	for (size_t i = 0; i < matrix.size(); ++i) {

		if (i != row) {

			size_t _i = i;

			if (i > row) {

				_i -= 1;

			}

			for (size_t j = 0; j < matrix[i].size(); ++j) {

				if (j != col) {

					size_t _j = j;

					if (j > col) {

						_j -= 1;

					}

					minorMatrix[_i][_j] = matrix[i][j];

				}

			}

		}

	}

	return minorMatrix;

}

double Determinant(const TMatrix& matrix) {

	double result = 0;

	if (matrix.size() == 1) {

		return matrix[0][0];

	}

	int sign = 1;

	for (size_t j = 0; j < matrix.size(); ++j) {

		TMatrix coFactor = GetCoFactor(matrix, 0, j);

		result += sign * matrix[0][j] * Determinant(coFactor);

		sign = -sign;

	}

	return result;

}

TMatrix GetTransposeMatrix(const TMatrix& matrix) {

	TMatrix transposed(matrix[0].size(), std::vector<double>(matrix.size()));

	for (size_t i = 0; i < matrix.size(); ++i) {

		for (size_t j = 0; j < matrix[i].size(); ++j) {

			transposed[j][i] = matrix[i][j];

		}

	}

	return transposed;

}

TMatrix GetCoFactorsDMatrix(const TMatrix& matrix) {

	TMatrix coFactorDMatrix = matrix;

	size_t k = 0;

	for (size_t i = 0; i < matrix.size(); ++i) {

		for (size_t j = 0; j < matrix.size(); ++j) {

			TMatrix coFactor = GetCoFactor(matrix, i, j);

			coFactorDMatrix[i][j] = Determinant(coFactor);

			++k;

			if ((k % 2) == 0) {

				coFactorDMatrix[i][j] = -coFactorDMatrix[i][j];

			}

		}

	}

	return coFactorDMatrix;

}

TMatrix MultiplyNumber(const TMatrix& matrix, double value) {

	TMatrix result = matrix;

	for (size_t i = 0; i < matrix.size(); ++i) {

		for (size_t j = 0; j < matrix[i].size(); ++j) {

			result[i][j] = matrix[i][j] * value;

		}

	}

	return result;

}

TMatrix MultiplyMatrix(const TMatrix& a, const TMatrix& b) {

	_ASSERT(a[0].size() == b.size());

	TMatrix result(a.size(), std::vector<double>(b[0].size()));

	for (size_t i = 0; i < result.size(); ++i) {

		for (size_t j = 0; j < result[i].size(); ++j) {

			result[i][j] = 0;

			for (size_t k = 0; k < a[i].size(); ++k) {

				result[i][j] += a[i][k] * b[k][j];

			}

		}

	}

	return result;

}

TMatrix GetReverseMatrix(const TMatrix& matrix) {

	double d = Determinant(matrix);

	return MultiplyNumber(GetTransposeMatrix(GetCoFactorsDMatrix(G_MATRIX)), 1.0 / d);

}

std::ostream& operator<<(std::ostream& out, const TMatrix& v) {

	for (auto&& i : v) {

		for (auto&& j : i) out << j << " ";

		out << std::endl;

	}

	return out;

}

std::ostream& operator<<(std::ostream& out, const std::vector<double>& v) {

	for (auto&& j : v) out << j << " ";

	out << std::endl;

	return out;

}

std::vector<double> GetRow(const TMatrix& matrix, size_t rowNumber) {

	std::vector<double> row(matrix[0].size());

	for (size_t j = 0; j < row.size(); ++j)

		row[j] = matrix[rowNumber][j];

	return row;

}

std::vector<double> GetColumn(const TMatrix& matrix, size_t colNumber) {

	std::vector<double> col(matrix.size());

	for (size_t i = 0; i < col.size(); ++i)

		col[i] = matrix[i][colNumber];

	return col;

}

size_t GetMinPos(const std::vector<double>& a) {

	double min = a[0];

	size_t pos = 0;

	for (size_t i = 1; i < a.size(); ++i) {

		if (a[i] < min) {

			min = a[i];

			pos = i;

		}

	}

	return pos;

}

size_t GetMaxPos(const std::vector<double>& a) {

	double max = a[0];

	size_t pos = 0;

	for (size_t i = 1; i < a.size(); ++i) {

		if (a[i] > max) {

			max = a[i];

			pos = i;

		}

	}

	return pos;

}

void OutputData(size_t k, size_t a, size_t b, const std::vector<double>& aw,

	const std::vector<double>& bl, double vup, double vdown, double eps, double v) {

	std::cout << "|" << std::setw(4) << k << "|"

		<< std::setw(4) << a + 1 << "|" << std::setw(4) << b + 1 << "|";

	for (size_t i = 0; i < aw.size(); ++i) {

		std::cout << std::setw(5) << (int)aw[i] << "|";

	}

	for (size_t i = 0; i < bl.size(); ++i) {

		std::cout << std::setw(5) << (int)bl[i] << "|";

	}

	std::cout << std::setw(8) << std::setprecision(4) << vup << "|"

		<< std::setw(8) << std::setprecision(4) << vdown << "|"

		<< std::setw(8) << std::setprecision(4) << eps << "|"

		<< std::setw(8) << std::setprecision(4) << v << "|" << std::endl;

}

void BrownRobinson(const TMatrix& game, TMatrix& as, TMatrix& bs, double& v, const double game_eps) {

	size_t k = 0; // счетчик

	std::vector<double> aw(game.size()); // выигрыш A

	as = TMatrix(1, std::vector<double>(game.size())); // частота стратегий A

	std::vector<double> bl(game[0].size()); // проигрыш B

	bs = TMatrix(1, std::vector<double>(game[0].size())); // частота стратегий B

	size_t a = 0; // стратегия игрока А

	size_t b = 0; // стратегия игрока B

	double eps = 1; // ошибка

	double vUp = 0;

	double vDown = 0;

	double minvUp = 1000000;

	double maxvDown = 0;

	v = 0;

	while (eps > game_eps) {

		++k;

		std::vector<double> curAw = GetColumn(game, b);

		std::vector<double> curBl = GetRow(game, a);

		for (size_t i = 0; i < curAw.size(); ++i)

			aw[i] += curAw[i];

		for (size_t i = 0; i < curBl.size(); ++i)

			bl[i] += curBl[i];

		size_t nA = GetMaxPos(aw);

		size_t nB = GetMinPos(bl);

		vUp = aw[nA] / (double)k;

		vDown = bl[nB] / (double)k;

		if (vUp < minvUp)

			minvUp = vUp;

		if (vDown > maxvDown)

			maxvDown = vDown;

		eps = minvUp - maxvDown;

		v = (minvUp + maxvDown) / 2.0;

		//v = (vUp + vDown) / 2.0;

		as[0][a] += 1;

		bs[0][b] += 1;

		//OutputData(k, a, b, aw, bl, vUp, vDown, eps, v);

		a = nA;

		b = nB;

	}

	as = MultiplyNumber(as, 1.0 / (double)k);

	bs = MultiplyNumber(bs, 1.0 / (double)k);

	std::cout << "x* = " << as << "y* = " << bs << "v = " << v << std::endl;

}

void Monotone(const TMatrix& game, const double eps, const double sub_game_eps)

{

	int m = game.size();

	int n = game[0].size();

	std::vector<double> x(m);

	std::vector<double> c(n);

	double v;

	x[0] = 1; //x0 = (1,0,....,0) - первая стратегия

	c = GetRow(game, 0); //c0 = a0 - первой строке в матрице игры

	double alpha = 1;

	int i = 1;

	std::cout << "Game:" << std::endl << game << std::endl;

	do

	{

		std::cout << "Iteration " << i << std::endl;

		std::vector<int> jk;

		v = c[0];

		for (int i = 1; i < n; ++i)

		{

			if (c[i] < v)

			{

				v = c[i];

			}

		}

		for (int i = 0; i < n; ++i)

		{

			if (c[i] == v)

			{

				jk.push_back(i);

			}

		}

		TMatrix sub_game;

		for (int i = 0; i < m; i++)

		{

			sub_game.push_back(std::vector<double>());

			for (int j = 0; j < jk.size(); ++j)

			{

				sub_game[i].push_back(game[i][jk[j]]);

			}

		}

		TMatrix x_sub, y_sub;

		double v_sub;

		std::cout << "Subgame1:" << std::endl << sub_game << std::endl;

		BrownRobinson(sub_game, x_sub, y_sub, v_sub, sub_game_eps);

		std::cout << "-----------------------" << std::endl;

		std::vector<double> c_new(n);

		for (int j = 0; j < n; ++j)

		{

			c_new[j] = 0;

			for (int i = 0; i < m; ++i)

			{

				c_new[j] += x_sub[0][i] * game[i][j];

			}

		}

		TMatrix sub_game2;

		sub_game2.push_back(std::vector<double>());

		sub_game2.push_back(std::vector<double>());

		for (int j = 0; j < n; ++j)

		{

			sub_game2[0].push_back(c_new[j]); // Перепутаны эти строки в учебнике Петросяна

		}

		for (int j = 0; j < n; ++j)

		{

			sub_game2[1].push_back(c[j]); // Перепутаны эти строки в учебнике Петросяна

		}

		TMatrix x_sub2, y_sub2;

		double v_sub2;

		std::cout << "Subgame2:" << std::endl << sub_game2 << std::endl;

		BrownRobinson(sub_game2, x_sub2, y_sub2, v_sub2, sub_game_eps);

		std::cout << "-----------------------" << std::endl;

		double alpha_new = x_sub2[0][0];

		for (int i = 0; i < m; ++i)

		{

			x[i] = (1 - alpha_new) * x[i] + alpha_new * x_sub[0][i];

		}

		for (int j = 0; j < n; ++j)

		{

			c[j] = (1 - alpha_new) * c[j] + alpha_new * c_new[j];

		}

		alpha = alpha_new;

		std::cout << "Iteration results:" << std::endl;

		std::cout << "alpha = " << alpha << std::endl;

		std::cout << "x* = " << x << "v* = " << v << std::endl;

		std::cout << "------------------------------------------------" << std::endl;

		i++;

	} while (alpha > eps);

	std::cout << "x = " << x << "v = " << v << std::endl;

}

int main()

{

	// монотонный метод

	std::cout << "Monotone method" << std::endl << std::endl;

	std::cout << "+--------------------------------------------------------------------------------------+" << std::endl;

	std::cout << std::fixed << std::setprecision(2); // точность округления вывода

	Monotone(G_MATRIX, 0.001, 0.0001); // передаем матрицу игры, точность решения в монотонном методе и точность решения подыгр методом Брауна-Робинсона

	system("Pause");
}