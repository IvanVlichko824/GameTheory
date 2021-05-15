#include <iostream>
#include <vector>
#include <functional>
#include <random>
#include <iomanip>
#include <chrono>

using namespace std;
using TMatrix = std::vector<std::vector<double>>;

class RandomDouble
{
public:
	RandomDouble(double low, double high):r(std::bind(std::uniform_real_distribution<>(low, high), std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count())))
	{}
	double operator()() {
		return r();
	}
private:
	std::function<double()> r;
};

std::ostream& operator<<(std::ostream& out, const TMatrix& v) {
	out << std::fixed << std::setprecision(3);
	for (auto&& i : v) {
		for (auto&& j : i) out << std::setw(7) << j << " ";
		out << std::endl;
	}
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

double SumVector(const std::vector<double>& vector) {
	double sum = 0;
	for (size_t i = 0; i < vector.size(); ++i) {
		sum += vector[i];
	}
	return sum;
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

TMatrix GetTransposeMatrix(const TMatrix& matrix) {
	TMatrix transposed(matrix[0].size(), std::vector<double>(matrix.size()));
	for (size_t i = 0; i < matrix.size(); ++i) {
		for (size_t j = 0; j < matrix[i].size(); ++j) {
			transposed[j][i] = matrix[i][j];
		}
	}
	return transposed;
}

void InfluenceMatrix(size_t agentsCount, TMatrix& result) {
	RandomDouble rd{ 0, 1 };
	result = TMatrix(agentsCount, std::vector<double>(agentsCount));
	for (size_t i = 0; i < agentsCount; ++i) {
		double rowSum = 0;
		for (size_t j = 0; j < agentsCount; ++j) {
			result[i][j] = rd();
			rowSum += result[i][j];
		}
		for (size_t j = 0; j < agentsCount; ++j) {
			result[i][j] /= rowSum;
		}
	}
}

void OpinionsVector(size_t agentsCount, TMatrix& result) {
	result = TMatrix(1, std::vector<double>(agentsCount));
	RandomDouble rd{ 1, 20 };
	for (size_t i = 0; i < agentsCount; ++i) {
		result[0][i] = rd();
	}
}

bool IsEnoughPrecision(const TMatrix& infMatrix, const TMatrix& opinVector, const double eps) {
	for (size_t i = 0; i < infMatrix.size(); ++i) {
		for (size_t j = 0; j < infMatrix.size(); ++j) {
			if (std::abs(infMatrix[i][j] - infMatrix[0][j]) > eps)
				return false;
		}
	}

	for (size_t i = 0; i < opinVector.size(); ++i) {
		if (std::abs(opinVector[0][0] - opinVector[i][0]) > eps)
			return false;
	}
	return true;
}

void BoughtAgentsVector(TMatrix& result, const double fO, const double sO) {
	RandomDouble rd{ 0, 3 };
	TMatrix firstBoughts = TMatrix(1, std::vector<double>());
	TMatrix secondBoughts = TMatrix(1, std::vector<double>());
	for (size_t i = 0; i < result[0].size(); ++i) {
		double o = rd();
		if (((int)o) % 3 == 0) {
			result[0][i] = fO;
			firstBoughts[0].push_back(i + 1);
		}
		else if (((int)o) % 3 == 1) {
			result[0][i] = sO;
			secondBoughts[0].push_back(i + 1);
		}
	}

	std::cout << "First player influence agents: " << firstBoughts;
	std::cout << "Second player influence agents: " << secondBoughts;
}

void Solution(const TMatrix& infMatrix, const TMatrix& opinVector, const double eps) {
	TMatrix inf = infMatrix;
	TMatrix opin = GetTransposeMatrix(opinVector);
	while (!IsEnoughPrecision(inf, opin, eps)) {
		inf = MultiplyMatrix(inf, inf);
		opin = MultiplyMatrix(infMatrix, opin);
	}

	std::cout << inf << endl;
	std::cout << "Agents final opinions vector" << endl;
	std::cout << GetTransposeMatrix(opin);
}

int main()

{
	TMatrix infMatrix;
	TMatrix opinVector;
	size_t numAgents = 10;
	double eps = 1e-6;
	InfluenceMatrix(numAgents, infMatrix);
	std::cout << "Relations matrix:" << std::endl << infMatrix << "-----------------------------" << std::endl;
	OpinionsVector(numAgents, opinVector);
	std::cout << "Opinions vector:" << std::endl << opinVector << "-----------------------------" << std::endl;
	std::cout << "Solution:" << std::endl;
	Solution(infMatrix, opinVector, eps);
	std::cout << "-----------------------------" << std::endl;
	std::cout << "Generating influence game:" << std::endl;
	RandomDouble fbf = { -100, 0 };
	RandomDouble sbf = { 0, 100 };
	int firstBoughtInf = (int)fbf();
	int secondBoughtInf = (int)sbf();
	BoughtAgentsVector(opinVector, firstBoughtInf, secondBoughtInf);
	std::cout << "New opinions vector:" << std::endl << opinVector << "-----------------------------" << std::endl;
	std::cout << "Solution:" << std::endl;
	Solution(infMatrix, opinVector, eps);
	std::cout << "-----------------------------" << std::endl;
	system("Pause");
}