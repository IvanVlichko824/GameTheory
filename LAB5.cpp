#include <iostream>
#include <set>
#include <algorithm>
#include <iterator>
#include <list>
#include "locale.h"

using namespace std;

int factorial(int i)
{
	if (i == 0)
		return 1;
	else
		return i * factorial(i - 1);
}

set<int> set0 = set<int>();
set<int> set1 = { 1 };
set<int> set2 = { 2 };
set<int> set3 = { 3 };
set<int> set4 = { 4 };
set<int> set12 = { 1, 2 };
set<int> set13 = { 1, 3 };
set<int> set14 = { 1, 4 };
set<int> set23 = { 2, 3 };
set<int> set24 = { 2, 4 };
set<int> set34 = { 3, 4 };
set<int> set123 = { 1, 2, 3 };
set<int> set124 = { 1, 2, 4 };
set<int> set134 = { 1, 3, 4 };
set<int> set234 = { 2, 3, 4 };
set<int> set1234 = { 1, 2, 3, 4 };

list<set<int>> coalitions = { set0, set1, set2, set3, set4, set12, set13,
	set14, set23, set24, set34, set123, set124, set134, set234, set1234 };

int v(set<int> s)
{
	if (s == set0)
		return 0;
	else if (s == set1)
		return 4;
	else if (s == set2)
		return 1;
	else if (s == set3)
		return 3;
	else if (s == set4)
		return 1;
	else if (s == set12)
		return 6;
	else if (s == set13)
		return 8;
	else if (s == set14)
		return 6;
	else if (s == set23)
		return 5;
	else if (s == set24)
		return 3;
	else if (s == set34)
		return 5;
	else if (s == set123)
		return 9;
	else if (s == set124)
		return 8;
	else if (s == set134)
		return 10;
	else if (s == set234)
		return 7;
	else if (s == set1234)
		return 11;
	return 0;
}

void printSet(set<int> S)
{
	if (S.empty())
	{
		cout << "{}";
		return;
	}
	cout << "{";
	auto it = S.begin();
	cout << *it;
	++it;
	for (; it != S.end(); ++it)
		cout << ", " << *it;
	cout << "}";
}


bool check_superadditivity()
{
	for (auto S : coalitions)
	{
		for (auto T : coalitions)
		{
			set<int> intersect_set;
			set_intersection(S.begin(), S.end(), T.begin(), T.end(),
				inserter(intersect_set, intersect_set.begin()));
			set<int> union_set = S;
			union_set.insert(T.begin(), T.end());
			if (intersect_set.empty() && v(union_set) < v(S) + v(T))
			{
				cout << "Игра не является супераддитивной. S: ";
				printSet(S);
				cout << ", T: ";
				printSet(T);
				cout << endl;
				return false;
			}
		}
	}
	cout << "Игра является супераддитивной" << endl;
	return true;
}


bool check_convexity()
{
	for (auto S : coalitions)
	{
		for (auto T : coalitions)
		{
			set<int> intersect_set;
			set_intersection(S.begin(), S.end(), T.begin(), T.end(),
				inserter(intersect_set, intersect_set.begin()));
			set<int> union_set = S;
			union_set.insert(T.begin(), T.end());
			if (v(union_set) + v(intersect_set) < v(S) + v(T))
			{
				cout << "Игра не является выпуклой. S: ";
				printSet(S);
				cout << ", T: ";
				printSet(T);
				cout << endl;
				return false;
			}
		}
	}
	cout << "Игра является выпуклой" << endl;
	return true;
}

void find_Shapley()
{
	const int N = 4;
	double X[N];

	for (int i = 1; i <= N; ++i)
	{
		int sum = 0;
		for (auto S : coalitions)
		{
			set<int> withoutI = S;
			auto it = withoutI.find(i);
			if (it != withoutI.end())
			{
				withoutI.erase(i);
				sum += factorial(S.size() - 1) * factorial(N - S.size()) * (v(S) - v(withoutI));
			}
		}
		X[i - 1] = double(sum) / factorial(N);
	}

	double sum = 0;
	cout << "Вектор Шепли: [";
	for (int i = 0; i < N - 1; ++i)
	{
		sum += X[i];
		cout << X[i] << ", ";
	}
	cout << X[N - 1] << "]" << endl;
	sum += X[N - 1];
	if (abs(sum - v(set1234)) < 0.000001)
	{
		cout << "Условие групповой рационализации выполняется" << endl;
	}
	else
	{
		cout << "Условие групповой рационализации не выполняется" << endl;
	}
	for (int i = 1; i <= N; ++i)
	{
		if (X[i - 1] < v({ i }))
		{
			cout << "Условие индивидуальной рационализации не выполняется" << endl;
			return;
		}
	}
	cout << "Условие индивидуальной рационализации выполняется" << endl;

}

int main()
{
	setlocale(LC_ALL, "Russian");
	check_superadditivity();
	check_convexity();
	find_Shapley();
	system("Pause");
}
