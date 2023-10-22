//реализация алгоритма Ланцоша - нахождения собственных значений симметричной матрицы
//Ход алгоритма Ланцоша может быть описан следующим образом:



#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace std;

class LanczosMethod
{
private:
    vector<vector<double>> matrix_;
public:
    LanczosMethod(vector<vector<double>> matrix) : matrix_(matrix) {
        ;
    };

    // Функция проверки на симметричность 
    bool checkSimmetry() const
    {
        if (matrix_.size() != matrix_[0].size()) return false;
        for (int i = 0; i < matrix_.size(); ++i)
            for (int j = i; j < matrix_.size(); ++j)
                if (matrix_[i][j] != matrix_[j][i]) return false;
        return true;
    }

    // Функция для умножения матрицы на вектор
    vector<double> matrixVectorMultiply(const vector<double>& vec) const {
        int n = matrix_.size();
        vector<double> result(n, 0.0);

        if (vec.size() != n)
        {
            cout << "\nIncorrect size of input vector in function matrixVectorMultiply, vector with null-values returned";
            return result;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i] += matrix_[i][j] * vec[j];
            }
        }

        return result;
    }

    // Функция для нормализации вектора
    void normalizeVector(vector<double>& vec) const {
        double norm = 0.0;
        for (double element : vec) {
            norm += element * element;
        }
        norm = sqrt(norm);

        for (double& element : vec) {
            element /= norm;
        }
    }

    // Функция для нахождения собственных значений методом Ланцоша
    vector<double> evaluateEigenvalues(int iterations) const {
        int n = matrix_.size();

        // Инициализация векторов eigenvalues и ortovalues
        vector<double> eigenvalues(iterations);
        vector<double> ortovalues(iterations - 1);

        if (!checkSimmetry())
        {
            cout << "\nInput matrix isn't simmetric, vector with null-values returned";
            return eigenvalues;
        }

        // Инициализация начального вектора
        vector<double> v(n, 1.0);
        normalizeVector(v);

        for (int i = 0; i < iterations; i++) {
            // Шаг Ланцоша
            vector<double> Av = matrixVectorMultiply(v);
            eigenvalues[i] = inner_product(v.begin(), v.end(), Av.begin(), 0.0);


            // Ортогонализация вектора Av из подпространства Крылова
            if (i > 0) {
                for (int j = 1; j < n; j++) {
                    Av[j] -= eigenvalues[i] * v[j] + ortovalues[i - 1] * v[j - 1];
                }
            }

            // Обновление вектора Крылова
            if (i < iterations - 1) {
                ortovalues[i] = sqrt(inner_product(Av.begin(), Av.end(), Av.begin(), 0.0));
                v = Av;
                normalizeVector(v);
            }
        }

        return eigenvalues;
    }
};

void printVector(const vector<double>& vec)
{
    for (int i = 0; i < vec.size(); ++i)
        cout << vec[i] << " ";
}

void printMatrix(const vector<vector<double>>& matrix)
{
    cout << "\n";
    for (int i = 0; i < matrix.size(); ++i)
    {
        for (int j = 0; j < matrix[i].size(); ++j)
            cout << matrix[i][j] << " ";
        cout << "\n";
    }
}

int main() {

    cout << "\nDemonstration of Lanczos algorithm: ";
    // Пример входной матрицы
    vector<vector<double>> matrix = {
        {2, -1, 0, -5},
        {-1, 2, -1, 3},
        {0, -1, 2,  4},
        {-5, 3, 4,  2}
    };

    vector<vector<double>> matrix1 = {
        {7, 8, 3},
        {8, 1, 0},
        {3, 0, 4}
    };

    // Задаем количество итераций
    int iterations = 3;

    // Запуск алгоритма Ланцоша
    cout << "\n1 example: ";
    printMatrix(matrix);
    LanczosMethod  test(matrix);
    vector<double> t;
    cout << "eigenvalues: ";
    t = test.evaluateEigenvalues(iterations + 1);
    printVector(t);

    cout << "\n\n\n2 example: ";
    printMatrix(matrix1);
    LanczosMethod test1(matrix1);
    cout << "eigenvalues: ";
    t = test1.evaluateEigenvalues(iterations);
    printVector(t);

    return 0;
}
