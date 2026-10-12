#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

const int n = 12;
const double EPS = 1e-4;
using Matrix = vector<vector<double>>;
using Vector = vector<double>;

// Генерация случайной симметричной матрицы (не диагональной и не трехдиагональной)
Matrix generateSymmetricMatrix(int size) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(-10, 10);

    Matrix A(size, Vector(size, 0.0));
    for (int i = 0; i < size; ++i) {
        for (int j = i; j < size; ++j) {
            if (i == j) {
                A[i][j] = dis(gen);
            }
            else {
                double val = dis(gen);
                A[i][j] = val;
                A[j][i] = val;
            }
        }
    }
    // Убедимся, что матрица не диагональная и не трехдиагональная
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (abs(i - j) > 1 && A[i][j] == 0) {
                A[i][j] = dis(gen);
                A[j][i] = A[i][j];
            }
        }
    }
    return A;
}

// Сохранение матрицы в файл
void saveMatrix(const Matrix& A, const string& filename) {
    ofstream fout(filename);
    fout << fixed << setprecision(6);
    for (const auto& row : A) {
        for (double val : row) {
            fout << setw(10) << val << " ";
        }
        fout << endl;
    }
    fout.close();
}

// Приведение к верхней форме Хессенберга (методом отражений Хаусхолдера)
void toHessenberg(Matrix& A) {
    for (int k = 0; k < n - 2; ++k) {
        Vector x(n - k - 1);
        for (int i = 0; i < n - k - 1; ++i)
            x[i] = A[k + 1 + i][k];
        double norm_x = sqrt(inner_product(x.begin(), x.end(), x.begin(), 0.0));
        if (norm_x < 1e-12) continue;
        if (x[0] >= 0) norm_x = -norm_x;
        x[0] -= norm_x;
        double norm_v = sqrt(inner_product(x.begin(), x.end(), x.begin(), 0.0));
        for (double& xi : x) xi /= norm_v;

        // Преобразование слева
        for (int i = k + 1; i < n; ++i) {
            double s = 0.0;
            for (int j = 0; j < n - k - 1; ++j)
                s += x[j] * A[k + 1 + j][i];
            for (int j = 0; j < n - k - 1; ++j)
                A[k + 1 + j][i] -= 2 * x[j] * s;
        }
        // Преобразование справа
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int j = 0; j < n - k - 1; ++j)
                s += x[j] * A[i][k + 1 + j];
            for (int j = 0; j < n - k - 1; ++j)
                A[i][k + 1 + j] -= 2 * x[j] * s;
        }
    }
}

// QR-разложение методом отражений Хаусхолдера
void qrDecomposition(const Matrix& A, Matrix& Q, Matrix& R) {
    int m = A.size();
    R = A;
    Q.assign(m, Vector(m, 0.0));
    for (int i = 0; i < m; ++i) Q[i][i] = 1.0;

    for (int k = 0; k < m - 1; ++k) {
        Vector x(m - k);
        for (int i = 0; i < m - k; ++i)
            x[i] = R[k + i][k];
        double norm_x = sqrt(inner_product(x.begin(), x.end(), x.begin(), 0.0));
        if (norm_x < 1e-12) continue;
        if (x[0] >= 0) norm_x = -norm_x;
        x[0] -= norm_x;
        double norm_v = sqrt(inner_product(x.begin(), x.end(), x.begin(), 0.0));
        for (double& xi : x) xi /= norm_v;

        // Преобразование R
        for (int j = k; j < m; ++j) {
            double s = 0.0;
            for (int i = 0; i < m - k; ++i)
                s += x[i] * R[k + i][j];
            for (int i = 0; i < m - k; ++i)
                R[k + i][j] -= 2 * x[i] * s;
        }
        // Преобразование Q
        for (int i = 0; i < m; ++i) {
            double s = 0.0;
            for (int j = 0; j < m - k; ++j)
                s += x[j] * Q[i][k + j];
            for (int j = 0; j < m - k; ++j)
                Q[i][k + j] -= 2 * x[j] * s;
        }
    }
    // Транспонируем Q
    for (int i = 0; i < m; ++i)
        for (int j = i + 1; j < m; ++j)
            swap(Q[i][j], Q[j][i]);
}

// QR-алгоритм (без сдвига)
vector<double> qrAlgorithm(Matrix A, int& iterations) {
    toHessenberg(A);
    int maxIter = 1000;
    int n = A.size();
    vector<double> eig(n);
    for (iterations = 0; iterations < maxIter; ++iterations) {
        Matrix Q, R;
        qrDecomposition(A, Q, R);
        // A = R * Q
        Matrix A1(n, Vector(n, 0.0));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k)
                    A1[i][j] += R[i][k] * Q[k][j];
        A = A1;
        // Проверка сходимости
        double off_diag = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j) off_diag += fabs(A[i][j]);
        if (off_diag < EPS) break;
    }
    for (int i = 0; i < n; ++i)
        eig[i] = A[i][i];
    return eig;
}

// Нормировка вектора
void normalize(Vector& v) {
    double norm = sqrt(inner_product(v.begin(), v.end(), v.begin(), 0.0));
    for (double& x : v) x /= norm;
}

// Умножение матрицы на вектор
Vector matVecMul(const Matrix& A, const Vector& x) {
    int n = A.size();
    Vector y(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            y[i] += A[i][j] * x[j];
    return y;
}

// Степенной метод
void powerMethod(const Matrix& A, double& lambda, Vector& eigenvector, int& iterations) {
    int n = A.size();
    Vector y(n, 1.0); // начальный вектор
    normalize(y);
    Vector y_next(n);
    lambda = 0.0;
    iterations = 0;
    while (iterations < 1000) {
        y_next = matVecMul(A, y);
        normalize(y_next);
        double lambda_next = inner_product(y_next.begin(), y_next.end(), matVecMul(A, y_next).begin(), 0.0);
        double diff = 0.0;
        Vector Ay = matVecMul(A, y_next);
        for (int i = 0; i < n; ++i)
            diff += fabs(Ay[i] - lambda_next * y_next[i]);
        if (diff < EPS) {
            lambda = lambda_next;
            eigenvector = y_next;
            break;
        }
        y = y_next;
        ++iterations;
    }
}

int main() {
    // 1. Генерация и сохранение матрицы
    Matrix A = generateSymmetricMatrix(n);
    saveMatrix(A, "matrix.txt");
    cout << "Симметричная матрица сохранена в matrix.txt\n";

    // 2. QR-алгоритм
    int qr_iters = 0;
    vector<double> eigs = qrAlgorithm(A, qr_iters);
    cout << "QR-алгоритм (без сдвига):\n";
    cout << "Количество итераций: " << qr_iters << endl;
    cout << "Собственные значения:\n";
    for (double v : eigs) cout << v << " ";
    cout << endl;

    // 3. Степенной метод
    double lambda1;
    Vector x1;
    int power_iters = 0;
    powerMethod(A, lambda1, x1, power_iters);
    cout << "\nСтепенной метод:\n";
    cout << "Количество итераций: " << power_iters << endl;
    cout << "Наибольшее по модулю собственное значение: " << lambda1 << endl;
    cout << "Собственный вектор:\n";
    for (double v : x1) cout << v << " ";
    cout << endl;

    // 4. Проверка невязки
    Vector Ax = matVecMul(A, x1);
    double residual = 0.0;
    for (int i = 0; i < n; ++i)
        residual += fabs(Ax[i] - lambda1 * x1[i]);
    cout << "Ax - λx = " << residual << endl;

    return 0;
}
