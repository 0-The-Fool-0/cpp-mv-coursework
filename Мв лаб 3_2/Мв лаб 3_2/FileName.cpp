#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace std;

const int N = 12;                
const double EPS = 1e-4;          
const int MAX_ITER = 2000;

using Matrix = vector<vector<double>>;
using Vector = vector<double>;

// Евклидова норма
double vecNorm(const Vector& v) {
    double s = 0.0;
    for (double x : v) s += x * x;
    return sqrt(s);
}

void normalize(Vector& v) {
    double n = vecNorm(v);
    if (n == 0) return; 
    for (double& x : v) x /= n;
}

Vector matVecMul(const Matrix& A, const Vector& x) {
    int n = (int)A.size();
    Vector y(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            y[i] += A[i][j] * x[j];
    return y;
}

Matrix matMul(const Matrix& A, const Matrix& B) {
    int n = (int)A.size();
    Matrix C(n, Vector(n, 0.0));
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < n; ++j)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

// --- Генерация симметричной матрицы (не диагональная и не трёхдиагональная) ---
Matrix generateMatrix(int n) {
    Matrix A = {
        {4, 1, 2, 0, 1, 3, 0, 2, 0, 0, 0, 1},
        {1, 5, 1, 2, 0, 1, 3, 0, 0, 2, 0, 0},
        {2, 1, 6, 0, 1, 2, 0, 3, 1, 0, 2, 0},
        {0, 2, 0, 7, 1, 0, 2, 1, 0, 3, 0, 2},
        {1, 0, 1, 1, 8, 2, 0, 0, 3, 0, 1, 0},
        {3, 1, 2, 0, 2, 9, 1, 0, 0, 2, 0, 1},
        {0, 3, 0, 2, 0, 1, 10, 2, 1, 0, 0, 3},
        {2, 0, 3, 1, 0, 0, 2, 11, 3, 1, 0, 0},
        {0, 0, 1, 0, 3, 0, 1, 3, 12, 2, 1, 0},
        {0, 2, 0, 3, 0, 2, 0, 1, 2, 13, 3, 1},
        {0, 0, 2, 0, 1, 0, 0, 0, 1, 3, 14, 2},
        {1, 0, 0, 2, 0, 1, 3, 0, 0, 1, 2, 15}
    };
    return A;
}   

void saveMatrix(const Matrix& A, const string& filename) {
    ofstream fout(filename);
    fout << fixed << setprecision(6);
    int n = (int)A.size();
    fout << n << "\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            fout << setw(12) << A[i][j] << " ";
        fout << "\n";
    }
    fout.close();
}

// --- Приведение к верхней Hessenberg (Householder) ---
void toHessenberg(Matrix& A) {
    int n = (int)A.size();
    for (int k = 0; k < n - 2; ++k) {
        // формируем вектор x (ниже главной поддиагонали в столбце k)
        Vector x(n - k - 1, 0.0);
        for (int i = 0; i < n - k - 1; ++i) x[i] = A[k + 1 + i][k];

        double sigma = vecNorm(x);
        if (sigma < 1e-14) continue;

        // вектор отражения v
        if (x[0] >= 0) sigma = -sigma;
        x[0] -= sigma;
        double vnorm = vecNorm(x);
        for (double& xi : x) xi /= vnorm;

        // применяем H = I - 2 v v^T к A слева: A = H*A
        for (int j = k; j < n; ++j) {
            double s = 0;
            for (int i = 0; i < (int)x.size(); ++i) s += x[i] * A[k + 1 + i][j];
            for (int i = 0; i < (int)x.size(); ++i) A[k + 1 + i][j] -= 2.0 * x[i] * s;
        }
        // и справа: A = A*H
        for (int i = 0; i < n; ++i) {
            double s = 0;
            for (int j = 0; j < (int)x.size(); ++j) s += x[j] * A[i][k + 1 + j];
            for (int j = 0; j < (int)x.size(); ++j) A[i][k + 1 + j] -= 2.0 * x[j] * s;
        }
    }
}

// --- QR-разложение методом отражений (возвращает Q и R: A = Q * R) ---
void qrDecomposition(const Matrix& A_in, Matrix& Q_out, Matrix& R_out) {
    int n = (int)A_in.size();
    Matrix R = A_in;
    Matrix Q(n, Vector(n, 0.0));
    for (int i = 0; i < n; ++i) Q[i][i] = 1.0;

    for (int k = 0; k < n - 1; ++k) {
        // строим вектор x = R[k..n-1][k]
        Vector x(n - k, 0.0);
        for (int i = 0; i < n - k; ++i) x[i] = R[k + i][k];

        double alpha = vecNorm(x);
        if (alpha < 1e-14) continue;
        if (x[0] >= 0) alpha = -alpha;
        x[0] -= alpha;
        double vnorm = vecNorm(x);
        for (double& xi : x) xi /= vnorm; // теперь x = v (unit)

        // сформировать H = I - 2 v v^T и применить к R слева: R = H * R
        for (int j = k; j < n; ++j) {
            double s = 0.0;
            for (int i = 0; i < n - k; ++i) s += x[i] * R[k + i][j];
            for (int i = 0; i < n - k; ++i) R[k + i][j] -= 2.0 * x[i] * s;
        }
        // применяем H справа к Q (чтобы Q = Q * H)
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int j = 0; j < n - k; ++j) s += Q[i][k + j] * x[j];
            for (int j = 0; j < n - k; ++j) Q[i][k + j] -= 2.0 * s * x[j];
        }
    }

    Q_out = Q;
    R_out = R;
}

// --- QR-алгоритм без сдвига (работаем с Hessenberg для ускорения) ---
vector<double> qrAlgorithm(Matrix A, int& iterations) {
    int n = (int)A.size();
    toHessenberg(A);

    iterations = 0;
    for (; iterations < MAX_ITER; ++iterations) {
        Matrix Q, R;
        qrDecomposition(A, Q, R);
        // A_{new} = R * Q
        A = matMul(R, Q);

        // критерий сходимости: сумма модулей элементов ниже главной диагонали (off-diagonal)
        double off = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j) off += fabs(A[i][j]);

        if (off < EPS) break;
    }

    vector<double> eig(n);
    for (int i = 0; i < n; ++i) eig[i] = A[i][i];
    return eig;
}

// --- Степенной метод (наибольшее по модулю собственное значение) ---
void powerMethod(const Matrix& A, double& lambda, Vector& eigenvector, int& iterations) {
    int n = (int)A.size();
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(-1.0, 1.0);

    Vector y(n);
    for (int i = 0; i < n; ++i) y[i] = dis(gen);
    normalize(y);

    iterations = 0;
    lambda = 0.0;

    while (iterations < MAX_ITER) {
        Vector Ay = matVecMul(A, y);
        // оценка собственного значения
        double lambda_next_num = inner_product(y.begin(), y.end(), Ay.begin(), 0.0);
        double lambda_next_den = inner_product(y.begin(), y.end(), y.begin(), 0.0);
        double lambda_next = lambda_next_num / lambda_next_den;

        // нормируем следующий вектор
        Vector y_next = Ay;
        normalize(y_next);

        // невязка 
        Vector r = matVecMul(A, y_next);
        for (int i = 0; i < n; ++i) r[i] -= lambda_next * y_next[i];
        double residual = vecNorm(r);

        // сохраняем результат если достигли точности
        if (residual < EPS) {
            lambda = lambda_next;
            eigenvector = y_next;
            ++iterations;
            break;
        }

        // подготовка к следующей итерации
        y = y_next;
        lambda = lambda_next;
        ++iterations;
    }
    // если не сошлось — вернуть текущее приближение
    if (eigenvector.empty()) {
        eigenvector = vector<double>(n, 0.0);
        Vector Ay = matVecMul(A, y);
        double lq = inner_product(y.begin(), y.end(), Ay.begin(), 0.0) / inner_product(y.begin(), y.end(), y.begin(), 0.0);
        lambda = lq;
        eigenvector = y;
    }
}
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    setlocale(LC_ALL, "Russian");
    // 1) Генерация матрицы
    Matrix A = generateMatrix(N);
    saveMatrix(A, "matrix.txt");
    cout << "Симметричная матрица (" << N << "x" << N << ") сохранена в matrix.txt\n";

    // 2) QR-алгоритм (без сдвига)
    int qr_iters = 0;
    vector<double> eigs = qrAlgorithm(A, qr_iters);

    cout << fixed << setprecision(8);
    cout << "\nQR-алгоритм (без сдвига):\n";
    cout << "Количество итераций: " << qr_iters << "\n";
    cout << "Собственные значения (диагонали после сходимости):\n";
    for (double val : eigs) cout << setw(14) << val << "\n";

    // 3) Степенной метод
    double lambda1;
    Vector x1;
    int p_iters = 0;
    powerMethod(A, lambda1, x1, p_iters);

    // вычислим невязку L2
    Vector Ax = matVecMul(A, x1);
    for (int i = 0; i < N; ++i) Ax[i] -= lambda1 * x1[i];
    double residual = vecNorm(Ax);

    cout << "\nСтепенной метод:\n";
    cout << "Количество итераций: " << p_iters << "\n";
    cout << "Наибольшее по модулю собственное значение (approx): " << lambda1 << "\n";
    cout << "Соответствующий собственный вектор (нормирован):\n";
    for (double v : x1) cout << setw(12) << v << " ";
    cout << "\n";
    cout << "||A*x - lambda*x||_2 = " << residual << "\n";

    return 0;
}
