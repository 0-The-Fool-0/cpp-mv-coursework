#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>
#include <cstdlib> 
#include <ctime> 

using namespace std;

const int N = 15;

// Генерация матрицы A
vector<vector<double>> generateMatrixA() {
    vector<vector<double>> A(N, vector<double>(N));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            A[i][j] = (i == j) ? 5.0 * (i + 1) : 0.1 * (i + 1) * (j + 1);
    return A;
}


vector<double> generateVectorB() {
    srand(static_cast<unsigned int>(time(nullptr))); 
    vector<double> b(N);
    for (int i = 0; i < N; ++i)
        b[i] = (rand() % 10); 
    return b;
}

bool choleskyDecomposition(const vector<vector<double>>& A, vector<vector<double>>& L) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (int k = 0; k < j; ++k)
                sum += L[i][k] * L[j][k];
            if (i == j) {
                double val = A[i][i] - sum;
                if (val <= 0.0) return false; // Матрица не положительно определённая
                L[i][j] = sqrt(val);
            }
            else {
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }
    return true;
}


// Прямой и обратный ход для решения системы
vector<double> forwardSubstitution(const vector<vector<double>>& L, const vector<double>& b) {
    vector<double> y(N);
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j)
            sum += L[i][j] * y[j];
        y[i] = (b[i] - sum) / L[i][i];
    }
    return y;
}

vector<double> backwardSubstitution(const vector<vector<double>>& L, const vector<double>& y) {
    vector<double> x(N);
    for (int i = N - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < N; ++j)
            sum += L[j][i] * x[j];
        x[i] = (y[i] - sum) / L[i][i];
    }
    return x;
}



// Определитель через LL^T 
double determinantCholesky(const vector<vector<double>>& L) {
    double det = 1.0;
    for (int i = 0; i < N; ++i)
        det *= L[i][i];
    return det * det;
}


// QR-разложение
void qrDecomposition(const vector<vector<double>>& A, vector<vector<double>>& Q, vector<vector<double>>& R) {
    R = A;
    Q.assign(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i)
        Q[i][i] = 1.0;

    for (int j = 0; j < N - 1; ++j) {
        for (int i = N - 1; i > j; --i) {
            double a = R[i - 1][j];
            double b = R[i][j];
            double r = sqrt(a * a + b * b);
            if (r == 0) continue;
            double c = a / r;
            double s = -b / r;

            // Применяем вращение к R
            for (int k = 0; k < N; ++k) {
                double temp1 = R[i - 1][k];
                double temp2 = R[i][k];
                R[i - 1][k] = c * temp1 - s * temp2;
                R[i][k] = s * temp1 + c * temp2;
            }
            // Применяем вращение к Q
            for (int k = 0; k < N; ++k) {
                double temp1 = Q[k][i - 1];
                double temp2 = Q[k][i];
                Q[k][i - 1] = c * temp1 - s * temp2;
                Q[k][i] = s * temp1 + c * temp2;
            }
        }
    }
}
vector<double> solveQR(const vector<vector<double>>& Q, const vector<vector<double>>& R, const vector<double>& b) {
    // y = Q^T * b
    vector<double> y(N, 0.0);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            y[i] += Q[j][i] * b[j];

    // Решаем R*x = y (R - верхнетреугольная)
    vector<double> x(N);
    for (int i = N - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < N; ++j)
            sum += R[i][j] * x[j];
        x[i] = (y[i] - sum) / R[i][i];
    }
    return x;
}
double residualNorm(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b) {
    double norm = 0.0;
    for (int i = 0; i < N; ++i) {
        double Ax = 0.0;
        for (int j = 0; j < N; ++j)
            Ax += A[i][j] * x[j];
        norm += (Ax - b[i]) * (Ax - b[i]);
    }
    return sqrt(norm); // Евклидова норма
}

int main() {
    setlocale(LC_ALL, "Russian");
    auto A = generateMatrixA();
    auto b = generateVectorB();

    // Метод квадратного корня (Холецкого, LLᵗ)
    vector<vector<double>> L(N, vector<double>(N, 0.0));
    if (choleskyDecomposition(A, L)) {
        auto y = forwardSubstitution(L, b);
        auto x_cholesky = backwardSubstitution(L, y);
        double detA = determinantCholesky(L);
        cout << "Решение методом квадратного корня (Холецкого, LL^T):\n";
        for (double xi : x_cholesky) cout << xi << " ";
        cout << "\nОпределитель матрицы A: " << detA << endl;
        cout << "Норма невязки: " << residualNorm(A, x_cholesky, b) << endl;
    }
    else {
        cout << "Матрица не является положительно определённой.\n";
    }

    // QR-разложение
    vector<vector<double>> Q, R;
    qrDecomposition(A, Q, R);
    auto x_qr = solveQR(Q, R, b);
    cout << "\nРешение методом вращений (QR):\n";
    for (double xi : x_qr) cout << xi << " ";
    cout << "\nНорма невязки: " << residualNorm(A, x_qr, b) << endl;

    // Норма разности A - Q*R
    double normAQR = 0.0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double qr = 0.0;
            for (int k = 0; k < N; ++k)
                qr += Q[i][k] * R[k][j];
            normAQR += (A[i][j] - qr) * (A[i][j] - qr);
        }
    cout << "Норма ||A - QR||: " << sqrt(normAQR) << endl;

    return 0;
}


