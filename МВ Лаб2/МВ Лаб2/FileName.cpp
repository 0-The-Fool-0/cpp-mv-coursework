#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

const int n = 15;
const double epsilon = 1e-5;

// Функция для вычисления максимальной нормы вектора
double vector_norm(const vector<double>& v) {
    double max_val = 0.0;
    for (double val : v) {
        if (abs(val) > max_val) {
            max_val = abs(val);
        }
    }
    return max_val;
}

// Функция для вычисления нормы матрицы (максимальная сумма по строкам)
double matrix_norm(const vector<vector<double>>& A) {
    double max_sum = 0.0;
    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            row_sum += abs(A[i][j]);
        }
        if (row_sum > max_sum) {
            max_sum = row_sum;
        }
    }
    return max_sum;
}

// Функция для вычисления погрешности
double compute_error(const vector<double>& x, const vector<double>& x_true) {
    vector<double> diff(n);
    for (int i = 0; i < n; i++) {
        diff[i] = x[i] - x_true[i];
    }
    return vector_norm(diff);
}

int main() {
    setlocale(LC_ALL, "Rus");
    // Точное решение
    vector<double> x_true(n);
    for (int i = 0; i < n; i++) {
        x_true[i] = i + 2;  // x* = (2, 3, 4, ..., 16)
    }

    // Построение матрицы A
    vector<vector<double>> A(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                A[i][j] = 15 * pow(i + 1, 2);
            }
            else {
                A[i][j] = -1.0 / ((i + 1) * (j + 1));
            }
        }
    }

    // Построение вектора b = A * x_true
    vector<double> b(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i] += A[i][j] * x_true[j];
        }
    }

    // Начальное приближение
    vector<double> x0(n, 0.0);

    cout << "Точное решение: ";
    for (int i = 0; i < min(5, n); i++) {
        cout << x_true[i] << " ";
    }
    cout << "... " << x_true[n - 1] << endl << endl;

    // 1. Метод простой итерации
    cout << "1. МЕТОД ПРОСТОЙ ИТЕРАЦИИ" << endl;
    cout << "=========================" << endl;

    // Подбор параметра tau для сходимости
    double tau = 0.0006;
    vector<vector<double>> B(n, vector<double>(n, 0.0));
    vector<double> c(n, 0.0);

    // Построение матрицы B = E - tau * A
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                B[i][j] = 1 - tau * A[i][j];
            }
            else {
                B[i][j] = -tau * A[i][j];
            }
        }
        c[i] = tau * b[i];
    }
    double B_norm = matrix_norm(B);
    cout << "Норма матрицы B: " << B_norm << endl;

    vector<double> x_prev = x0;
    vector<double> x_curr(n, 0.0);
    int iterations_simple = 0;
    double error_estimate;

    do {
        iterations_simple++;

        // x_curr = B * x_prev + c
        for (int i = 0; i < n; i++) {
            x_curr[i] = c[i];
            for (int j = 0; j < n; j++) {
                x_curr[i] += B[i][j] * x_prev[j];
            }
        }

        // Вычисление оценки погрешности
        vector<double> diff(n);
        for (int i = 0; i < n; i++) {
            diff[i] = x_curr[i] - x_prev[i];
        }
        double diff_norm = vector_norm(diff);
        error_estimate = (B_norm / (1 - B_norm)) * diff_norm;

        x_prev = x_curr;

    } while (error_estimate > epsilon && iterations_simple < 10000);

    cout << "Количество итераций: " << iterations_simple << endl;
    cout << "Полученное решение: ";
    for (int i = 0; i < min(5, n); i++) {
        cout << fixed << setprecision(4) << x_curr[i] << " ";
    }
    cout << "... " << x_curr[n - 1] << endl;
    cout << "Погрешность: " << compute_error(x_curr, x_true) << endl << endl;

    // Априорная оценка количества итераций
    vector<double> x1(n, 0.0);
    for (int i = 0; i < n; i++) {
        x1[i] = c[i];
        for (int j = 0; j < n; j++) {
            x1[i] += B[i][j] * x0[j];
        }
    }

    vector<double> diff1(n);
    for (int i = 0; i < n; i++) {
        diff1[i] = x1[i] - x0[i];
    }
    double diff1_norm = vector_norm(diff1);

    double k_apriori = log(epsilon * (1 - B_norm) / diff1_norm) / log(B_norm);
    cout << "Априорная оценка итераций: " << ceil(k_apriori) << endl << endl;

    // 2. Метод Зейделя для системы x = Bx + c (альтернативная реализация)
    cout << "2. МЕТОД ЗЕЙДЕЛЯ ДЛЯ x = Bx + c" << endl;
    cout << "================================" << endl;

    vector<double> x_seidel_B = x0;
    int iterations_seidel_B = iterations_simple; // столько же, сколько в простой итерации

    // Создаем матрицы L_B и U_B
    vector<vector<double>> L_B(n, vector<double>(n, 0.0));
    vector<vector<double>> U_B(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j <= i) {
                L_B[i][j] = B[i][j];
            }
            else {
                U_B[i][j] = B[i][j];
            }
        }
    }

    for (int iter = 0; iter < iterations_seidel_B; iter++) {
        vector<double> x_old = x_seidel_B;
        for (int i = 0; i < n; i++) {
            double sum_L = 0.0;
            for (int j = 0; j < i; j++) {
                sum_L += L_B[i][j] * x_seidel_B[j];
            }
            double sum_U = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum_U += U_B[i][j] * x_old[j];
            }
            x_seidel_B[i] = (c[i] + sum_U + sum_L) / (1 - B[i][i]);
        }
    }

    cout << "Количество итераций: " << iterations_seidel_B << endl;
    cout << "Полученное решение: ";
    for (int i = 0; i < min(5, n); i++) {
        cout << fixed << setprecision(4) << x_seidel_B[i] << " ";
    }
    cout << "... " << x_seidel_B[n - 1] << endl;
    cout << "Погрешность: " << compute_error(x_seidel_B, x_true) << endl << endl;

    // 3. Метод релаксации
    cout << "3. МЕТОД РЕЛАКСАЦИИ" << endl;
    cout << "===================" << endl;

    // ω = 0.5
    vector<double> x_relax1 = x0;
    double omega1 = 0.5;

    for (int iter = 0; iter < iterations_simple; iter++) {
        for (int i = 0; i < n; i++) {
            double sum1 = 0.0;  // sum_{j=1}^{i-1} a_ij * x_j^{k+1}
            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * x_relax1[j];
            }

            double sum2 = 0.0;  // sum_{j=i+1}^{n} a_ij * x_j^k
            for (int j = i + 1; j < n; j++) {
                sum2 += A[i][j] * x_relax1[j];
            }

            x_relax1[i] = (1 - omega1) * x_relax1[i] +
                (omega1 / A[i][i]) * (b[i] - sum1 - sum2);
        }
    }

    cout << "ω = 0.5:" << endl;
    cout << "Количество итераций: " << iterations_simple << endl;
    cout << "Полученное решение: ";
    for (int i = 0; i < min(5, n); i++) {
        cout << fixed << setprecision(4) << x_relax1[i] << " ";
    }
    cout << "... " << x_relax1[n - 1] << endl;
    cout << "Погрешность: " << compute_error(x_relax1, x_true) << endl;

    // ω = 1.5
    vector<double> x_relax2 = x0;
    double omega2 = 1.5;

    for (int iter = 0; iter < iterations_simple; iter++) {
        for (int i = 0; i < n; i++) {
            double sum1 = 0.0;  // sum_{j=1}^{i-1} a_ij * x_j^{k+1}
            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * x_relax2[j];
            }

            double sum2 = 0.0;  // sum_{j=i+1}^{n} a_ij * x_j^k 
            for (int j = i + 1; j < n; j++) {
                sum2 += A[i][j] * x_relax2[j];
            }

            x_relax2[i] = (1 - omega2) * x_relax2[i] +
                (omega2 / A[i][i]) * (b[i] - sum1 - sum2);
        }
    }

    cout << "\nω = 1.5:" << endl;
    cout << "Количество итераций: " << iterations_simple << endl;
    cout << "Полученное решение: ";
    for (int i = 0; i < min(5, n); i++) {
        cout << fixed << setprecision(4) << x_relax2[i] << " ";
    }
    cout << "... " << x_relax2[n - 1] << endl;
    cout << "Погрешность: " << compute_error(x_relax2, x_true) << endl << endl;

    // 4. Сравнение методов
    cout << "4. СРАВНЕНИЕ МЕТОДОВ" << endl;
    cout << "====================" << endl;
    cout << "Метод простой итерации: " << compute_error(x_curr, x_true) << endl;
    cout << "Метод Зейделя: " << compute_error(x_seidel_B, x_true) << endl;
    cout << "Метод релаксации (ω=0.5): " << compute_error(x_relax1, x_true) << endl;
    cout << "Метод релаксации (ω=1.5): " << compute_error(x_relax2, x_true) << endl;


    return 0;
}