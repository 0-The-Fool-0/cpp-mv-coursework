#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x) {
    return pow(x, 4) - x - 1;
}

double df(double x) {
    return 4 * pow(x, 3) - 1;
}

double phi(double x) {  
    return pow(x + 1, 0.25);
}

double dphi(double x) {
    return 0.25 * pow(x + 1, -0.75);
}

int main() {
    setlocale(LC_ALL, "Russian");

    const double eps = 1e-6;

    cout << fixed << setprecision(6);

    // ===============================
    // Метод простой итерации
    // ===============================

    double x0 = 1.0;
    double x1;
    int iter_simple = 0;

    double q = fabs(dphi(1));  
    cout << "Метод простой итерации:\n";

    do {
        x1 = phi(x0);
        iter_simple++;
        cout << "Итерация " << iter_simple
            << ": x = " << x1
            << ", |x_n - x_{n-1}| = " << fabs(x1 - x0)
            << endl;
        x0 = x1;
    } while ((q / (1 - q)) * fabs(x1 - phi(x1)) > eps);

    cout << "Корень (простая итерация): " << x1 << endl;
    cout << "Количество итераций: " << iter_simple << endl;


    // ===============================
    // Метод Ньютона
    // ===============================

    double xn = 1.2;
    double xn1;
    int iter_newton = 0;

    cout << "\nМетод Ньютона:\n";

    do {
        xn1 = xn - f(xn) / df(xn);
        iter_newton++;
        cout << "Итерация " << iter_newton
            << ": x = " << xn1
            << ", |x_n - x_{n-1}| = " << fabs(xn1 - xn)
            << endl;
        xn = xn1;
    } while (fabs(xn1 - (xn - f(xn) / df(xn))) > eps);

    cout << "Корень (Ньютон): " << xn << endl;
    cout << "Количество итераций: " << iter_newton << endl;

    return 0;
}