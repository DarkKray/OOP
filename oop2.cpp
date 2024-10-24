#include <cmath>
#include <iostream>
#include <random>
#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <chrono>

#include "Primary.h"
#include "Test_Primaryh.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

double input_value(const double& bottom, double top) {
    double number;
    if (top < bottom) {
        top = numeric_limits<float>::max();
    }
    do {
        cin >> number;
        if (cin.fail()) {
            cin.clear(); // Сбросить ошибку ввода
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // Очистить буфер
            cout << "Ошибка: введеный символ не является числом. Попробуйте снова." << endl;
        }
        else if (number < bottom || number > top) {
            cout << "Ошибка: число находится вне диапазона. Попробуйте снова." << endl;
        }
    } while (number < bottom || number > top);
    return number;
}
double input_value_int(const double& bottom, double top) {
    int number;
    if (top < bottom) {
        top = numeric_limits<float>::max();
    }
    do {
        cin >> number;
        if (cin.fail()) {
            cin.clear(); // Сбросить ошибку ввода
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // Очистить буфер
            cout << "Ошибка: введеный символ не является числом. Попробуйте снова." << endl;
        }
        else if (number < bottom || number > top) {
            cout << "Ошибка: число находится вне диапазона. Попробуйте снова." << endl;
        }
    } while (number < bottom || number > top);
    return number;
}


void print_actions() {
    cout << "Выберите действие: " << endl;
    cout << "1 - Расчет плотности распределения" << endl;
    cout << "2 - Вычисление математического ожидания, дисперсии и коэффициентов асимметрии и эксцесса" << endl;
    cout << "3 - Моделирование случайной величины из распределения" << endl;
    cout << "4 - Переход в главное меню" << endl;

    cout << endl;
    return;
}

double mixture(const double& x, const double& p, Primary laplace1, Primary laplace2) {
    return p * laplace1.Laplace_function(x) + (1.0 - p) * laplace2.Laplace_function(x);
}

double* mixture_metrics(const double& p, Primary laplace1, Primary laplace2) {
    double* result = new double[4];
    double* mixture1 = new double[4];
    double* mixture2 = new double[4];
    mixture1 = laplace1.get_metrics();
    mixture2 = laplace2.get_metrics();

    result[0] = p * mixture1[0] + (1 - p) * mixture2[0];

    result[1] = (p * (pow(mixture1[0], 2) + mixture1[1])) +
        (1 - p) * (pow(mixture2[0], 2) + mixture2[1]) - pow(result[0], 2);

    result[2] = (1 / pow(result[1], 1.5)) *
        (p * (pow(mixture1[0] - result[0], 3) +
            3 * (mixture1[0] - result[0]) * mixture1[1] +
            pow(mixture1[1], 1.5) * mixture1[2]) +
            (1 - p) * (pow(mixture2[0] - result[0], 3) +
                3 * (mixture2[0] - result[0]) * mixture2[1] +
                pow(mixture2[1], 1.5) * mixture2[2]));

    result[3] = (p * (pow(mixture1[0] - result[0], 4) +
        6 * pow(mixture1[0] - result[0], 2) * mixture1[1] +
        4 * (mixture1[0] - result[0]) * pow(mixture1[1], 3 / 2) * mixture1[2] +
        mixture1[1] * mixture1[1] * (mixture1[3] + 3)) +
        (1 - p) * (pow(mixture2[0] - result[0], 4) +
            6 * pow(mixture2[0] - result[0], 2) * mixture2[1] +
            4 * (mixture2[0] - result[0]) * pow(mixture2[1], 3 / 2) * mixture2[2] +
            mixture2[1] * mixture2[1] * (mixture2[3] + 3))) / pow(result[1], 2);
    return result;
}

double mixture_random_value(const double& p,
    Primary laplace1, Primary laplace2) {
    double p1 = laplace1.pseudo_number();

    if (p1 > p) {
        return laplace1.Laplace_random_value();
    }
    else {
        return laplace2.Laplace_random_value();
    }
}

double* empirical_metrics(const double* data, int dataSize) {
    double* result = new double[4];
    double sum = 0.0;

    for (int i = 0; i < dataSize; i++) { // мат. ожидание
        sum += data[i];
    }
    result[0] = sum / dataSize;
    sum = 0.0;

    for (int i = 0; i < dataSize; i++) {
        sum += pow(data[i] - result[0], 2);// дисперсия
    }
    result[1] = sum / dataSize;
    sum = 0.0;

    for (int i = 0; i < dataSize; i++) {
        sum += pow(data[i] - result[0], 3);
    }
    result[2] = sum / (dataSize * pow(result[1], 1.5)); // коэф. ассиметрии
    sum = 0.0;

    for (int i = 0; i < dataSize; i++) {
        sum += pow(data[i] - result[0], 4);
    }
    result[3] = sum / (dataSize * pow(result[1], 2)) - 3; // коэф. эксцесса

    return result;
}

double* shift_scaling(const double* data, double shift, double scale, int dataSize) {
    double* result = new double[dataSize];
    for (int i = 0; i < dataSize; i++) {
        result[i] = data[i] / scale - shift;
    }
    return result;
}

double empirical_random_value(const double* data, const int& dataSize) {
    double scale = 10;
    double curr = -1000000000000;

    // Поиск максимального значения в массиве
    for (int i = 0; i < dataSize; i++) {
        if (data[i] > curr) {
            curr = data[i];
        }
    }

    double shift = curr / scale + 1;
    double* TransformedData = new double[dataSize];

    // Масштабирование и сдвиг данных
    TransformedData = shift_scaling(data, shift, scale, dataSize);

    // Генерация случайного индекса
    static std::mt19937 generator(static_cast<unsigned int>(std::time(0)));
    std::uniform_real_distribution<> dis(0, dataSize - 1);
    int index = dis(generator);

    return (TransformedData[index] + shift) * 10;
}

double** intervals_and_densities(const double* info, const int& dataSize) {
    double** result = new double* [3];
    int k = log2(dataSize) + 1;
    double* density = new double[k];
    double* lower_borders = new double[k];
    double* higher_borders = new double[k];
    double mx = -1e15;
    double mn = 1e15;
    for (int i = 0; i < dataSize; i++) {
        mx = max(info[i], mx);
        mn = min(info[i], mn);
    }
    mx = static_cast<int>(ceil(mx));
    mn = static_cast<int>(ceil(mn));
    double step = fabs(mx - mn) / k;
    int counter = 0;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < dataSize; j++) {
            if (info[j] <= (mn + step * (i + 1)) && info[j] > (mn + step * i)) {
                counter++;
            }
        }
        lower_borders[i] = mn + step * i;
        higher_borders[i] = mn + step * (i + 1);
        density[i] = (double)counter / dataSize;
        counter = 0;
    }
    result[0] = lower_borders;
    result[1] = higher_borders;
    result[2] = density;
    return result;
}

double empirical_distribution(const int& x, const double* data, const int& dataSize) {
    double result;
    double** borders_density = new double* [3];
    borders_density = intervals_and_densities(data, dataSize);
    int k = log2(dataSize) + 1;
    for (int i = 0; i < k; i++) {
        if (borders_density[0][i] < x && borders_density[1][i] >= x) {
            result = borders_density[2][i];
            return result;
        }
    }
    return 0;
}
double create_noise() {
    return -0.001 + static_cast<double>(rand()) / RAND_MAX * (0.002);
}

double* creating_new_emp_data(double* data, const int& dataSize) {
    Primary Laplace(1, 0, 1.0);
    double* result = new double[dataSize];
    sort(data, data + dataSize);
    int index;
    for (int i = 0; i < dataSize; i++) {
        index = static_cast<int> (Laplace.pseudo_number() * dataSize);
        if (index >= dataSize) index = dataSize - 1;
        result[i] = data[index] + create_noise();
    }
    return result;
}

void start_main_distribution() {
    Primary Laplace(1, 0, 1.0);

    bool flagIsWorking = true; //Состояние программы
    double x, mu, lambda;
    double* metrics = new double[4];
    int select, select2;
    int n;
    cout << "Введите способ ввода параметров. " << endl;
    cout << "1. Ввод вручную. " << endl;
    cout << "2. Ввод из файла. " << endl;
    cin >> select2;
    switch (select2) {
    case 1:
        cout << "Введите параметр n > 0 для распределения:";
        n = input_value_int(0, 100000);
        cout << "Введите параметр mu для распределения:";
        mu = input_value(-100000, 100000);
        cout << "Введите параметр lambda для распределения:";
        lambda = input_value(0, 100000);
        cout << endl;
        Laplace.set_n(n);
        Laplace.set_mu(mu);
        Laplace.set_lambda(lambda);
        break;
    case 2:
        cout << "Введите путь к файлу с данными: " << endl;
        string source;
        cin >> source;
        Laplace.load(source);
        break;
    }
    if (select2 < 1 or select2 > 2) {
        cout << "Ошибка: номер дейтсвия введен неверно. Повторите попытку." << endl;
    }
    while (flagIsWorking) {
        string output;
        print_actions();
        cout << "5 - Вывод в файл" << endl;
        cin >> select;
        switch (select) {
        case 1:
            cout << "Введите x0:";
            x = input_value(-100000, 100000);
            cout << "Плотность распределения Лапласа в точке x0 = " << x << " :  " <<
                Laplace.Laplace_function(x) << endl;
            break;
        case 2:
            metrics = Laplace.get_metrics();
            cout << "Математическое ожидание: " << metrics[0] << endl;
            cout << "Дисперсия: " << metrics[1] << endl;
            cout << "Коэффициент асимметрии: " << metrics[2] << endl;
            cout << "Коэффициент эксцесса: " << metrics[3] << endl;
            cout << endl;
            break;
        case 3:
            cout << "Смоделированная случайная величина : " << Laplace.Laplace_random_value() << endl;
            break;
        case 5:
            cout << "Введите путь к файлу для записи " << endl;
            cin >> output;
            Laplace.save(output);
            cout << "Запись успешно завершена. " << endl;
            break;
        case 4:
            flagIsWorking = false;
            break;
        }
        if (select < 1 or select > 5) {
            cout << "Ошибка: номер дейтсвия введен неверно. Повторите попытку." << endl;
        }
    }
}

void start_mixture_distribution() {
    Primary Laplace(1, 0.0, 1.1);

    double x, mu1, lambda1, p;
    double mu2, lambda2;
    int n1, n2, select, select2;
    bool flagIsWorking = true; //Состояние программы
    double* metrics = new double[4];

    cout << "Введите параметр n > 0 для первого распределения:";
    n1 = input_value_int(0, 100000);
    cout << "Введите параметр mu для первого распределения:";
    mu1 = input_value(-100000, 100000);
    cout << "Введите параметр lambda для первого распределения:";
    lambda1 = input_value(0, 100000);
    cout << endl;
    cout << "Введите параметр n > 0 для второго распределения:";
    n2 = input_value_int(0, 100000);
    cout << "Введите параметр mu для второго распределения:";
    mu2 = input_value(-100000, 100000);
    cout << "Введите параметр lambda для второго распределения:";
    lambda2 = input_value(0, 100000);
    cout << "Введите параметр смеси p в интервале [0, 1]:";
    p = input_value(0, 1);
    cout << endl;
    Primary Laplace_1(n1, mu1, lambda1);
    Primary Laplace_2(n2, mu2, lambda2);

    while (flagIsWorking) {
        print_actions();
        cin >> select;
        switch (select) {
        case 1:
            cout << "Введите x0:";
            x = input_value(-100000, 100000);
            cout << "Плотность смеси распределений Лапласа в точке x0 = " << x <<
                " :  " << mixture(x, p, Laplace_1, Laplace_2) << endl;
            break;
        case 2:
            metrics = mixture_metrics(p, Laplace_1, Laplace_2);
            cout << "Математическое ожидание: " << metrics[0] << endl;
            cout << "Дисперсия: " << metrics[1] << endl;
            cout << "Коэффициент асимметрии: " << metrics[2] << endl;
            cout << "Коэффициент эксцесса: " << metrics[3] << endl;
            cout << endl;
            break;
        case 3:
            cout << "Смоделированная случайная величина : " <<
                mixture_random_value(p, Laplace_1, Laplace_2) << endl;
            break;
        case 4:
            flagIsWorking = false;
        }
        if (select < 1 or select > 4) {
            cout << "Ошибка: номер дейтсвия введен неверно. Повторите попытку." << endl;
        }
    }
}

double* creating_emp_sample(const int& dataSize) {

    double x, mu1, lambda1, p;
    double mu2, lambda2, mu, lambda;
    int n1, n2, select, choice, n;
    bool flagIsWorking = true; //Состояние программы
    double* metrics = new double[4];
    Primary Laplace(1, 0.0, 1.1);
    Primary Laplace_1(1, 0, 1.0);
    Primary Laplace_2(1, 0, 1.0);
    cout << "Введите способ создания выборки: " << endl;
    cout << "1 - Генерация выборки на основе основного распределения" << endl;
    cout << "2 - Генерация выборки на основе смеси основных распределений" << endl;
    cout << "3 - Ввод выборки с клавиатуры" << endl;
    cin >> choice;
    cout << "Введите количество элементов в выборке: ";
    if (choice < 1 or choice > 4) {
        cout << "Ошибка: номер дейтсвия введен неверно. Повторите попытку." << endl;
        return NULL;
    }
    double* emp_sample = new double[dataSize];
    switch (choice) {
    case 1:
        cout << "Введите параметр n > 0 для распределения:";
        n = input_value_int(0, 100000);
        cout << "Введите параметр mu для распределения:";
        mu = input_value(-100000, 100000);
        cout << "Введите параметр lambda для распределения:";
        lambda = input_value(0, 100000);
        cout << endl;

        for (int i = 0; i < dataSize; i++) {
            emp_sample[i] = Laplace.Laplace_random_value();
        }
        break;
    case 2:
        cout << "Введите параметр n1 > 0 для первого распределения:";
        n1 = input_value_int(0, 100000);
        cout << "Введите параметр mu1 для первого распределения:";
        mu1 = input_value(-100000, 100000);
        cout << "Введите параметр lambda1 для первого распределения:";
        lambda1 = input_value(0, 100000);
        cout << endl;
        cout << "Введите параметр n1 > 0 для второго распределения:";
        n2 = input_value_int(0, 100000);
        cout << "Введите параметр mu2 для второго распределения:";
        mu2 = input_value(-100000, 100000);
        cout << "Введите параметр lambda3 для второго распределения:";
        lambda2 = input_value(0, 100000);
        cout << "Введите параметр смеси p в интервале [0, 1]:";
        p = input_value(0, 1);
        cout << endl;
        Laplace_1.set_n(n1);
        Laplace_1.set_mu(mu1);
        Laplace_1.set_lambda(lambda1);
        Laplace_2.set_n(n2);
        Laplace_2.set_mu(mu2);
        Laplace_2.set_lambda(lambda2);
        for (int i = 0; i < dataSize; i++) {
            emp_sample[i] = mixture_random_value(p, Laplace_1, Laplace_2);
        }
        break;
    case 3:
        for (int i = 0; i < dataSize; i++) {
            cout << "Введите " << i + 1 << " элемент: ";
            cin >> emp_sample[i];
        }
    }
    return emp_sample;
}


void start_empirical_distribution() {
    Primary laplace(1, 0.0, 1.0);
    bool flagIsWorking = true; //Состояние программы
    int dataSize;
    double* metrics = new double[4];
    cout << "Введите количество элементов в выборке: ";
    dataSize = input_value_int(2, 10000000);
    double* data = creating_emp_sample(dataSize);
    int select, n, n1, n2;
    double x, mu, lambda, lambda1, lambda2, mu1, mu2, p;
    double** borders_density;
    int k;
    if (data == NULL) {
        return;
    };
    while (flagIsWorking) {
        print_actions();
        cout << "5 - Пересоздание выборки основываясь на текущей" << endl;
        cout << "6 - Моделирование теоритической выбороки" << endl;
        cin >> select;
        switch (select) {
        case 1:
            cout << "Введите x0:";
            x = input_value(-100000, 100000);
            cout << "Плотность эмпирического распределения в точке x0 = " << x <<
                " :  " << empirical_distribution(x, data, dataSize) << endl;
            break;
        case 2:
            metrics = empirical_metrics(data, dataSize);
            cout << "Математическое ожидание: " << metrics[0] << endl;
            cout << "Дисперсия: " << metrics[1] << endl;
            cout << "Коэффициент асимметрии: " << metrics[2] << endl;
            cout << "Коэффициент эксцесса: " << metrics[3] << endl;
            cout << endl;
            break;
        case 3:
            cout << "Смоделированная случайная величина : " << empirical_random_value(data, dataSize) << endl;
            break;
        case 4:
            flagIsWorking = false;
            break;
        case 5:
            data = creating_new_emp_data(data, dataSize);
            break;
        case 6: 
            borders_density = intervals_and_densities(data, dataSize);
            k = log2(dataSize) + 1;
            cout << "Введите параметр n > 0 для распределения:";
            n = input_value_int(0, 100000);
            cout << "Введите параметр mu для распределения:";
            mu = input_value(-100000, 100000);
            cout << "Введите параметр lambda для распределения:";
            lambda = input_value(0, 100000);
            cout << endl;

            for (int j = 0; j < k; j++) {
                cout << "( " << borders_density[0][j] << ", " << borders_density[1][j] << " ) = " << borders_density[2][j] << endl;
            }
            const auto start{ chrono::steady_clock::now() };

            for (double i = borders_density[0][0]; i <= borders_density[1][k - 1]; i += 0.1) {
                x = i;
                laplace.set_n(n);
                laplace.set_mu(mu);
                laplace.set_lambda(lambda);
                cout << laplace.Laplace_function(x) << endl;
            }
            const auto end{ chrono::steady_clock::now() };
            const chrono::duration<double> elapsed_seconds{ end - start };
            cout << "time " << elapsed_seconds.count() << endl;
            break;
        }

        if (select < 1 or select > 6) {
            cout << "Ошибка: номер дейтсвия введен неверно. Повторите попытку." << endl;
        }
    }
}

int main()
{
    setlocale(LC_ALL, "Russian");
    double* parametrs; //Дисперсия, Мат ожидание, Асимметрия и Эксцесс
    double* inputed_parametrs; //Входные данные
    double* data; //Выборка для эмпирического
    double x; //Икс для плотностей
    bool flagIsWorking = true; //Состояние программы
    int select; //Выбранное действие в программе

    bool IsWorkingMain = true;
    Primary Laplace(1, 0.0, 1.0);
    if (testPrimary(Laplace)) {
        cout << "Все Unit тесты выполнены. " << endl;
        while (IsWorkingMain) {
            cout << "Выберите желаемое распределение:" << endl;
            cout << "1 - Нормальное распределение Гаусса" << endl;
            cout << "2 - Смесь двух нормальных распределений Гаусса" << endl;
            cout << "3 - Эмпирическое распределение построенное по выборке" << endl;
            cout << "4 - Завершение работы программы" << endl;
            cin >> select;
            flagIsWorking = true;

            switch (select) {
            case 1:
                start_main_distribution();
                break;
            case 2:
                start_mixture_distribution();
                break;
            case 3:
                start_empirical_distribution();
                break;
            case 4:
                IsWorkingMain = false;
                break;
            }
            if (select < 1 or select > 5) {
                cout << "Ошибка: номер дейтсвия введен неверно. Повторите попытку." << endl;
            }
        }
    }

}


