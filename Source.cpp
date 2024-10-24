#include <cmath>
#include <iostream>
#include <random>
#include <math.h>
#include <cstdlib>
#include <chrono>

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

unsigned long long factorial(const int& n) {
    unsigned long long result = 1;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}


double main_function(const double& x, const int& n, const double& mu, const double& lambda) {
    double value = exp(-abs(x - mu)) / (factorial(n - 1) * pow(2, n) * lambda);
    //cout << value << endl;
    double sum = 0;
    for (int i = 0; i <= n - 1; i++) {
        sum += (factorial(n - 1 + i) / (factorial(n - 1 - i) * factorial(i))) * (pow(abs(x - mu), n - 1 - i) / pow(2, i));
        //cout << "sum = " << sum << endl;
        //cout << "s " << (pow(abs(x - mu), n - 1 - i)) << endl;
    }
    //cout << sum << endl;
    value = value * sum;
    return value;
}

double calculateThirdMoment(const double& a, const double& b, const int& n_intervals,
    const int& n, const double& mu, const double& lambda, const double& mean) {
    double h = (b - a) / n_intervals;  // Шаг интегрирования
    double sum = 0.5 * (pow(a - mean, 3) * main_function(a, n, mu, lambda) + pow(b - mean, 3) * main_function(b, n, mu, lambda));  // Начальные значения

    for (int i = 1; i < n_intervals; ++i) {
        double x = a + i * h;
        sum += pow(x - mean, 3) * main_function(x, n, mu, lambda);  // (x - mean)^3 * f(x)
    }

    return sum * h;
}

double calculateMean(const double& a, const double& b, const int& n_intervals,
    const int& n, const double& mu, const double& lambda) {
    double h = (b - a) / n_intervals;  // Шаг интегрирования
    double sum = 0.5 * (a * main_function(a, n, mu, lambda) + b * main_function(b, n, mu, lambda));  // Начальные значения

    for (int i = 1; i < n_intervals; ++i) {
        double x = a + i * h;
        sum += x * main_function(x, n, mu, lambda);  // x * f(x)
    }

    return sum * h;
}



double* main_metrics(const int& n, const double& mu, const double& lambda) {
    double n1 = n;
    double* result = new double[4];
    result[1] = 2 * n1; // дисперсия
    result[3] = 3 / n1; // эксцесс
    double a = -25;
    double b = 25;
    int n_intervals = 10000;

    double mean = calculateMean(a, b, n_intervals, n, mu, lambda); // мат ожидание
    result[0] = mean;
    double thirdMoment = calculateThirdMoment(a, b, n_intervals, n, mu, lambda, mean);
    double stddev = sqrt(result[1]);
    result[2] = thirdMoment / pow(stddev, 3); // коэф. ассиметрии

    return result;
}

double pseudo_number() {
    double random;
    do random = (double)rand() / RAND_MAX;
    while (random == 0. || random == 1.);
    return random;
}

double Laplace_random_value(const int& n, const double& mu, const double& lambda) {
    double p1 = 1;
    double p2 = 1;
    for (int i = 1; i <= n; i++) {
        double rand = pseudo_number();
        if (rand <= 0.5) {
            p1 = p1 * 2 * rand;
        }
        else {
            p2 = p2 * 2 * (1 - rand);
        }
    }
    return log(p1 / p2) * lambda + mu;
}

double mixture(const double& x, const double& p,
    const int& n1, const double& mu1, const double& lambda1,
    const int& n2, const double& mu2, const double& lambda2) {
    return p * main_function(x, n1, mu1, lambda1) + (1.0 - p) * main_function(x, n2, mu2, lambda2);
}

double* mixture_metrics(const double& p, const int& n1, const double& mu1, const double& lambda1,
    const int& n2, const double& mu2, const double& lambda2) {
    double* result = new double[4];
    double* mixture1 = new double[4];
    double* mixture2 = new double[4];
    mixture1 = main_metrics(n1, mu1, lambda1);
    mixture2 = main_metrics(n2, mu2, lambda2);

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
    const int& n1, const int& mu1, const double& lambda1,
    const double& n2, const double& mu2, const double& lambda2) {
    double p1 = pseudo_number();

    if (p1 > p) {
        return Laplace_random_value(n1, mu1, lambda1);
    }
    else {
        return Laplace_random_value(n2, mu2, lambda2);
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
    double* result = new double[dataSize];
    sort(data, data + dataSize);
    int index;
    for (int i = 0; i < dataSize; i++) {
        index = static_cast<int> (pseudo_number() * dataSize);
        if (index >= dataSize) index = dataSize - 1;
        result[i] = data[index] + create_noise();
    }
    return result;
}

void start_main_distribution() {
    bool flagIsWorking = true; //Состояние программы
    double x, mu, lambda;
    double* metrics = new double[4];
    int select;
    int n;
    cout << "Введите параметр n > 0 для распределения:";
    n = input_value_int(1, 100000);
    cout << "Введите параметр mu для распределения:";
    mu = input_value(-100000, 100000);
    cout << "Введите параметр lambda для распределения:";
    lambda = input_value(0, 100000);
    cout << endl;

    while (flagIsWorking) {
        print_actions();
        cin >> select;
        switch (select) {
        case 1:
            cout << "Введите x0:";
            x = input_value(-100000, 100000);
            cout << "Плотность распределения Лапласа в точке x0 = " << x << " :  " <<
                main_function(x, n, mu, lambda) << endl;
            break;
        case 2:
            metrics = main_metrics(n, mu, lambda);
            cout << "Математическое ожидание: " << metrics[0] << endl;
            cout << "Дисперсия: " << metrics[1] << endl;
            cout << "Коэффициент асимметрии: " << metrics[2] << endl;
            cout << "Коэффициент эксцесса: " << metrics[3] << endl;
            cout << endl;
            break;
        case 3:
            cout << "Смоделированная случайная величина : " << Laplace_random_value(n, mu, lambda) << endl;
            break;
        case 4:
            flagIsWorking = false;
            break;
        }
        if (select < 1 or select > 4) {
            cout << "Ошибка: номер дейтсвия введен неверно. Повторите попытку." << endl;
        }
    }
}

void start_mixture_distribution() {
    double x, mu1, lambda1, p;
    double mu2, lambda2;
    int n1, n2, select;
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

    while (flagIsWorking) {
        print_actions();
        cin >> select;
        switch (select) {
        case 1:
            cout << "Введите x0:";
            x = input_value(-100000, 100000);
            cout << "Плотность смеси распределений Лапласа в точке x0 = " << x <<
                " :  " << mixture(x, p, n1, mu1, lambda1, n2, mu2, lambda2) << endl;
            break;
        case 2:
            metrics = mixture_metrics(p, n1, mu1, lambda1, n2, mu2, lambda2);
            cout << "Математическое ожидание: " << metrics[0] << endl;
            cout << "Дисперсия: " << metrics[1] << endl;
            cout << "Коэффициент асимметрии: " << metrics[2] << endl;
            cout << "Коэффициент эксцесса: " << metrics[3] << endl;
            cout << endl;
            break;
        case 3:
            cout << "Смоделированная случайная величина : " <<
                mixture_random_value(p, n1, mu1, lambda1, n2, mu2, lambda2) << endl;
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


    cout << "Введите способ создания выборки: " << endl;
    cout << "1 - Генерация выборки на основе основного распределения" << endl;
    cout << "2 - Генерация выборки на основе смеси основных распределений" << endl;
    cout << "3 - Ввод выборки с клавиатуры" << endl;
    cin >> choice;
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
            emp_sample[i] = Laplace_random_value(n, mu, lambda);

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
        for (int i = 0; i < dataSize; i++) {
            emp_sample[i] = mixture_random_value(p, n1, mu1, lambda1, n2, mu2, lambda2);
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
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
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
    //get_emp_3_3_2(data, dataSize);
    while (flagIsWorking) {
        print_actions();
        cout << "5 - Пересоздание выборки основываясь на текущей" << endl;
        cout << "6 - Построение интервалов и плотностей для графиков на основе распределения Лапласа" << endl;
        cout << "7 - Построение интервалов и плотностей для графиков на основе смесей" << endl;
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
            start = std::chrono::system_clock::now();

            for (double i = borders_density[0][0]; i <= borders_density[1][k - 1]; i += 0.1) {
                x = i;
                cout << main_function(x, n, mu, lambda) << endl;
            }
            end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end - start;
            

            std::cout << "finished computation at " << std::ctime(&end_time)
                << "elapsed time: " << elapsed_seconds.count() << "s\n";
            cout << "time " << elapsed_seconds.count() << endl;
            break;
        case 7:
            borders_density = intervals_and_densities(data, dataSize);
            k = log2(dataSize) + 1;
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

            cout << endl;
            for (int j = 0; j < k; j++) {
                cout << "Плотность на интервале [ " << borders_density[0][j] << ", " << borders_density[1][j] << " ] = " << borders_density[2][j] << endl;
            }
            cout << "Начало для подсчета плотностей: " << borders_density[0][0] << endl;
            for (double i = borders_density[0][0]; i <= borders_density[1][k - 1]; i += 0.1) {
                x = i;
                cout << mixture(x, p, n1, mu1, lambda1, n2, mu2, lambda2) << endl;
            }
            break;

        }

        if (select < 1 or select > 7) {
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
    while (IsWorkingMain) {
        cout << "Выберите желаемое распределение:" << endl;
        cout << "1 - Распределение Лапласа" << endl;
        cout << "2 - Смесь двух нормальных распределений Лапласа" << endl;
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


