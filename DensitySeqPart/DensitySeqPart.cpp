#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <algorithm>

using namespace std;

double const PI = 3.1415;

// Решаемая задача
struct ro_problem
{
    string quantity = "Плотность";
    int method = 0;
    double D = 0.1, // Внешний диаметр трубы, м
        d = 0.008, // Толщина стенок трубы, м
        Q = 0.1,   // Объемный расход, м3/c
        maxT = 200, // Время моделирования
        L = 100, // Длина участка трубы
        xStep = 10,// Шаг сетки по координате
        dt, // Шаг сетки по времени
        ro_left = 840, // Граничное условие при (0, t)
        ro_right = 860; // Граничное условие при (L, t)
   vector<double> result; // Расчетный слой

};

// Вывод вектора слоев в файл (используется для построения графиков в Python)
// problem - указатель на структуру ro_problem
// filename - имя файла с результатами расчета
// формат вывода - "t, x, ro"
void print_layer_to_file(ro_problem& problem, int step, string filename = "data.txt")
{
    std::ofstream out;
    if (step == 0)
        out.open(filename);
    else
        out.open(filename, ios::app);
    if (out.is_open())
        for (size_t coord = 0; coord < problem.result.size(); coord++)
        {
            out << step * problem.dt << ',' << coord * problem.xStep << ',' << problem.result[coord] << endl;
        }
    out.close();
}

// Функция расчета задачи
// problem - указатель на структуру ro_problem
// verbose - флаг вывода результатов в консоль
void calc_ro_problem(ro_problem& problem, bool verbose = false)
{
    // Скорость потока
    double speed = abs(4.0 * problem.Q / (PI * (problem.D - problem.d) * (problem.D - problem.d)));

    bool isReverse = problem.Q < 0;


    // Шаг  во времени
    problem.dt = problem.xStep / speed;

    // Число узлов по координате и по времени
    int numStepsCoord = (problem.L / problem.xStep) + 1;
    int numStepsTime = (problem.maxT / problem.dt);

    if (verbose)
    {
        cout << "Рассчитываемая величина = " << problem.quantity << endl;
        cout << "Длина трубы = " << problem.L << endl;
        cout << "Время моделирования = " << problem.maxT << endl;
        cout << "Скорость потока = " << speed << endl;
        cout << "Шаг по координате = " << problem.xStep << endl;
        cout << "Шаг во времени = " << problem.dt << endl;
        cout << "Число узлов по времени = " << numStepsTime << endl;
        cout << "Число узлов по координате = " << numStepsCoord << endl;
    }

    // Заполнение начального слоя граничным условием
    problem.result = vector<double>(numStepsCoord, isReverse ? problem.ro_left : problem.ro_right);

    // Расчет слоев 
    for (size_t step = 0; step < problem.result.size() - 1; step++)
    {
        vector<double> temp = vector<double>(numStepsCoord);
        switch (problem.method)
        {
            // Расчет на основе свойства самоподобия
            case 0:
            {
                if (!isReverse)
                {
                    transform(
                        problem.result.begin(), // Текущая точка
                        problem.result.end() - 1,   // Предыдущая точка
                        temp.begin() + 1,
                        [](double a) { return a; }
                    );
                    temp[0] = problem.ro_left; // Дополнение слоя начальным условием
                }
                else
                {
                    transform(
                        problem.result.begin() + 1, // Текущая точка
                        problem.result.end(),   // Предыдущая точка
                        temp.begin(),
                        [](double a) { return a; }
                    );
                    problem.result.pop_back();
                    problem.result.push_back(problem.ro_right); // Дополнение слоя начальным условием
                }
                break;
            }
        // Расчет уголком
            case 1:
            {
                double sigma = problem.dt / problem.xStep;
                if (!isReverse)
                {
                    transform(
                        problem.result.begin(), // начало L1
                        problem.result.end() - 1,   // Вектор точек на одну назад; Размер совпадает с размером результата
                        problem.result.begin() + 1,
                        temp.begin() + 1,
                        [&sigma](double a, double b) { return (sigma * a) + ((1 - sigma) * b); }
                    );
                    temp[0] = problem.ro_left;
                }
                else
                {
                    transform(
                        problem.result.begin() + 1, // начало L1
                        problem.result.end(),   // Вектор точек на одну назад; Размер совпадает с размером результата
                        problem.result.begin(),
                        temp.begin(),
                        [&sigma](double a, double b) { return (sigma * a) + ((1.0 - sigma) * b); }
                    );
                    temp.pop_back();
                    temp.push_back(problem.ro_right);
                }
                break;
            }       
        }
        temp.swap(problem.result);

        if (verbose)
            for (size_t i = 0; i < problem.result.size(); i++)
                cout << printf("Шаг = %d    t = %.3f    x = %.3f    ro = %.3f", 
                    step,
                    step * problem.dt,
                    i * problem.xStep,
                    problem.result[i]) << endl;

        print_layer_to_file(problem, step);

    }
}




int main()
{
    setlocale(LC_ALL, "Russian"); // Установка кодировки для вывода кириллицы в консоль

    clock_t tStart = clock(); // Для замера времени работы

    // Создание структуры для решаемой задачи
    ro_problem problem = ro_problem();

    // Решение задачи
    calc_ro_problem(problem, true);

    printf("Время расчета: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    // Вызов скрипта построения графиков
    //system("py PlotPrecomputedGraph.py");
    return 0;
}