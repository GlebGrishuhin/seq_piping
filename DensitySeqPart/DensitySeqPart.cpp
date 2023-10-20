#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <algorithm>

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

using namespace std;

typedef composite_layer_t<profile_collection_t<1>, moc_solver<1>::specific_layer> single_var_moc_t;

double const PI = 3.1415;


// Решаемая задача
struct ro_problem
{
    string quantity = "Плотность";
    int method = 0;
    double D = 0.1; // Внешний диаметр трубы, м
    double d = 0.008; // Толщина стенок трубы; м
    double Q = 0.01;   // Объемный расход; м3/c
    double maxT = 3600; // Время моделирования
    double L = 10000; // Длина участка трубы
    double xStep = 100;// Шаг сетки по координате
    double dt = xStep / abs(4.0 * Q / (PI * (D - d) * (D - d))); // Шаг сетки по времени
    double ro_left = 860; // Граничное условие при (0; t)
    double ro_right = 860; // Граничное условие при (L; t)
    double ro_init = 850;

    void set_initial_layer(custom_buffer_t<single_var_moc_t>& buffer)
    {
        buffer.current().vars.point_double[0] = vector<double>((L / xStep) + 1, ro_init);
    }
};


/// @brief Расчет следующего слоя по свойству самоподобия
/// @param problem Ссылка на структуру ro_problem
/// @param curr Ссылка на текущий слой в буфере
/// @param prev Ссылка на предыдущий слой в буфере
void step_self_sim(ro_problem& problem, single_var_moc_t& curr, single_var_moc_t& prev)
{
    vector<double>& layer = curr.vars.point_double[0];
    vector<double>& prev_layer = prev.vars.point_double[0];
    if (problem.Q > 0)
            {
                transform(
                    prev_layer.begin(), // Текущая точка
                    prev_layer.end() - 1,   // Предыдущая точка
                    layer.begin() + 1,
                    [](double a) { return a; }
                );
                layer[0] = problem.ro_left; // Дополнение слоя начальным условием
            }
            else
            {
                transform(
                    prev_layer.begin() + 1, // Текущая точка
                    prev_layer.end(),   // Предыдущая точка
                    layer.begin(),
                    [](double a) { return a; }
                );
                layer.pop_back();
                layer.push_back(problem.ro_right); // Дополнение слоя начальным условием
            }        
}


/// @brief Расчет следующего слоя уголком
/// @param problem Ссылка на структуру ro_problem
/// @param curr Ссылка на текущий слой в буфере
/// @param prev Ссылка на предыдущий слой в буфере
void step_corner(ro_problem& problem, single_var_moc_t& curr, single_var_moc_t& prev)
{
    vector<double>& layer = curr.vars.point_double[0];
    vector<double>& prev_layer = prev.vars.point_double[0];
    
    double sigma = problem.dt / problem.xStep;
    if (problem.Q > 0)
    {
        transform(
            prev_layer.begin(), // Текущая точка
            prev_layer.end() - 1,   // Предыдущая точка
            prev_layer.begin() + 1,
            layer.begin() + 1,
            [&sigma](double a, double b) { return (sigma * a) + ((1 - sigma) * b); }
        );
        layer[0] = problem.ro_left; // Дополнение слоя начальным условием
    }
    else
    {
        transform(
            prev_layer.begin() + 1, // Текущая точка
            prev_layer.end(),   // Предыдущая точка
            prev_layer.begin(),
            layer.begin(),
            [&sigma](double a, double b) { return (sigma * a) + ((1.0 - sigma) * b); }
        );
        layer.pop_back();
        layer.push_back(problem.ro_right); // Дополнение слоя начальным условием
    }
}


/// @brief Вывод текущего слоя в консоль
/// @param problem Ссылка на структуру ro_problem
/// @param curr Ссылка на текущий слой в буфере
/// @param step Шаг моделирования
void layer_to_console(ro_problem& problem, single_var_moc_t& curr, int step)
{
    vector<double>& layer = curr.vars.point_double[0];
    for (size_t i = 0; i < layer.size(); i++)
            printf("Шаг = %d    t = %.3F    x = %.3F    ro = %.3F\n",
                step,
                step * problem.dt,
                i * problem.xStep,
                layer[i]);
}


/// @brief Вывод описания задачи в консоль
/// @param problem Ссылка на структуру ro_problem
void task_to_console(ro_problem& problem)
{
    cout << "===== Описание задачи ===== " << endl;
    cout << "Рассчитываемая величина = " << problem.quantity << endl;
    cout << "Длина трубы = " << problem.L << endl;
    cout << "Время моделирования = " << problem.maxT << endl;
    cout << "Скорость потока = " << abs(4.0 * problem.Q / (PI * (problem.D - problem.d) * (problem.D - problem.d))) << endl;
    cout << "Шаг по координате = " << problem.xStep << endl;
    cout << "Шаг во времени = " << problem.dt << endl;
    cout << "Число узлов по времени = " << (int)(problem.maxT / problem.dt) << endl;
    cout << "Число узлов по координате = " << (problem.L / problem.xStep) + 1 << endl;
    cout << "=========================== " << endl;
}


/// @brief Запись текущего слоя в файл
/// @param problem Ссылка на структуру ro_problem
/// @param curr Ссылка на текущий слой в буфере
/// @param step Шаг моделирования
/// @param filename Имя файла
void layer_to_file(ro_problem& problem, single_var_moc_t& curr, int step, string filename = "data.txt")
{
    vector<double>& layer = curr.vars.point_double[0];
    std::ofstream out;
    if (step == 0)
        out.open(filename);
    else
        out.open(filename, ios::app);
    if (out.is_open())
        for (size_t i = 0; i < layer.size(); i++)
        {
            out << step * problem.dt << ',' << i * problem.xStep << ',' << layer[i] << endl;
        }
    out.close();
}


/// @brief Расчет плотности
/// @param problem Ссылка на структуру ro_problem
/// @param buffer Ссылка на буфер
/// @param verbose Флаг вывода в консоль
void calc_ro_problem(ro_problem& problem, custom_buffer_t<single_var_moc_t>& buffer, int method, bool verbose = false)
{
    // Начальное значение параметра
    problem.set_initial_layer(buffer);
    
    int numStepsTime = (problem.maxT / problem.dt);

    // Расчет слоев 
    for (size_t step = 0; step <= numStepsTime; step++)
    {
        if (verbose)
        {
            if (step == 0)
            {
                task_to_console(problem);
            }
            layer_to_console(problem, buffer.current(), step);
        }
        single_var_moc_t curr = buffer.current();
        single_var_moc_t prev = buffer.previous();

        layer_to_file(problem, buffer.current(), step);
        buffer.advance(1);

        switch (method)
        {
        case 0:
        {
            step_self_sim(problem, buffer.current(), buffer.previous());
            break;
        }
        case 1:
        {
            step_corner(problem, buffer.current(), buffer.previous());
            break;
        }
        }
        
    }
}


int main()
{
    // Установка кодировки для вывода кириллицы в консоль
    setlocale(LC_ALL, "Russian"); 

    // Для замера времени работы
    clock_t tStart = clock(); 

    // Создание структуры для решаемой задачи
    ro_problem problem = ro_problem();

    // Создание буфера для решаемой задачи
    custom_buffer_t<single_var_moc_t> buffer(2, int(problem.L / problem.xStep) + 1);

    // Решение задачи
    calc_ro_problem(problem, buffer, 1, true);

    printf("Время расчета: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    // Вызов скрипта построения графиков
    //system("py PlotPrecomputedGraph.py");

    return 0;
}