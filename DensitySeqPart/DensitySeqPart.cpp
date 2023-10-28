#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <algorithm>

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

using namespace std;

const int num_params = 2;
typedef profile_collection_t<num_params> many_profiles;
typedef vector<double> layer;

/// @brief  Физические и модельные характеристики трубы
struct PipeModel
{
    double D = 0.1; // Внешний диаметр трубы, м
    double d = 0.008; // Толщина стенок трубы; м
    double Q = 0.01;   // Объемный расход; м3/c
    double L = 10000; // Длина участка трубы
    double V; // Скорость потока
    double maxT = 600; // Время моделирования
    double dx = 100;// Шаг по координате
    double dt; // Шаг по времени

    PipeModel(double D, double d, double Q, double L, double maxT, double dx)
    {
        this->D = D;
        this->d = d;
        this->Q = Q;
        this->L = L;
        this->V = abs(4 * Q / (M_PI * pow(D - d, 2)));
        this->maxT = maxT;
        this->dx = dx;
        this->dt = dx / V;
    }
};

/// @brief Начальные условия задачи расчета определенного параметра
struct Problem
{
    PipeModel& pipe;
    string quantity = "Плотность"; // Рассчитываемый параметр
    double left_cond = 860; // Граничное условие (0; t)
    double right_cond = 860; // Граничное условие (L; t)
    double init_cond = 850; // Начальное условие (x, 0)

    Problem(PipeModel& pipeModel, string quantity, double left_cond, double right_cond, double init_cond)
        : pipe{ pipeModel }
    {
        this->quantity = quantity;
        this->left_cond = left_cond;
        this->right_cond = right_cond;
        this->init_cond = init_cond;
    }
};


class SolverMOC
{
private:
    layer& previousLayer;
    layer& nextLayer;

public:
    SolverMOC(layer& prev, layer& next)
        : previousLayer{ prev },
        nextLayer{ next }
    {
    }

/// @brief Расчет следующего слоя по свойству самоподобия
/// @param problem Ссылка на структуру задачи расчета параметра
    void step(Problem& problem)
    {
        if (problem.pipe.Q > 0)
        {
            transform(
                previousLayer.begin(), // Текущая точка
                previousLayer.end() - 1,   // Предыдущая точка
                nextLayer.begin() + 1,
                [](double a) { return a; }
            );
            nextLayer[0] = problem.left_cond; // Дополнение слоя начальным условием
        }
        else
        {
            transform(
                previousLayer.begin() + 1, // Текущая точка
                previousLayer.end(),   // Предыдущая точка
                nextLayer.begin(),
                [](double a) { return a; }
            );
            nextLayer.pop_back();
            nextLayer.push_back(problem.right_cond); // Дополнение слоя начальным условием
        }
    }
};


class SolverCorner
{
private:
    layer& previousLayer;
    layer& nextLayer;

public:
    SolverCorner(layer& prev, layer& next)
        : previousLayer{ prev },
        nextLayer{ next }
    {
    }

    /// @brief Расчет следующего слоя уголком
    /// @param problem Ссылка на структуру задачи расчета параметра
    void step(Problem& problem)
    {
        double sigma = 1 / problem.pipe.V;
        if (problem.pipe.Q > 0)
        {
            transform(
                previousLayer.begin(), // Текущая точка
                previousLayer.end() - 1,   // Предыдущая точка
                previousLayer.begin() + 1,
                nextLayer.begin() + 1,
                [&sigma](double a, double b) { return (sigma * a) + ((1 - sigma) * b); }
            );
            this->nextLayer[0] = problem.left_cond; // Дополнение слоя начальным условием
        }
        else
        {
            transform(
                previousLayer.begin() + 1, // Текущая точка
                previousLayer.end(),   // Предыдущая точка
                previousLayer.begin(),
                nextLayer.begin(),
                [&sigma](double a, double b) { return (sigma * a) + ((1.0 - sigma) * b); }
            );
            nextLayer.pop_back();
            nextLayer.push_back(problem.right_cond); // Дополнение слоя начальным условием
        }
    }
};


/// @brief Вывод текущего слоя в консоль
/// @param problem Ссылка на структуру задачи расчета параметра
/// @param currentLayer Ссылка на текущий слой в буфере
/// @param step Шаг моделирования
void layerToConsole(vector<Problem>& problems, many_profiles& currentLayer, int step)
{
    for (size_t i = 0; i < currentLayer.point_double[0].size(); i++)
    {
        printf("Шаг = %d    t = %.3F    x = %.3F    ",
            step,
            step * problems[0].pipe.dt,
            i * problems[0].pipe.dx);
        for (size_t problemNumber = 0; problemNumber < problems.size(); problemNumber++)
        {
            printf("    %s = %.3F",
                problems[problemNumber].quantity.c_str(),
                currentLayer.point_double[problemNumber][i]);
        }
        printf("\n");
    }
}

/// @brief Вывод описания задачи в консоль
/// @param problem Ссылка на структуру задачи расчета параметра
void problemToConsole(Problem& problem)
{
    cout << "======== Описание задачи ======== " << endl;
    cout << "Рассчитываемая величина = " << problem.quantity << endl;
    cout << "Длина трубы = " << problem.pipe.L << endl;
    cout << "Время моделирования = " << problem.pipe.maxT << endl;
    cout << "Скорость потока = " << problem.pipe.V << endl;
    cout << "Шаг по координате = " << problem.pipe.dx << endl;
    cout << "Шаг во времени = " << problem.pipe.dt << endl;
    cout << "Число узлов по времени = " << (int)(problem.pipe.maxT / problem.pipe.dt) << endl;
    cout << "Число узлов по координате = " << (problem.pipe.L / problem.pipe.dx) + 1 << endl << endl;
    cout << "Начальное условие = " << problem.init_cond << endl;
    cout << "Левое граничное условие = " << problem.left_cond << endl;
    cout << "Правое граничное условие = " << problem.right_cond << endl;
    cout << "==================================" << endl;
}


/// @brief Вывод текущего слоя в консоль
/// @param problem Ссылка на структуру задачи расчета параметра
/// @param currentLayer Ссылка на текущий слой в буфере
/// @param step Шаг моделирования
/// @param filename Имя файла
void layerToFile(vector<Problem>& problems, many_profiles& currentLayer, int step, string filename = "data.txt")
{
    std::ofstream out;
    if (step == 0)
    {
        out.open(filename);
        out << "время" << ',' << "координата";
        for (size_t problemNumber = 0; problemNumber < problems.size(); problemNumber++)
            out << ',' << problems[problemNumber].quantity;
        out << endl;
    }
    else
        out.open(filename, ios::app);
    if (out.is_open())
        for (size_t i = 0; i < currentLayer.point_double[0].size(); i++)
        {
            out << step * problems[0].pipe.dt << ',' << i * problems[0].pipe.dx;
            for (size_t problemNumber = 0; problemNumber < problems.size(); problemNumber++)
                out << ',' << currentLayer.point_double[problemNumber][i];
            out << endl;
        }
    out.close();
}


/// @brief Оболочка для цикла расчета плотности за время maxT
/// @param Problem Ссылка на структуру Problem
/// @param buffer Ссылка на буфер
/// @param verbose Флаг вывода в консоль
void сalculateProblems(vector<Problem>& problems, custom_buffer_t<many_profiles>& buffer)
{
    int numStepsTime = (problems[0].pipe.maxT / problems[0].pipe.dt);
    
    for (size_t step = 0; step <= numStepsTime; step++)
    {
        for (size_t problemIndex = 0; problemIndex < problems.size(); problemIndex++)
        {
            Problem& problem = problems[problemIndex];
            layer& previousLayer = buffer.previous().point_double[problemIndex];
            layer& currentLayer = buffer.current().point_double[problemIndex];

            if (step == 0)
            {
                previousLayer = layer(previousLayer.size(), problem.init_cond);
                problemToConsole(problem);
            }

            SolverCorner solver(previousLayer, currentLayer);
            solver.step(problem); 

            if (problemIndex + 1 == problems.size())
            {
                layerToConsole(problems, buffer.previous(), step);
                layerToFile(problems, buffer.previous(), step);
            }
        }
        buffer.advance(1);
    }  
}


int main()
{
    double D = 0.1;
    double d = 0.008;
    double Q = 0.01;
    double L = 10000;
    double maxT = 5000;
    double dx = 100;

    PipeModel pipe(D, d, Q, L, maxT, dx);

    double leftDensity = 860;
    double rightDensity = leftDensity;
    double initialDensity = 850;

    double leftSulfur = 0.9;
    double rightSulfur = leftSulfur;
    double initialSulfur = 0.6;

    vector<Problem> problems{
        Problem(pipe, "плотность", leftDensity, rightDensity, initialDensity),
        Problem(pipe, "сера", leftSulfur, rightSulfur, initialSulfur)
    };
  
    // Установка кодировки для вывода кириллицы в консоль
    setlocale(LC_ALL, "Russian"); 

    // Для замера времени работы
    clock_t tStart = clock(); 


    // Создание буфера для решаемой задачи
    custom_buffer_t<many_profiles> buffer(2, int(problems[0].pipe.L / problems[0].pipe.dx) + 1);
    // Решение задачи
    сalculateProblems(problems, buffer);

    printf("Время расчета: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    // Вызов скрипта построения графиков
    //system("py PlotPrecomputedGraph.py");

    return 0;
}