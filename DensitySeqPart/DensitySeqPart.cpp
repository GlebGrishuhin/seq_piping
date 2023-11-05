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
typedef profile_collection_t<num_params> many_profiles_t;
typedef vector<double> layer_t;

/// @brief  Физические и модельные характеристики трубы
struct pipe_model_t
{
    /// @brief Внешний диаметр трубы, м
    double D = 0.1;
    /// @brief Толщина стенок трубы; м
    double d = 0.008;
    /// @brief Объемный расход; м3/c
    double Q = 0.01; 
    /// @brief Длина участка трубы
    double L = 10000;
    /// @brief Скорость потока
    double V;
    /// @brief Время моделирования
    double maxT = 600;
    /// @brief Шаг по координате
    double dx = 100;
    /// @brief Шаг по времени
    double dt; 

    pipe_model_t(double D, double d, double Q, double L, double maxT, double dx)
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
struct problem_t
{
    pipe_model_t& pipe;
    /// @brief Рассчитываемый параметр
    string quantity = "Плотность"; 
    /// @brief Граничное условие (0; t)
    double left_condition = 860;
    /// @brief Граничное условие (L; t)
    double right_condition = 860;
    /// @brief Начальное условие (x, 0)
    double initial_condition = 850; 

    problem_t(pipe_model_t& pipe_model_t, string quantity, 
        double left_condition, double right_condition, double initial_condition)
        : pipe{ pipe_model_t }
    {
        this->quantity = quantity;
        this->left_condition = left_condition;
        this->right_condition = right_condition;
        this->initial_condition = initial_condition;
    }
};


class solver_t
{
private:
    layer_t& previous_layer;
    layer_t& next_layer;

public:
    solver_t(layer_t& prev, layer_t& next)
        : previous_layer{ prev },
        next_layer{ next }
    {
    }

/// @brief Расчет следующего слоя по свойству самоподобия
/// @param problem_t Ссылка на структуру задачи расчета параметра
    void step_moc(problem_t& problem)
    {
        if (problem.pipe.Q > 0)
        {
            /// @brief Формирование нового слоя путем сдвига предыдущего
            transform(
                previous_layer.begin(),
                previous_layer.end() - 1,
                next_layer.begin() + 1,
                [](double a) { return a; }
            );
            
            next_layer[0] = problem.left_condition; // 
        }
        else
        {
            /// @brief Формирование нового слоя путем сдвига предыдущего
            transform(
                previous_layer.begin() + 1,
                previous_layer.end(),
                next_layer.begin(),
                [](double a) { return a; }
            );
            /// @brief Дополнение слоя начальным условием
            next_layer.pop_back();
            next_layer.push_back(problem.right_condition);
        }
    }

    /// @brief Расчет следующего слоя уголком
    /// @param problem_t Ссылка на структуру задачи расчета параметра
    void step_corner(problem_t& problem)
    {
        double sigma = 1 / problem.pipe.V;
        if (problem.pipe.Q > 0)
        {
            transform(
                previous_layer.begin(), // Текущая точка
                previous_layer.end() - 1,   // Предыдущая точка
                previous_layer.begin() + 1,
                next_layer.begin() + 1,
                [&sigma](double a, double b) { return (sigma * a) + ((1 - sigma) * b); }
            );
            this->next_layer[0] = problem.left_condition; // Дополнение слоя начальным условием
        }
        else
        {
            transform(
                previous_layer.begin() + 1, // Текущая точка
                previous_layer.end(),   // Предыдущая точка
                previous_layer.begin(),
                next_layer.begin(),
                [&sigma](double a, double b) { return (sigma * a) + ((1.0 - sigma) * b); }
            );
            next_layer.pop_back();
            next_layer.push_back(problem.right_condition); // Дополнение слоя начальным условием
        }
    }
};


/// @brief Вывод текущего слоя в консоль
/// @param problem Ссылка на структуру задачи расчета параметра
/// @param currentLayer Ссылка на текущий слой в буфере
/// @param step Шаг моделирования
void layer_to_console_(vector<problem_t>& problem, many_profiles_t& currentLayer, int step)
{
    for (size_t i = 0; i < currentLayer.point_double[0].size(); i++)
    {
        printf("Шаг = %d    t = %.3F    x = %.3F    ",
            step,
            step * problem[0].pipe.dt,
            i * problem[0].pipe.dx);
        for (size_t problem_index = 0; problem_index < problem.size(); problem_index++)
        {
            printf("    %s = %.3F",
                problem[problem_index].quantity.c_str(),
                currentLayer.point_double[problem_index][i]);
        }
        printf("\n");
    }
}

/// @brief Вывод описания задачи в консоль
/// @param problem_t Ссылка на структуру задачи расчета параметра
void problem_to_console_(problem_t& problem)
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
    cout << "Начальное условие = " << problem.initial_condition << endl;
    cout << "Левое граничное условие = " << problem.left_condition << endl;
    cout << "Правое граничное условие = " << problem.right_condition << endl;
    cout << "==================================" << endl;
}


/// @brief Вывод текущего слоя в консоль
/// @param problem Ссылка на структуру задачи расчета параметра
/// @param currentLayer Ссылка на текущий слой в буфере
/// @param step Шаг моделирования
/// @param filename Имя файла
void layer_to_file(vector<problem_t>& problem, many_profiles_t& currentLayer, int step, string filename = "data.txt")
{
    std::ofstream out;
    if (step == 0)
    {
        out.open(filename);
        out << "время" << ',' << "координата";
        for (size_t problem_index = 0; problem_index < problem.size(); problem_index++)
            out << ',' << problem[problem_index].quantity;
        out << endl;
    }
    else
        out.open(filename, ios::app);
    if (out.is_open())
        for (size_t i = 0; i < currentLayer.point_double[0].size(); i++)
        {
            out << step * problem[0].pipe.dt << ',' << i * problem[0].pipe.dx;
            for (size_t problem_index = 0; problem_index < problem.size(); problem_index++)
                out << ',' << currentLayer.point_double[problem_index][i];
            out << endl;
        }
    out.close();
}


/// @brief Оболочка для цикла расчета плотности за время maxT
/// @param problems Ссылка на вектор структур problem
/// @param buffer Ссылка на буфер
/// @param method Используемый метод: 1 - метод хар-ик; 2 - уголок;
void сalculate_problem_(vector<problem_t>& problems, ring_buffer_t<many_profiles_t>& buffer, int method = 1)
{
    int numstep_mocsTime = (problems[0].pipe.maxT / problems[0].pipe.dt);
    
    for (size_t step_moc = 0; step_moc <= numstep_mocsTime; step_moc++)
    {
        for (size_t problem_index = 0; problem_index < problems.size(); problem_index++)
        {
            problem_t& problem = problems[problem_index];
            layer_t& previous_layer = buffer.previous().point_double[problem_index];
            layer_t& currentLayer = buffer.current().point_double[problem_index];

            if (step_moc == 0)
            {
                previous_layer = layer_t(previous_layer.size(), problem.initial_condition);
                problem_to_console_(problem);
            }

            solver_t solver(previous_layer, currentLayer);
            
            switch (method)
            {
            case 1:
            {
                solver.step_moc(problem);
                break;
            }
            case 2:
            {
                solver.step_corner(problem);
                break;
            }
            }
            
            if (problem_index + 1 == problems.size())
            {
                layer_to_console_(problems, buffer.previous(), step_moc);
                layer_to_file(problems, buffer.previous(), step_moc);
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

    pipe_model_t pipe(D, d, Q, L, maxT, dx);

    double left_condition_density = 860;
    double right_condition_density = left_condition_density;
    double initial_condition_density = 850;

    double left_condition_sulfur = 0.9;
    double right_condition_sulfur = left_condition_sulfur;
    double initial_condition_sulfur = 0.6;

    vector<problem_t> problem{
        problem_t(pipe, "плотность", left_condition_density, right_condition_density, initial_condition_density),
        problem_t(pipe, "сера", left_condition_sulfur, right_condition_sulfur, initial_condition_sulfur)
    };
  
    // Установка кодировки для вывода кириллицы в консоль
    setlocale(LC_ALL, "Russian"); 

    // Для замера времени работы
    clock_t tStart = clock(); 


    // Создание буфера для решаемой задачи
    ring_buffer_t<many_profiles_t> buffer(2, int(problem[0].pipe.L / problem[0].pipe.dx) + 1);
    // Решение задачи
    сalculate_problem_(problem, buffer);

    printf("Время расчета: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    // Вызов скрипта построения графиков
    system("py PlotPrecomputedGraph.py");

    return 0;
}