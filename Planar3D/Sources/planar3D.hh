#pragma once

#define _USE_MATH_DEFINES

#include <cmath>
#include <string>
#include <vector>

#if !defined M_PI
    #define M_PI 3.14159265358979323846264338327950288
#endif

//глобальные переменные, необходимые в разных фрагментах расчётной программы

/// Константа для инициализации неизвестных значений
#define INITIAL -666

/// Тип внешних для трещины сеточных элементов
#define OUTSIDE 0

/// Тип внутренних для трещины сеточных элементов, близких к фронту трещины
#define RIBBON  1

/// Тип внутренних для трещины сеточных элементов
#define CHANNEL 2

/// Тип x-производной сеточного элемента с одним активным соседом
#define LEFT   -1

/// Тип x-производной сеточного элемента без активных соседей
#define ALONE   0

/// Тип x-производной сеточного элемента с одним активным соседом
#define RIGHT   1

/// Тип x-производной сеточного элемента с двумя активными соседями
#define REGULAR 2

/// Тип y-производной сеточного элемента с одним активным соседом
#define BOTTOM -1

/// Тип y-производной сеточного элемента без активных соседей
#define ALONE   0

/// Тип y-производной сеточного элемента с одним активным соседом
#define TOP     1

/// Тип y-производной сеточного элемента с двумя активными соседями
#define REGULAR 2

extern int regime;

extern double initialRadius;
extern double dx;
extern double dy;
extern double bn;
extern double mu;
extern double Q;
extern double E;
extern double C;
extern double Kic;
extern double alpha;
extern double Amu;
extern double epsilon;
extern double dt_step;		///<Переменная для сохранения шага по времени для передачи в расчет проппанта

/// Координата источника по абсциссе
extern size_t i00;
/// Координата источника по ординате
extern size_t j00;

/// \brief Класс граничного элемента
/// \details Класс, задающий граничный элемент (координаты – индексы ячейки) с
/// оператором сравнения двух граничных элементов
class Ribbon{
    public:
        size_t i;
        size_t j;

        Ribbon(const size_t i, const size_t j): i(i), j(j){};

        bool operator == (const Ribbon ribbon) const{
            return this->i == ribbon.i && this->j == ribbon.j;
        }
};

/// \brief Класс ячейки сетки
/// \details Класс, задающий ячейку расчётной сетки (координаты центра, тип,
/// тип производной вдоль обоих орт, время активации)
class Cell{
    public:
        double x = 0;
        double y = 0;

        int type = INITIAL;
        int xDerivative = INITIAL;
        int yDerivative = INITIAL;

        double activationTime = -1.;

        void setCoordinates(const double x, const double y){
            this->x = x;
            this->y = y;
        }
};

/// \brief Функция возвращает название режима
std::string regimeName();

/// \brief Функция выводит лого программы
void printLogo();

/// \brief Основная расчётная функция
int planar3D(
    double T1,
    const double ts,
    double cellSize,
    double meshSize,
    int meshScalingCounter,
    std::string pathToLayersFile,
    std::string pathToInjectionFile,
    const std::string pathToImportFolder,
    const bool runningFromGUI,
    const bool saveSteps,
    const std::size_t numberOfThreads
);



/// \brief Функция экспорта данных в json формате
/// 
void SaveJson(
	const std::string filename,		//Имя файла (без разрешения) для сохранения выходного файла json
	std::vector<double> &Wk,
	std::vector<double> &pressure,
	std::vector<double>&concentration,
	std::vector< std::vector<size_t> > &index, 
	double Time,
	double xSize,
	double ySize,
	double wn,
	double DproppantInjection						//Величина мнгновенной плотности проппанта (кг/(м3*dt)
	);