/*!
/// \file planar3D.cc
/// \author Egor Starobinskii, Sergey Khlopin
/// \date 17 Spt 2018


файл с функцией planar3D - главной функцией расчета и вызова функций различных режимов
Данный файл содержит в себе определения основных
классов, используемых в демонстрационной программе
*/

#include <vector>
#include <iostream>
#include <algorithm>

#include "Include/nlohmann/json.hpp"

#include "ai.hh"
#include "rhs3D.hh"
#include "initialData.hh"
#include "findRibbons.hh"
#include "findNewRibbons.hh"
#include "setDerivatives.hh"
#include "findActiveElements.hh"
#include "buildInfluenceMatrix.hh"
#include "buildPartialInfluenceMatrix.hh"

#ifdef OPENMP
    #include <omp.h>
    int getNumberOfThreads(){
        int counter = 0;
        #pragma omp parallel
        {
            counter = omp_get_num_threads();
        }
        return counter;
    }
#else
    int getNumberOfThreads(){return 1;}
    void omp_set_num_threads(std::size_t numberOfThreads){};
#endif

///глобальные переменные, необходимые в разных фрагментах расчётной программы

int regime;

double initialRadius;
/*!
double dx шаг вычислений по координате x
\brief шаг по координате х
*/
double dx;
double dy;
double bn;
double mu;
double Q;
double E;
double C;
double Kic;
double alpha;
double Amu;
double epsilon;
double dt_step;

size_t i00;
size_t j00;

std::vector<double> x;
std::vector<double> y;

/// \details Функция возвращает название текущего режима распространения
/// трещины
std::string regimeName(){
    switch(regime){
        case -1:
            return "leak-off dominated";
        case 0:
            return "toughness dominated";
        case 1:
            return "viscosity dominated";
        default:
            return "unknown";
    }
}











void calculateCrackGeometry(
    std::vector< std::vector<Cell> > &mesh,
    std::vector< std::vector<double> > &distances,
    double &length,
    double &height
){
    length = 0.;
    height = 0.;
    double heightUp = 0.;
    double heightDown = 0.;

    const std::size_t xSize = distances.size();
    const std::size_t ySize = distances[0].size();

    for(size_t i = 0; i < xSize; ++i){
        if(RIBBON == mesh[i][j00].type){
            length = distances[i][j00]
                + std::sqrt(std::pow(mesh[i][j00].x, 2)
                + std::pow(mesh[i][j00].y, 2));
            length *= 2.;

            break;
        }
    }
    for(size_t j = 0; 2 * j < ySize; ++j){
        if(
            RIBBON == mesh[i00][j].type
            || RIBBON == mesh[i00][ySize - j - 1].type
        ){
            if(0. == heightUp && RIBBON == mesh[i00][j].type ){
                heightUp = distances[i00][j]
                    + std::sqrt(std::pow(mesh[i00][j].x, 2)
                    + std::pow(mesh[i00][j].y, 2));
            }

            if(0. == heightDown && RIBBON == mesh[i00][ySize - j - 1].type){
                heightDown = distances[i00][ySize - j - 1]
                    + std::sqrt(std::pow(mesh[i00][ySize - j - 1].x, 2)
                    + std::pow(mesh[i00][ySize - j - 1].y, 2));
            }

            if(0. != heightUp && 0 != heightDown){
	        height = heightUp + heightDown;

                break;
	    }
        }
    }
}
void printEfficiency(
    std::vector<double> &opening
){
    /// \todo Вынести расчёт эффективности в эту функцию
}




















/*!
\brief __printLogo__ - функция отображения логотипа программы в консоли
\details Функция выводит текстовое лого программы
*/
void printLogo(){
    std::cout << " ____  _                       _____ ____     ____ _     "
        << "___ " << std::endl;
    std::cout << "|  _ \\| | __ _ _ __   __ _ _ _|___ /|  _ \\   / ___| |   "
        << "|_ _|" << std::endl;
    std::cout << "| |_) | |/ _` | '_ \\ / _` | '__||_ \\| | | | | |   | |    "
        << "| | " << std::endl;
    std::cout << "|  __/| | (_| | | | | (_| | |  ___) | |_| | | |___| |___ | "
        << "| " << std::endl;
    std::cout << "|_|   |_|\\__,_|_| |_|\\__,_|_| |____/|____/   "
        << "\\____|_____|___|" << std::endl;
    std::cout << std::endl << "Developed by REC \"Gazpromneft-Polytech\"."
        << std::endl << std::endl;
}

/*!
\detailed основная расчётная функция, считающая распространение трещины в
режиме доминирующей вязкости по планарной модели и сохраняющая данные о
результатах расчёта в виде текстовых файлов


\param[in] modelingTime Время окончания расчета
\param[in] ts Масштаб времени (число реальных секунд в расчетном времени)
\param[in] cellSize	Размер стороны квадратной ячейки в метрах
\param[in] meshSize	Число ячеек в расчётной области вдоль горизонтальной оси
\f$meshSize=\frac{numCels}{2*\left( cellSize-1 \right)}\f$
\param[in] meshScalingCounter Сколько раз допускается масштабировать сетку
\param[in] pathToLayersFile Путь (относительный или абсолютный) к расположению
файла слоёв
\param[in] pathToInjectionFile Путь (относительный или абсолютный) к
расположению файла параметров закачки
\param[in] runningFromGUI Флаг запуска программы через графический интерфейс
\param[in] saveSteps Флаг сохранения значений промежуточных расчетов на каждом
временном шаге
\param[in] numberOfThreads флаг сохранения значений промежуточных расчетов на
каждом временном шаге
*/
int planar3D(
    double modelingTime,
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
){

    if(!runningFromGUI){
        printLogo();
    }

    /// Проверяем соответствие структуры рабочей директории
    if(!ai::folderExists("./Results")){
        std::cerr << "Cannot find './Results', create directory and restart."
            << std::endl;

        return 21;
    }

    if(saveSteps){
        if(!ai::folderExists("./Results/Concentration")){
            std::cerr << "Cannot find './Results/Concentration', create "
                << "directory and restart." << std::endl;

            return 21;
        }
        if(!ai::folderExists("./Results/Opening")){
            std::cerr << "Cannot find './Results/Opening', create directory "
                << "and restart." << std::endl;

            return 21;
        }
        if(!ai::folderExists("./Results/Pressure")){
            std::cerr << "Cannot find './Results/Pressure', create directory "
                << "and restart." << std::endl;

            return 21;
        }
    }

    if(
        std::string() != pathToImportFolder
        && !ai::folderExists(pathToImportFolder)
    ){
        std::cerr << "Cannot find '" << pathToImportFolder << "' to import "
            << "data, create directory and restart." << std::endl;

        return 21;
    }

    /// Подгружаем входные данные

    std::vector< std::vector<double> > modelSolution;
    std::vector< std::vector<double> > layers;
    std::vector< std::vector<double> > injection;

    nlohmann::json importData;

    if(!importInitialData(pathToImportFolder, importData)){
        return 22;
    }

    if(0 < importData.size()){
        /// \todo Настроить полный импорт даныых в соответсвии с новой
        /// разметкой файла parameters.json
        bn = importData["fluid"]["0"]["rheology index"].get<double>();
        mu = importData["fluid"]["0"]["dynamic viscosity"].get<double>();
        modelingTime = importData["time"].get<double>();
        pathToLayersFile = std::string();
        pathToInjectionFile = std::string();
    }

    if(
        !setInitialData(bn, initialRadius, modelSolution,
            pathToLayersFile, layers, pathToInjectionFile, injection
        )
    ){
        return 22;
    }

    for(std::size_t i = 0; i < layers.size(); ++i){
        layers[i][2] *= std::pow(10, 6);
        layers[i][3] *= std::pow(10, 9);
        layers[i][5] *= std::pow(10, -6);
    }

    for(std::size_t i = 0; i < layers.size(); ++i){
        if(0. > layers[i][0] * layers[i][1]){
            E = layers[i][3] / (1 - std::pow(layers[i][4], 2));
            /// \todo Добавить учёт слоистости по коэффициенту Картера
            C = layers[i][5];

            break;
        }
    }

    /// Проверяем параметры закачки, добавляем данные о свойствах жидкости

    if(!recalculateInjection(injection, modelingTime)){
        return 23;
    }

    size_t injectionIndex = 0;

    double timeToChangeInjection = modelingTime * 1.1;
    double proppantInjection = injection[injectionIndex][3];
    Q = injection[injectionIndex][2];
    bn = injection[injectionIndex][4];
    mu = injection[injectionIndex][5];

    ++injectionIndex;

    if(injection.size() > injectionIndex){
        timeToChangeInjection = injection[injectionIndex][0];
    }


    double T0 = 1.;
    const double wn = std::pow(
        mu * 2. * std::pow(2. * (2. * bn + 1.) / bn, 1. / bn)
        / (E * std::pow(ts, bn)), 1. / (bn + 2.)
    );
    const double gammaR = 1. / 3. * (1. + bn / (bn + 2.));

    /// Выводим входные параметры
    if(!runningFromGUI){
        /// \todo Переосмыслить вывод начальных параметров
        std::cout << std::endl;
        std::cout << "Incoming parameters:" << std::endl
            << "  Q = " << Q << ", mu = " << mu << ", "
            << "n = " << bn << ";" << std::endl
            << "  E\' = " << E << ", " << "C = " << C << ","
            << " Kic = " << Kic << ";" << std::endl
            << "  time = " << modelingTime << ", ts = " << ts << ";"
            << std::endl;
    }

    /// Масштабируем величины

    Q *= ts / (wn * 60.);

    //// Определяем момент времени, вплоть до которого будет использоваться
    /// автомодельное решение для радиальной плоской трещины в однородной среде
    /// (диаметр трещины – 90% от мощности центрального слоя)

    double requiredRadius = -1.;

    for(std::size_t i = 0; i < layers.size(); ++i){
        if(0. > layers[i][0] * layers[i][1]){
            requiredRadius = 0.45 * std::abs(layers[i][0] - layers[i][1]);
        }
    }

    T0 = std::pow(
        requiredRadius / (initialRadius * std::pow(Q, 1. / 3.)),
        1. / gammaR
    );

    T0 = ai::min(T0, timeToChangeInjection);


    initialRadius *= std::pow(Q, 1. / 3.) * std::pow(T0, gammaR);

    /// Пересчитываем автомодельное решение в соответсвии с параметрами закачки

    std::vector<double> zP;
    std::vector<double> openingAtTheStart;

    for(size_t i = 0; i < modelSolution.size(); ++i){
        modelSolution[i][1] *= std::pow(Q, 1. / 3.);
        openingAtTheStart.push_back(
            modelSolution[i][1] * std::pow(T0, (1. - 2. * bn / (bn + 2.)) / 3.)
        );
        zP.push_back(modelSolution[i][0]);
    }

    modelSolution.clear();

    const double initialVelocity = initialRadius * gammaR
        * std::pow(T0, gammaR - 1);

    // if(0 < importData.size()){
    //     cellSize = round(
    //         initialRadius / importData["mesh"]["cell"]["length"].get<double>()
    //     );
    //     meshSize = floor(
    //         (importData["mesh"]["length"].get<double>() - 1.) / cellSize
    //     );
    // }

    /// \todo проверить вычисление dt(mu)

    const double dt = 0.0001 * std::pow(5. / floor(initialRadius / cellSize), 2)
        * std::pow(mu / 0.4, 0.5);
    dt_step = dt;                       //сохранили шаг по времени для использования в расчете проппанта (Света)

    /// Выводим входные параметры (продолжение)
    if(!runningFromGUI){
        std::cout << "  mesh = " << meshSize << ", cell = " << cellSize << ","
            << " scaling = " << meshScalingCounter << "." << std::endl;
        std::cout << "  required radius = " << requiredRadius << ", it's time = "
            << T0 << ";"<< std::endl;
        std::cout << "Initial regime: " << regimeName() << "." << std::endl;

        std::cout << "Number of cell per initial radius:"
            << floor(initialRadius / cellSize) << std::endl;

        std::cout << "Time step:"
            << dt << std::endl;
    }

    /// Задаём расчётную область
    dx = cellSize;
    dy = dx;

    double axMax = meshSize * dx;

    x.clear();
    y.clear();

    for(double i = 0; i <= axMax + epsilon; i += dx){
        x.push_back(i);
    }
    for(double i = -axMax; i <= axMax + epsilon; i += dy){
        y.push_back(i);
    }

    const size_t xSize = x.size();
    const size_t ySize = y.size();

    std::vector< std::vector<Cell> > mesh(xSize);

    for(size_t i = 0; i < xSize; ++i){
        /// \todo Отказаться от сетки как класса, использовать вектора
        mesh[i].resize(ySize);

        for(size_t j = 0; j < ySize; ++j){
            mesh[i][j].setCoordinates(x[i], y[j]);
        }
    }

    /// Задаём положение точечного источника закачки

    i00 = 0;
    j00 = floor(0.5 * ySize);

    /// Сохраняем параметры расчёта

    saveInitialData(
        "./Results/parameters",
        modelingTime,
        ySize,
        xSize,
        dy,
        dx,
        injection,
        layers
    );

    /// Пересчитываем контраст напряжений по построенной расчётной сетке

    std::vector<double> stress;

    if(!recalculateStressContrast(layers, stress, y)){
        return 31;
    }

    /// Масштабируем выличины (продолжение)
    for(std::size_t i = 0; i < stress.size(); ++i){
        stress[i] /= (wn * E);
    }
    for(std::size_t i = 0; i < injection.size(); ++i){
        injection[i][2] *= ts / (wn * 60.);
    }

    C *= std::sqrt(ts) / wn;
    Kic /= std::pow(mu * E * std::pow(E / ts, bn), 1. / (bn + 2.));








    alpha = 2. / (bn + 2.);
    const double BAlpha = 0.25 * alpha * tan(0.5 * M_PI - M_PI * (1 - alpha));
    Amu = std::pow(BAlpha * (1 - alpha), -0.5 * alpha);

    std::vector<double> zeroVectorX(xSize, 0.);
    std::vector<double> zeroVectorY(ySize, 0.);
    std::vector<size_t> zeroSizeTVectorY(ySize, 0);
    std::vector<double> zeroVectorXY(xSize * ySize, 0.);
    std::vector< std::vector<double> > zeroMatrixXY(xSize, zeroVectorY);

    std::cout << std::endl;
    std::cout << "Building influence matrix... ";

    std::vector< std::vector<double> > influenceMatrix(xSize * ySize,
        zeroVectorXY);

    buildInfluenceMatrix(influenceMatrix, xSize, ySize);

    std::cout << "OK." << std::endl;

    std::vector<double> Wk = zeroVectorXY;
    std::vector< std::vector<size_t> > index(xSize, zeroSizeTVectorY);


    size_t p = 0;
    size_t m = 0;
    for(size_t i = 0; i < xSize; ++i){
        for(size_t j = 0; j < ySize; ++j){
            p = 0;

            for(size_t k = 0; k < xSize; ++k){
                for(size_t l = 0; l < ySize; ++l){
                    ++p;
                }
            }

            index[i][j] = m;
            Wk[m] = getInitialOpening(mesh[i][j].x, mesh[i][j].y,
                initialRadius, zP, openingAtTheStart);

            ++m;
        }
    }
	ai::saveVector("./iniopen",Wk);
	
    const double dMin1 = std::sqrt(std::pow(0.5 * dx, 2)
        + std::pow(1.5 * dy, 2));
    const double dCenter1 = dx;
    const double dMax1 = std::sqrt(std::pow(0.5 * dx, 2)
        + std::pow(0.5 * dy, 2));
    const double dMin2 = std::sqrt(std::pow(1.5 * dx, 2)
        + std::pow(1.5 * dy, 2));
    const double dCenter2 = std::sqrt(2) * dx;

    std::cout << "Mesh: " << mesh.size() << "x" << mesh[0].size() << "."
        << std::endl;

    std::vector< std::vector<size_t> >
        activeElements = findActiveElements(mesh, initialRadius);

    std::cout << "Active elements: " << activeElements.size() << "."
        << std::endl;

    setDerivatives(mesh, activeElements);

    for(size_t k = 0; k < activeElements.size(); ++k){
        const size_t i = activeElements[k][0];
        const size_t j = activeElements[k][1];

        if(
            RIBBON != mesh[i][j].type
            && REGULAR == mesh[i][j].xDerivative
            && REGULAR == mesh[i][j].yDerivative
        ){
            mesh[i][j].type = CHANNEL;
        }

        if(
            RIBBON != mesh[i][j].type
            && LEFT == mesh[i][j].xDerivative
            && REGULAR == mesh[i][j].yDerivative
        ){
            mesh[i][j].type = CHANNEL;
        }

        if(OUTSIDE < mesh[i][j].type){
            mesh[i][j].activationTime = T0 - dt;
        }
    }

    std::vector< std::vector<double> > distances;
    std::vector<Ribbon> ribbons = findRibbons(mesh, activeElements,
        distances, initialRadius, dMax1, 2. * dx);

    std::cout << "Ribbons: " << ribbons.size() << "." << std::endl;

    std::vector<double> WkNew;
    std::vector< std::vector<double> > partialInfluenceMatrix;

    std::cout << "Building partial influence matrix... ";

    buildPartialInfluenceMatrix(influenceMatrix, activeElements, Wk, WkNew,
        partialInfluenceMatrix, index);

    std::cout << "OK." << std::endl;

    double height = 0.;
    double length = 0.;
    std::vector< std::vector<double> > fracture;


    size_t stepToCheck = round(dx / (dt * initialVelocity) / 20.);
    size_t step = 0;

    double T = T0;
    double savedTime = T0;

    bool meshIsNotExhausted = true;

    std::vector<double> rhsCE;
    std::vector<double> rhsCEs;
    std::vector<double> pressure;
    std::vector<double> savedDistances;

    std::vector< std::vector<double> > velocities = zeroMatrixXY;
    std::vector< std::vector<double> > openingAtTheStep = zeroMatrixXY;
    std::vector< std::vector<double> > pressureAtTheStep = zeroMatrixXY;

	std::vector<double> concentration = zeroVectorXY;            //Вектор концентраций проппанта Света
    std::vector<double>  concentration_temp = zeroVectorXY;        //Вектор новых концентраций проппанта Света
	std::vector< std::vector<double> > concentrationAtTheStep = zeroMatrixXY;

    omp_set_num_threads(numberOfThreads);

    std::cout << std::endl << "Running threads... OK. Size: "
        << getNumberOfThreads() << "." << std::endl;

    std::cout << std::endl << "Starting calculations..." << std::endl;
    if(runningFromGUI){
        std::cout << "Progress: 0.0" << std::endl;
    }else{
        ai::showProgressBar(0.);
    }

    auto startTime = ai::time();

    while(modelingTime >= T && meshIsNotExhausted){
        calculatePressure(
            pressure,
            index,
            activeElements,
            partialInfluenceMatrix,
            WkNew,
            stress,
            Wk.size()
        );

         calculateOpeningSpeedProp(
             Wk,
             rhsCE,
             rhsCEs,
             pressure,
             mesh,
             index,
             activeElements,
             concentration,
             concentration_temp,
             proppantInjection,
             T
        );

        // if (T > 17)
        //  SaveJson("Results/JSON/rezult", Wk, pressure, concentration, index, T, xSize, ySize, wn, proppantInjection);


 /*       calculateOpeningSpeed(Wk, rhsCE, rhsCEs, pressure, mesh, index,
            activeElements, T);*/

        if(0 == regime){
            double maxDeltaDistance = 0.;

            for(size_t k = 0; k < ribbons.size(); ++k){
                const size_t i = ribbons[k].i;
                const size_t j = ribbons[k].j;

                 distances[i][j] = 0.25 *  sqrt(0.5 * M_PI)
                    * Wk[index[i][j]] / Kic;

                if(epsilon > distances[i][j]){
                    
					distances[i][j] = epsilon;
                }else{
                    distances[i][j] *= distances[i][j];
					std::cout<<" dx = "<<dx<<std::endl;
					if (distances[i][j] > 2. * sqrt(2.) * dx)
					{
						std::cout <<"dist is LARGE"<<std::endl;
					}
                }
				if(savedDistances.size() > 0)
                maxDeltaDistance = ai::max(
                    distances[i][j] - savedDistances[k],
                    maxDeltaDistance
                );
            }
			// ai::saveMatrix("./dist", distances);
			// std::cout<<"maxDeltaDistance = "<<maxDeltaDistance<<std::endl;
            if(epsilon < maxDeltaDistance){
                // stepToCheck = std::round(
                    // dx * (T - savedTime) / (20 * dt * maxDeltaDistance)
                // );
				if(maxDeltaDistance > 2. * sqrt(2.) * dx)
				stepToCheck = step;
			
				if (savedTime > 0.)
				stepToCheck = std::round(dx * (T-savedTime)/ (20. * dt * maxDeltaDistance));
            }else{
                
				stepToCheck = 200;
            }
			
			// std::cout<<"stepToCheck = "<<stepToCheck<<std::endl;
        }else{
            calculateVelocity(
                velocities,
                index,
                ribbons,
                distances,
                Wk
            );

            for(size_t i = 0; i < xSize; ++i){
                for(size_t j = 0; j < ySize; ++j){
                    distances[i][j] += velocities[i][j] * dt;
                }
            }
        }

        for(size_t i = 0; i < rhsCE.size(); ++i){
            Wk[i] = ai::max(Wk[i] + rhsCE[i] * dt, 0.);
        }

        for(size_t i = 0; i < rhsCEs.size(); ++i){
            WkNew[i] = ai::max(WkNew[i] + rhsCEs[i] * dt, 0.);
        }

        const size_t savedSize = activeElements.size();

        if(step == stepToCheck){
            if(runningFromGUI){
                std::cout << "Progress: " << (T - T0) / (modelingTime - T0)
                    << std::endl;
            }else{
                ai::showProgressBar((T - T0) / (modelingTime - T0));
            }

            step = 0;

            if(0 != regime){
                stepToCheck = std::round(dx / (20 * dt * ai::max(velocities)));
            }

            std::vector<Ribbon> oldRibbons = ribbons;

            for(size_t k = 0; k < ribbons.size(); ++k){
                const size_t i = ribbons[k].i;
                const size_t j = ribbons[k].j;

                if(1 == j || ySize - 1 == j || xSize - 1 == i){
                    meshIsNotExhausted = false;
                    break;
                }

                const double d = distances[i][j];

                if(0 < i){
                    findNewRibbons(i - 1, j, d, dMin1, dCenter1, Wk, mesh,
                        index, activeElements, ribbons, oldRibbons, distances,
                        velocities, T, Wk[index[i - 1][j]]);

                    findNewRibbons(i - 1, j - 1, d, dMin2, dCenter2, Wk, mesh,
                        index, activeElements, ribbons, oldRibbons, distances,
                        velocities, T, Wk[index[i - 1][j - 1]]);

                    findNewRibbons(i - 1, j + 1, d, dMin2, dCenter2, Wk, mesh,
                        index, activeElements, ribbons, oldRibbons, distances,
                        velocities, T, Wk[index[i - 1][j + 1]]);
                }

                findNewRibbons(i + 1, j, d, dMin1, dCenter1, Wk, mesh, index,
                    activeElements, ribbons, oldRibbons, distances,
                    velocities, T, Wk[index[i + 1][j]]);

                findNewRibbons(i, j - 1, d, dMin1, dCenter1, Wk, mesh, index,
                    activeElements, ribbons, oldRibbons, distances,
                    velocities, T, Wk[index[i][j - 1]]);

                findNewRibbons(i, j + 1, d, dMin1, dCenter1, Wk, mesh, index,
                    activeElements, ribbons, oldRibbons, distances,
                    velocities, T, Wk[index[i][j + 1]]);

                findNewRibbons(i + 1, j - 1, d, dMin2, dCenter2, Wk, mesh,
                    index, activeElements, ribbons, oldRibbons, distances,
                    velocities, T, Wk[index[i + 1][j - 1]]);

                findNewRibbons(i + 1, j + 1, d, dMin2, dCenter2, Wk, mesh,
                    index, activeElements, ribbons, oldRibbons, distances,
                    velocities, T, Wk[index[i + 1][j + 1]]);

                if(dMin2 < d){
                    mesh[i][j].type = CHANNEL;
                }
            }

            for(size_t k = 0; k < ribbons.size(); ++k){
                const size_t i = ribbons[k].i;
                const size_t j = ribbons[k].j;

                if(1 != mesh[i][j].type){
                    ribbons.erase(ribbons.begin() + k);
                    distances[i][j] = 0;
                }
            }

            if(savedSize != activeElements.size() && meshIsNotExhausted){
                setDerivatives(mesh, activeElements);

                buildPartialInfluenceMatrix(influenceMatrix, activeElements,
                    Wk, WkNew, partialInfluenceMatrix, index
                );
            }

            /// Сохраняем параметры трезины

            calculateCrackGeometry(mesh, distances, length, height);

            fracture.push_back(
                std::vector<double>{
                    T,
                    1000 * wn * Wk[index[i00][j00]],
                    pressure[index[i00][j00]],
                    length,
                    height,
                    length / height
                }
            );

            if(saveSteps){
                for(size_t i = 0; i < xSize; ++i){
                    for(size_t j = 0; j < ySize; ++j){
                        openingAtTheStep[i][j] = 1000 * wn * Wk[index[i][j]];
                        pressureAtTheStep[i][j] = pressure[index[i][j]];
                        concentrationAtTheStep[i][j] = concentration[index[i][j]];
                    }
                }
                ai::saveMatrix(ai::string("./Results/Opening/time=")
                    + ai::string(T), openingAtTheStep, true);
                ai::saveMatrix(ai::string("./Results/Pressure/time=")
                    + ai::string(T), pressureAtTheStep, true);
                ai::saveMatrix(ai::string("./Results/Concentration/time=")
                    + ai::string(T), concentrationAtTheStep, true);
            }

            if(0 == regime){
                savedDistances.resize(ribbons.size());

                for(size_t k = 0; k < ribbons.size(); ++k){
                    const size_t i = ribbons[k].i;
                    const size_t j = ribbons[k].j;

                    savedDistances[k] = distances[i][j];
                }

				ai::printVector(savedDistances);
				
                savedTime = T;
				ai::saveMatrix("./distan",distances, true);
				break;
            }
        }
		
        T += dt;
        ++step;

        /// Меняем скорость закачки в соответсвии с планом

        if(T > timeToChangeInjection){
            Q = injection[injectionIndex][2];
            proppantInjection = injection[injectionIndex][3];
            ++injectionIndex;

            if(injection.size() > injectionIndex){
                timeToChangeInjection = injection[injectionIndex][0];
            }else{
                timeToChangeInjection = modelingTime + 10 * dt;
            }

            std::string line("Time = ");
            line += ai::string(injection[injectionIndex - 1][0])
                + ai::string(". Injection changed.");
            ai::printLine(line);
        }



        if(
            0 <= meshScalingCounter
            && length > 2 * initialRadius
            && height > 2 * initialRadius
        ){
            --meshScalingCounter;
            /*
            удвоение шага сетки (cellMesh)

            пересчёт давлений
            */
            ai::printLine(ai::string("Doubling the mesh at time = ") + ai::string(T));
            axMax *= 2.;
            dx *= 2.;
            dy = dx;
            // dt *= 2.;
            x.clear();
            y.clear();
            for(double i = 0; i <= axMax + epsilon; i += dx){
                x.push_back(i);
            }
            for(double i = -1.0*axMax; i <= 1.0*axMax + epsilon; i += dy){
                y.push_back(i);
            }
            for(size_t i = 0; i < xSize; ++i){
                for(size_t j = 0; j < ySize; ++j){
                    mesh[i][j].setCoordinates(x[i], y[j]);
                }
            }
            std::vector<double> openingNew = zeroVectorXY;
            std::vector< std::vector<double> > distancesNew = zeroMatrixXY;
            ribbons.clear();
            activeElements.clear();
            for(int j = 0; 2 * j < ySize - j00; ++j){
                for(int i = 0; 2 * i < xSize; ++i){
                    openingNew[index[i00 + i][j00 + j]] = Wk[index[i00 + 2 * i][j00 + 2 * j]];
                    openingNew[index[i00 + i][j00 - j]] = Wk[index[i00 + 2 * i][j00 - 2 * j]];
                    distancesNew[i00 + i][j00 + j] = distances[i00 + 2 * i][j00 + 2 * j];
                    distancesNew[i00 + i][j00 - j] = distances[i00 + 2 * i][j00 - 2 * j];

                    // if(
                    //     (0 <= 2 * i - 1 && 0 <= j00 + 2 * j && RIBBON == mesh[2 * i - 1][j00 + 2 * j].type)
                    //     || (0 <= 2 * i - 1 && 0 <= j00 + 2 * j - 1 && RIBBON == mesh[2 * i - 1][j00 + 2 * j - 1].type)
                    //     || (0 <= 2 * i && 0 <= j00 + 2 * j - 1 && RIBBON == mesh[2 * i][j00 + 2 * j - 1].type)
                    //     || (0 <= 2 * i && 0 <= j00 + 2 * j && RIBBON == mesh[2 * i][j00 + 2 * j].type)
                    // ){
                    //     ribbons.push_back(Ribbon(i, j00 + j));
                    //     distancesNew[i][j00 + j] = 1000000000000;
                    //     if(0 <= 2 * i - 1 && 0 <= j00 + 2 * j){
                    //         distancesNew[i][j00 + j] = ai::min(distancesNew[i][j00 + j], distances[2 * i - 1][j00 + 2 * j]);
                    //     }
                    //     if(0 <= 2 * i - 1 && 0 <= j00 + 2 * j - 1){
                    //         distancesNew[i][j00 + j] = ai::min(distancesNew[i][j00 + j], distances[2 * i - 1][j00 + 2 * j - 1]);
                    //     }
                    //     if(0 <= 2 * i && 0 <= j00 + 2 * j - 1){
                    //         distancesNew[i][j00 + j] = ai::min(distancesNew[i][j00 + j], distances[2 * i][j00 + 2 * j - 1]);
                    //     }
                    //     if(0 <= 2 * i && 0 <= j00 + 2 * j){
                    //         distancesNew[i][j00 + j] = ai::min(distancesNew[i][j00 + j], distances[2 * i][j00 + 2 * j]);
                    //     }
                    //     distancesNew[i][j00 + j] += 0.25 * sqrt(2.) * dx;
                    // }
                    // if(
                    //     (0 <= 2 * i - 1 && 0 <= j00 - 2 * j && RIBBON == mesh[2 * i - 1][j00 - 2 * j].type)
                    //     || (0 <= 2 * i - 1 && 0 <= j00 - 2 * j + 1 && RIBBON == mesh[2 * i - 1][j00 - 2 * j + 1].type)
                    //     || (0 <= 2 * i && 0 <= j00 - 2 * j + 1 && RIBBON == mesh[2 * i][j00 - 2 * j + 1].type)
                    //     || (0 <= 2 * i && 0 <= j00 - 2 * j && RIBBON == mesh[2 * i][j00 - 2 * j].type)
                    // ){
                    //     ribbons.push_back(Ribbon(i, j00 - j));
                    //     distancesNew[i][j00 - j] = 1000000000000;
                    //     if(0 <= 2 * i - 1 && 0 <= j00 - 2 * j){
                    //         distancesNew[i][j00 - j] = ai::min(distancesNew[i][j00 - j], distances[2 * i - 1][j00 - 2 * j]);
                    //     }
                    //     if(0 <= 2 * i - 1 && 0 <= j00 - 2 * j + 1){
                    //         distancesNew[i][j00 - j] = ai::min(distancesNew[i][j00 - j], distances[2 * i - 1][j00 - 2 * j + 1]);
                    //     }
                    //     if(0 <= 2 * i && 0 <= j00 - 2 * j + 1){
                    //         distancesNew[i][j00 - j] = ai::min(distancesNew[i][j00 - j], distances[2 * i][j00 - 2 * j + 1]);
                    //     }
                    //     if(0 <= 2 * i && 0 <= j00 - 2 * j){
                    //         distancesNew[i][j00 - j] = ai::min(distancesNew[i][j00 - j], distances[2 * i][j00 - 2 * j]);
                    //     }
                    //     distancesNew[i][j00 - j] += 0.25 * sqrt(2.) * dx;
                    // }
                }
            }
            /// CHECK: нет ли продублированных ribbon-элементов?
            for(size_t i = 0; i < xSize; ++i){
                for(size_t j = 0; j < ySize; ++j){
                    if(distancesNew[i][j] > epsilon){
                        ribbons.push_back(Ribbon(i, j));
                        mesh[i][j].type = RIBBON;
                    }else{
                        mesh[i][j].type = OUTSIDE;
                    }
                    if(epsilon < openingNew[index[i][j]] || RIBBON == mesh[i][j].type){
                        activeElements.push_back(std::vector<size_t>{i, j});
                    }
                }
            }

            std::cout << "Active elements: " << activeElements.size() << "."
                << std::endl;

            setDerivatives(mesh, activeElements);

            for(size_t k = 0; k < activeElements.size(); ++k){
                const size_t i = activeElements[k][0];
                const size_t j = activeElements[k][1];


                if(
                    RIBBON != mesh[i][j].type
                    && REGULAR == mesh[i][j].xDerivative
                    && REGULAR == mesh[i][j].yDerivative
                ){
                    mesh[i][j].type = CHANNEL;
                }

                if(
                    RIBBON != mesh[i][j].type
                    && LEFT == mesh[i][j].xDerivative
                    && REGULAR == mesh[i][j].yDerivative
                ){
                    mesh[i][j].type = CHANNEL;
                }

                // CHECK: правильно считать закачки
                // if(OUTSIDE < mesh[i][j].type){
                //     mesh[i][j].activationTime = T0 - dt;
                // }
            }
            initialRadius *= 2;
            std::cout << "Ribbons: " << ribbons.size() << "." << std::endl;

            Wk = openingNew;
            distances = distancesNew;

            // CHECK: правильно считать барьеры
            // m = 0;
            // for(size_t i = 0; i < xSize; ++i){
            //     for(size_t j = 0; j < ySize; ++j){
            //         ai::assignFromVectorByIntervalCondition(sigmaCol[m], y[j],
            //             barriers);
            //         ++m;
            //     }
            // }

            ai::printLine("types");
            for(std::size_t j = 0; j < ySize; ++j){
                for(std::size_t i = 0; i < xSize; ++i){
                    std::cout << mesh[i][j].type;// << " ";
                }
                std::cout << std::endl;
            }
            ai::printLine("types");
        }
    }

    auto finishTime = ai::time();

    if(meshIsNotExhausted){
        if(runningFromGUI){
            std::cout << "Progress: 1.0" << std::endl;
        }else{
            ai::showProgressBar(1.);
        }

        std::cout << std::endl;
    }else{
        if(runningFromGUI){
            std::cout << "Progress: " << (T - T0) / (modelingTime - T0) << std::endl;
        }else{
            ai::showProgressBar((T - T0) / (modelingTime - T0));
        }

        std::cout << std::endl;
        std::cout << "Attention. Mesh is exhausted! "
            << "Cannot continue calculation." << std::endl;
        std::cout << "Finished at time = " << T << "." << std::endl;
    }

    std::cout << "Time used: " << ai::duration(startTime, finishTime, "s")
        << "s" << std::endl;

    std::vector< std::vector<double> > openingAtTheEnd = zeroMatrixXY;
    std::vector< std::vector<double> > pressureAtTheEnd = zeroMatrixXY;
    std::vector< std::vector<double> > concentrationAtTheEnd = zeroMatrixXY;
    for(size_t i = 0; i < xSize; ++i){
        for(size_t j = 0; j < ySize; ++j){
            openingAtTheEnd[i][j] = 1000 * wn * Wk[index[i][j]];
            pressureAtTheEnd[i][j] = pressure[index[i][j]];
            concentrationAtTheEnd[i][j] = concentration[index[i][j]];
        }
    }

    std::cout << "Saving results... ";
    ai::saveMatrix("./Results/opening", openingAtTheEnd, true);
    ai::saveMatrix("./Results/pressure", pressureAtTheEnd, true);
    ai::saveMatrix("./Results/concentration", concentrationAtTheEnd, true);

    calculateCrackGeometry(mesh, distances, length, height);

    fracture.push_back(
        std::vector<double>{
            T,
            1000 * wn * Wk[index[i00][j00]],
            pressure[index[i00][j00]],
            length,
            height,
            length / height
        }
    );
    std::stringstream comment;
    comment << std::right << std::setw(14) << "Time [min]\t"
        << std::right << std::setw(14) << "Opening [mm]\t"
        << std::right << std::setw(14) << "Pressure [Pa]\t"
        << std::right << std::setw(14) << "Length [m]\t"
        << std::right << std::setw(14) << "Height [m]\t"
        << std::right << std::setw(14) << "Aspect raio [n/d]\t";
    ai::saveMatrix("./Results/fracture", fracture, false, comment.str());
    comment.str(std::string());

    std::cout << "OK." << std::endl;

    /// Вычисляем эффективность

    double fluidVolumeOut = 0.;
    double fluidVolumeIn = 0.;
    double proppantVolumeOut = 0.;
    double proppantVolumeIn = 0.;
    for(std::size_t i = 0; i < injection.size(); ++i){
        double volumeOut = 0.;

        if(injection[i][1] < T){
            volumeOut = (injection[i][1] - injection[i][0]) * injection[i][2];
        }else{
            if(injection[i][0] < T){
                volumeOut = (T - injection[i][0]) * injection[i][2];
            }
        }

        fluidVolumeOut += volumeOut;
        proppantVolumeOut += volumeOut * injection[i][3];
    }

    /// Поправка на автомодельное решение

    proppantVolumeOut -= T0 * injection[0][2] * injection[0][3];

    fluidVolumeOut *= 60 * wn / ts;
    proppantVolumeOut *= 60 * wn / ts;

    for(std::size_t i = i00 + 1; i < xSize; ++i){
        for(std::size_t j = 0; j < ySize; ++j){
            fluidVolumeIn += Wk[index[i][j]];
            proppantVolumeIn += concentration[index[i][j]] * Wk[index[i][j]];
        }
    }

    fluidVolumeIn *= 2.;
    proppantVolumeIn *= 2.;

    for(std::size_t j = 0; j < ySize; ++j){
        fluidVolumeIn += Wk[index[i00][j]];
        proppantVolumeIn += concentration[index[i00][j]] * Wk[index[i00][j]];
    }

    fluidVolumeIn *= wn * dx * dy;
    proppantVolumeIn *= wn * dx * dy * 2.65 * 1000;

    std::cout << std::endl;
    std::cout << "Efficiency: " << 100 * fluidVolumeIn / fluidVolumeOut
        << "% [fluid]";
    if(0. < proppantVolumeOut){
        std::cout << ", " << 100 * proppantVolumeIn / proppantVolumeOut
            << "% [proppant]";
    }
    std::cout << "." << std::endl;

    return 0;
}


void SaveJson(const std::string filename,
	std::vector<double> &Wk,
	std::vector<double> &pressure,
	std::vector<double>&concentration,
	std::vector< std::vector<size_t> > &index,
	double Time,
	double xSize,
	double ySize,
	double wn,
	double DproppantInjection						//Величина мнгновенной плотности проппанта (кг/(м3*dt)
	)
{





	//////////////////////////////////////////////////////////////////////////////////////////
	//Эти данные мне непонятно откуда брать
	//////////////////////////////////////////////////////////////////////////////////////////
	double smin = 288.4;								//инимальное значиние напряжений в интервале перфорации, [атм]
	long double fluid_efficiency = 48.43787629043374;	//эффективность жидкости в момент времени "time", [%]
	double net_pressure = 41.505082768865805;			//значение чистого давления на устье трещины в слое интеравала перфорации с напряжением "Smin", [атм]
	double prop_density = 666.666;						// mass density of components - плотность компонент(в той же последовательности, что и в "concentrations"), [кг / куб.м]
	double fluid_density = 999.666;						// mass density of components - плотность компонент(в той же последовательности, что и в "concentrations"), [кг / куб.м]


	//////////////////////////////////////////////////////////////////////////////////////////
	//Эти данные пока МФТИ просит заполнять нулями
	//////////////////////////////////////////////////////////////////////////////////////////
	double azimuth = 0;								//азимут расчетной ячейки в пространстве(задается 0)
	int id = 0;										//(задается 0)
	int stage_id = 0;								//(задается 0)
	int num_fluids = 1;								//number of fluids - общее количество жидкостей
	int num_proppants = 1;							//number of proppants - общее количество пропантов





	std::string extension(".json");
	std::ofstream output(filename+ ai::string(Time) + extension);

	if (!output.good()) {
		throw std::runtime_error(
			ai::string("Exception while saving the matrix into the file: ")
			+ filename
		);
	}


	output << "{\n";
	//////////////////////////////////////////////////////////////////////////////////////////
	//Results
	//////////////////////////////////////////////////////////////////////////////////////////
	output << std::setprecision(14) << "\"Results\": {\n";
	//Smin - минимальное значиние напряжений в интервале перфорации, [атм]
	output << "\t\"Smin\": " << smin << ",\n";
	output << "\t\"accumulated data\": {\n";
	//	fluid efficiency - эффективность жидкости в момент времени "time", [%]
	output << "\t  \"fluid efficiency\": " << fluid_efficiency << ",\n";
	//	net pressure - значение чистого давления на устье трещины в слое интеравала перфорации с напряжением "Smin", [атм]
	output << "\t  \"net pressure\": " << net_pressure << ",\n";
	//	proppant concentration - конентрация пропанта в момент времени "time" на входе в трещину, [кг / куб.м.]
	output << "\t  \"proppant concentration\": " << DproppantInjection << ",\n";
	//	rate - расход смеси в момент времени "time" на входе в трещину, [куб.м / мин]
	output << "\t  \"rate\": " << Q << ",\n";
	//	time - момент времени на который записан результат расчета
	output << "\t  \"time\": " << Time << "\n";
	output << "\t},\n"; //close accumulated data

	//////////////////////////////////////////////////////////////////////////////////////////
	//	geometry
	//////////////////////////////////////////////////////////////////////////////////////////
	output << "\t\"geometry\": {\n";

	//	branches - обозначение полукрыла трещины(0)
	output << "\t\t\"branches\": [\n\t\t {\n";
	//	cells - расчетные ячейки записанные в последовательности : сверху - вниз(j индексы), слева - направо(i индексы)
	output << "\t\t \"cells\": [\n";
	int activeElements = 0;			//Число активных элементов (те элементы, где вектор раскрытия не нулевой)

	//////////////////////////////////////////////////////////////////////////////////////////
	//	Расчет числа активных элементов
	//////////////////////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < xSize; i++)
	{
		for (int j = 0; j < ySize; j++)
		{
			if (Wk[index[i][j]] != 0) activeElements++; // если ячейка не раскрылась, то ее не выводим!!!
		}
	}
	int Index=0;								//Индексная перменная - номер выводимого массива данных

	for (int i = 0; i < xSize; i++)
	{
		for (int j = 0; j < ySize; j++)
		{
			if (Wk[index[i][j]] == 0) continue; // если ячейка не раскрылась, то ее не выводим!!!
			Index++;
			output << "\t\t{ \n";
			//	azimuth - азимут расчетной ячейки в пространстве(задается 0)
			output << "\t\t  \"azimuth\": " << azimuth << ",\n";
			//	concentrations - объемная доля концентрации фаз, [д.ед.]
			output << "\t\t  \"concentrations\": [\n ";
				//for (int index=0; index< num_proppants; index++) // вывод нескольких проппантов (когда будет несколько проппантов)
			output << "\t\t  " << concentration[index[i][j]] << ",\n";			//Вывод концентрации проппанта
			output << "\t\t  " << 1-concentration[index[i][j]] << "\n";						//Вывод концентрации жидкости

			output << "\t\t ],\n"; //close concentrations
								   //	dx, dy, dz - размеры расчетной ячейки(dy - раскрытие), [м]
			output << "\t\t  \"dx\": " << dx << ",\n";
			output << "\t\t  \"dy\": " << Wk[index[i][j]]*wn << ",\n";
			output << "\t\t  \"dz\": " << dy << ",\n";
			//	id - 0; stage id - 0
			output << "\t\t  \"i\": " << i << ",\n";
			output << "\t\t  \"id\": " << id << ",\n";
			output << "\t\t  \"j\": " << j << ",\n";
			output << "\t\t  \"stage id\": " << stage_id << ",\n";
			//	x, y, z - координаты центра расчетной ячеки в пространстве, [м]
			output << "\t\t  \"x\": " << i*dx << ",\n";
			output << "\t\t  \"y\": " << j*dy << ",\n";
			output << "\t\t  \"z\": " << 0 << "\n";
			if (i != xSize - 1 && j != ySize - 1 && activeElements > Index)
				output << "\t\t },\n";
			else
				output << "\t\t }\n";
		}
	}

	output << "\n\t\t]\n"; //close cells
	output << "\t\t }\n\t\t]\n"; //close branches
	output << "\t},\n"; //close geometry

	//////////////////////////////////////////////////////////////////////////////////////////
	//	grid - параметры расчетной области
	//////////////////////////////////////////////////////////////////////////////////////////
	output << "\t\"grid\": {\n";
	//	nx - общее количество расчетных ячеек по Х
	output << "\t\t  \"nx\": " << xSize << ",\n";
	//	nz - общее количество расчетных ячеек по Z
	output << "\t\t  \"nz\": " << ySize << "\n";
	output << "\t},\n"; //close grid

	//////////////////////////////////////////////////////////////////////////////////////////
	//	slurry - свойства смеси в закачке
	//////////////////////////////////////////////////////////////////////////////////////////
	output << "\t\"slurry\": {\n";

	output << "\t\t\"mass density of components\": [\n";
	//	mass density of components - плотность компонент(в той же последовательности, что и в "concentrations"), [кг / куб.м]
	output << "\t\t " << prop_density <<",\n";
	output << "\t\t " << fluid_density <<"\n";

	output << "\t\t ],\n"; //close mass density of components

	output << "\t\t\"name of components\": [\n";
	//	name of components - названия компонент(в той же последовательности, что и в "concentrations")
	output << "\t\t \"Propant 1\",\n";
	output << "\t\t \"Fluid 1\"\n";
	output << "\t\t ],\n"; //close name of componentss

	//	number of fluids - общее количество флюидов
	output << "\t\t  \"number of fluids\": " << num_fluids << ",\n";
	//	number of proppants - общее количество пропантов
	output << "\t\t  \"number of proppants\": " << num_proppants << "\n";
	output << "\t}\n"; //close slurry

	output << "},\n"; //close Results


	//////////////////////////////////////////////////////////////////////////////////////////
	//coordinates - координаты центра интервала перфорации, [м]
	//////////////////////////////////////////////////////////////////////////////////////////
	output << "\"coordinates\": {\n";
	output << "\t\t  \"x\": " << i00 * dx << ",\n";
	output << "\t\t  \"y\": " << j00 * dy << ",\n";
	output << "\t\t  \"z\": " << 0 << "\n";

	output << "}\n";
	output << "}\n";


	output.close();
}
