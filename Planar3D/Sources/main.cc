#include <vector>
#include <iostream>

#include "ai.hh"
#include "planar3D.hh"

/// \todo Исправить учёт утечек при масштабировании сетки

/// \brief Версия программы
/// \details Текущая версия программы
std::string version("3.1.0");

#if defined(_MSC_VER)
    std::string compiler(ai::string("msc") + ai::string(_MSC_VER));
#elif defined(_COMPILER_)
    std::string compiler(TO_STRING(_COMPILER_));
#else
    std::string compiler("%compiler%");
#endif
#if defined(_TARGET_)
    std::string target(TO_STRING(_TARGET_));
#else
    std::string target("%target%");
#endif
#if defined(_TIMESTAMP_)
    std::string timestamp(TO_STRING(_TIMESTAMP_));
#else
    std::string timestamp("%timestamp%");
#endif
#if defined(_PLATFORM_)
    std::string buildVersion(
        version + ai::string("+") + TO_STRING(_PLATFORM_)
    );
#else
    std::string buildVersion(version + ai::string("+%platform%"));
#endif

#if defined(OPENMP)
    #include <omp.h>
#else
    std::size_t omp_get_max_threads(){return 1;}
#endif

/// \brief Main-функция
/// \details Первичная функция, считывающая параметры, переданные при запуске
/// программы, и запускающая расчёт распространения трещины по планарной
/// модели
int main(const int argc, const char *argv[]){
    bool saveSteps = false;
    bool runningFromGUI = false;

    std::string pathToLayersFile("./InitialConditions/layers.txt");
    std::string pathToInjectionFile("./InitialConditions/injection.txt");
    std::string pathToImportFolder = std::string();

    std::size_t numberOfThreads = (std::size_t) omp_get_max_threads();

    int meshScalingCounter = -1;

    double cellSize = 4.;
    double meshSize = 30.;
    double time = 25.0;
    double ts = 60.;

    E = 32 * pow(10, 9);
    epsilon = std::pow(10., -8);
    C = 0. * std::pow(10., -6);
    Kic = 0. * std::pow(10., 6);
    Q = 2.2 / 60.;//4.4/60.
    bn = 1.;
    mu = 0.4;

    /// \todo Добавить define для режимов
    regime = 1;

    for(int i = 1; i < argc; ++i){
        if(
            std::string("-v") == std::string(argv[i])
            || std::string("--version") == std::string(argv[i])
        ){
            std::cout << "  Build: "  << buildVersion << "." << std::endl;
            std::cout << "  Compiler: "  << compiler << "." << std::endl;
            std::cout << "  Target: "  << target << "." << std::endl;
            std::cout << "  AiLibrary: " << ai::getVersion() << "."
                << std::endl;
            std::cout << "  Compilation timestamp: "  << timestamp << "."
                << std::endl;

            return 0;
        }

        if(std::string("--name") == std::string(argv[i])){
            std::cout << "Planar3D-LL" << std::endl;

            return 0;
        }

        if(
            std::string("-h") == std::string(argv[i])
            || std::string("--help") == std::string(argv[i])
        ){
            std::cout << "usage: planar3D [options]"
                << std::endl
                << "    -h  --help            print this usage and exit"
                << std::endl
                << "    -v  --version         print build info and exit"
                << std::endl
                << "    --name                print program name and exit"
                << std::endl
                << "    --list-errors         print possible errors ans exit"
                << std::endl << std::endl

                << "  Liquid parameters" << std::endl
                << "    --Q=<value>           injection [double, m^3 / min]"
                << std::endl
                << "    --mu=<value>          viscosity [double, Pa * s]"
                << std::endl
                << "    --n=<value>           index of the liquid's rheology "
                << "[double, n/d]"
                << std::endl << std::endl

                << "  Shelf parameters" << std::endl
                << "    --E=<value>           flat Young's modulus "
                << "[double, GPa]"
                << std::endl
                << "    --C=<value>           Carter coefficient [double, "
                << "um / s^0.5]"
                << std::endl
                << "    --Kic=<value>         stress intensity factor "
                << "[double, MPa/m^0.5]"
                << std::endl << std::endl

                << "  Time parameters" << std::endl
                << "    --time=<value>        modeling time [double, min]"
                << std::endl
                << "    --ts=<value>          time scale [double, s]"
                << std::endl << std::endl

                << "  Mesh parameters" << std::endl
                << "    --mesh-size=<value>    mesh size in cells [uint, n/d]" 
                << std::endl
                << "    --cell-size=<value>    cell size in meters [double, m]"
                << std::endl
                << "    --mesh-scale=<value> how many times mesh scaling is "
                << "allowed [uint, n/d]"
                << std::endl
                << "    --mesh-scale          allow infinite mesh scaling"
                << std::endl
                << std::endl << std::endl

                << "  Fracture regime" << std::endl
                << "    --viscosity           viscosity dominated regime "
                << "{default}"
                << std::endl
                << "    --toughness           toughness dominated regime "
                << std::endl
                << "    --leak-off            leak-off dominated regime "
                << std::endl << std::endl

                << "  Compressive stress barriers" << std::endl
                << "    --layers=<path>       path to a txt-file with layers "
                << "description [string]"
                << std::endl
                << "    --layers              default load from "
                << "./InitialConditions/layers.txt"
                << std::endl << std::endl

                << "  Dynamic Injection" << std::endl
                << "    --injection=<path>     path to a txt-file with a plan "
                << "of injections [string]"
                << std::endl
                << "    --injection            default load from "
                << "./InitialConditions/injection.txt"
                << std::endl << std::endl

                << "  Other flags" << std::endl
                << "    --import=<path>       import data from folder [string]"
                << std::endl
                << "    --import              default import from ./Results"
                << std::endl
                << "    --save-steps          save concentration and opening "
                << "at every check"
                << std::endl
                << "    --env=gui             flag for FractureGUI"
                << std::endl
                << "    --threads=<value>     number of parallel openmp"
                << " threads [uint, n/d]"
                << std::endl;

            return 0;
        }

        if("--list-errors" == std::string(argv[i])){
            std::cout << "  User input errors" << std::endl
                << "    Code 11. Cell size is less than 1 meter."
                << std::endl
                << "    Code 12. Mesh size is less than 10 cells."
                << std::endl
                << "    Code 13. Mesh size is not a positive integer."
                << std::endl
                << "    Code 14. Carter coefficient isn't positive in a "
                << "leak-off dominated regime."
                << std::endl
                << "    Code 15. Stress intensity factor isn't positive in a "
                << "toughness dominated "
                << std::endl
                << "             regime."
                << std::endl
                << "    Code 16. Viscosity isn't positive in a "
                << "viscosity dominated regime."
                << std::endl << std::endl

                << "  Enviroment errors" << std::endl
                << "    Code 21. Cannot find %folderName%."
                << std::endl
                << "    Code 22. Cannot open %fileName%."
                << std::endl
                << "    Code 23. Format error in %fileName%."
                << std::endl << std::endl

                << "  Calculating errors" << std::endl
                << "    Code 31. Calculations failed due to unmet "
                << "mathematical condition."
                << std::endl
                << "    Code 32. Automodel solution cannot be applied. "
                << std::endl;

            return 0;
        }

        if(
            ai::assignAbsDoubleParameter(argv[i], "--Q=", Q)
            || ai::assignAbsDoubleParameter(argv[i], "--n=", bn)
            || ai::assignAbsDoubleParameter(argv[i], "--mu=", mu)
            || ai::assignAbsDoubleParameter(argv[i], "--time=", time)
            || ai::assignAbsDoubleParameter(argv[i], "--ts=", ts)
            || ai::assignAbsDoubleParameter(argv[i], "--cell-size=", cellSize)
            || ai::assignAbsDoubleParameter(argv[i], "--mesh-size=", meshSize)
            || ai::assignStringParameter(
                argv[i],
                "--layers=",
                pathToLayersFile
            )
            || ai::assignStringParameter(
                argv[i],
                "--injection=",
                pathToInjectionFile
            )
            || ai::assignParameter(
                argv[i],
                "--mesh-scaling=",
                meshScalingCounter
            )
            || ai::assignParameter(argv[i], "--threads=", numberOfThreads)
        ){
            continue;
        }

        if(ai::assignAbsDoubleParameter(argv[i], "--E=", E)){
            E *= std::pow(10., 9);

            continue;
        }

        if(ai::assignAbsDoubleParameter(argv[i], "--C=", C)){
            C *= std::pow(10., -6);

            continue;
        }

        if(ai::assignAbsDoubleParameter(argv[i], "--Kic=", Kic)){
            Kic *= std::pow(10., 6);

            continue;
        }

        if("--viscosity" == std::string(argv[i])){
            regime = 1;

            continue;
        }

        if("--toughness" == std::string(argv[i])){
            regime = 0;

            continue;
        }

        if("--leak-off" == std::string(argv[i])){
            regime = -1;

            continue;
        }

        if("--layers" == std::string(argv[i])){
            pathToLayersFile = std::string("./InitialConditions/layers.txt");

            continue;
        }

        if("--injection" == std::string(argv[i])){
            pathToInjectionFile = std::string(
                "./InitialConditions/injection.txt"
            );

            continue;
        }

        if("--mesh-scale" == std::string(argv[i])){
            meshScalingCounter = 1000;

            continue;
        }

        if("--import" == std::string(argv[i])){
            pathToImportFolder = "./Results";

            continue;
        }

        if("--save-steps" == std::string(argv[i])){
            saveSteps = true;

            continue;
        }

        if("--env=gui" == std::string(argv[i])){
            runningFromGUI = true;

            continue;
        }
    }


    if(1. > cellSize){
        std::cerr << "Cell size should be at least 1 meter" << std::endl;

        return 11;
    }

    if(10. > meshSize){
        std::cerr << "Mesh size should be at least 10 cells." << std::endl;

        return 12;
    }

    if(meshSize != (int) meshSize){
        std::cerr << "Mesh size should be a positive integer." << std::endl;

        return 13;
    }

    if(-1 == regime && epsilon > C){
        std::cerr << "Carter coefficient should be positive in a leak-off "
            << "dominated regime."
            << std::endl;

        return 14;
    }

    if(0 == regime && epsilon > Kic){
        std::cerr << "Stress intensity factor should be positive in a "
            << "toughness dominated regime."
            << std::endl;

        return 15;
    }

    if(1 == regime && epsilon > mu){
        std::cerr << "Viscosity should be positive in a viscosity dominated"
            << "regime."
            << std::endl;

        return 16;
    }
    
    if(0 != (int) meshSize % 2){
        meshSize = round(meshSize + 1);
    }

    return planar3D(
        time, ts, cellSize, meshSize, meshScalingCounter,
        pathToLayersFile, pathToInjectionFile, pathToImportFolder,
        runningFromGUI, saveSteps, numberOfThreads
    );
}
