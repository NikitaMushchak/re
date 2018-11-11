#pragma once

#include "planar3D.hh"

/// \brief Функция строит профиль раскрытия по автомодельному решению
double getInitialOpening(
     const double x,
     const double y,
     const double radius,
     const std::vector<double> &zP,
     const std::vector<double> &W0
);

/// \brief Функция считывает начальные данные из файлов
bool setInitialData(
    double &bn,
    double &radius,
    std::vector< std::vector<double> > &modelSolution,
    const std::string pathToLayersFile,
    std::vector< std::vector<double> > &layers,
    const std::string pathToInjectionFile,
    std::vector< std::vector<double> > &injection
);

/// \brief Функция загружает начальные параметры из папки с сохранёнными
/// предыдущими расчётами
bool importInitialData(
    const std::string pathToImportFolder,
    nlohmann::json &importData
);

/// \brief Функция пересчитывает параметры закачки в требуемый формат
bool recalculateInjection(
    std::vector< std::vector<double> > &injection,
    const double modelingTime
);

/// \brief Функция пересчитывает контраст напряжений в требуемый формат 
/// с помощью линейного оператора
bool recalculateStressContrast(
    const std::vector< std::vector<double> > &layers,
    std::vector<double> &stress,
    const std::vector<double> &y
);

/// \brief Функция сохранят начальные параметры в файл формата JSON
void saveInitialData(
    const std::string filename,
    const double modelingTime,
    const std::size_t meshHeight,
    const std::size_t meshLength,
    const double cellHeight,
    const double cellLength,
    std::vector< std::vector<double> > &injection,
    std::vector< std::vector<double> > &layers,
    const std::string tab = std::string("    ")
);
