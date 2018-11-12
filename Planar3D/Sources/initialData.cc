#include <vector>
#include <iostream>

#include "Include/nlohmann/json.hpp"

#include "ai.hh"
#include "planar3D.hh"
#include "linearOperator.hh"

/// \details Функция рассчитывает начальное раскрытие в точке по 
/// автомодельному решению для раскрытия и радиуса трещины
double getInitialOpening(
    const double x,
    const double y,
    const double radius,
    const std::vector<double> &zP,
    const std::vector<double> &W0
){
    const double distance = std::sqrt(std::pow(x, 2) + std::pow(y, 2)) 
        / radius;
    
    double opening = 0;
    
    if(1 > distance){
        for(size_t i = 0; i < zP.size() - 1; ++i){
            if(distance >= zP[i] && distance < zP[i + 1]){
                opening = W0[i] * zP[i + 1] - W0[i + 1] * zP[i]
                    - distance * (W0[i] - W0[i + 1]);
                opening /= (zP[i + 1] - zP[i]);
            }
        }
    }
    
    return opening;
}

/// \details Функция подгружает данные о контрасте напряжений в пласте и 
/// автомодельное решение (радиус трещины, давление и раскрытие) для жидкости
/// заданной реологии из соответствующих текстовых файлов
bool setInitialData(
    double &bn,
    double &radius,
    std::vector< std::vector<double> > &modelSolution,
    const std::string pathToLayersFile,
    std::vector< std::vector<double> > &layers,
    const std::string pathToInjectionFile,
    std::vector< std::vector<double> > &injection
){
    double approximateBn = 1;
    
    const std::vector< std::vector<double> > bnIntervals = {
        {0, 0.2, 0.7296},
        {0.2, 0.4, 0.7236},
        {0.4, 0.6, 0.7162},
        {0.6, 1, 0.6975754},
        {1, 1, 0.6975754}
    };
    
    ai::assignFromVectorByIntervalCondition(
        approximateBn, 
        radius, 
        bn, 
        bnIntervals
    );
    
    /// Загружаем автомодельное решение из файла
    try{
        /// \todo Заменить загрузку файла на вызов функции
        std::cout << "Model solution (" << approximateBn << ")... ";
        
        ai::parseFileInMatrix(
            "./InitialConditions/A" + ai::string(approximateBn) + ".txt", 
            '\t', 
            modelSolution
        );
        
        std::cout << "OK. Size: " << modelSolution.size() << "." << std::endl;
    }catch(const std::exception &error){
        std::cout << "Fail!" << std::endl;
        
        std::cerr << "Cannot open './InitialConditions/A"
            << ai::string(approximateBn) + ".txt', provide file and restart." 
            << std::endl;
        
        std::cerr << error.what() << std::endl;
        
        return false;
    }
    
    if(std::string() != pathToInjectionFile){
        try{
            std::cout << "Injection... ";
            
            ai::parseFileInMatrix(pathToInjectionFile, ' ', injection);
            
            std::cout << "OK. Size: " << injection.size() << "." << std::endl;
        }catch(const std::exception &error){
            std::cout << "Fail!" << std::endl;
            
            std::cerr << "Cannot open '" << pathToInjectionFile << "', "
                << "provide file and restart." << std::endl;
            
            std::cerr << error.what() << std::endl;
            
            return false;
        }
    }else{
        std::cout << "Injection... No file." << std::endl;
    }
    
    if(std::string() != pathToLayersFile){
        try{
            std::cout << "Layers... ";
            
            ai::parseFileInMatrix(pathToLayersFile, ' ', layers);
            
            std::cout << "OK. Size: " << layers.size() << "." << std::endl;
        }catch(const std::exception &error){
            std::cout << "Fail!" << std::endl;
            
            std::cerr << "Cannot open '" << pathToLayersFile << "', "
                << "provide file and restart." << std::endl;
            
            std::cerr << error.what() << std::endl;
            
            return false;
        }
    }else{
        std::cout << "Layers... No file." << std::endl;
    }
    
    return true;
}

bool importInitialData(
    const std::string pathToImportFolder,
    nlohmann::json &importData
){
    std::cout << "Reading initial conditions from text files:" << std::endl;
    
    if(std::string() != pathToImportFolder){
        try{
            std::cout << "Parameters (JSON)... ";
            
            std::string importDataString = std::string();
            
            ai::parseFileIntoString(
                pathToImportFolder + std::string("/parameters.json"), 
                importDataString
            );
            
            importData = nlohmann::json::parse(importDataString);
            
            std::cout << "OK. Size: " << importData.size() << "." << std::endl;
        }catch(const std::exception &error){
            std::cout << "Fail!" << std::endl;
            
            std::cerr << "Cannot open '" << pathToImportFolder 
                << "/parameters.json', provide file and restart." << std::endl;
            
            std::cerr << error.what() << std::endl;
            
            return false;
        }
    }else{
        std::cout << "Parameters... No file." << std::endl;
    }
    
    return true;
}


bool recalculateInjection(
    std::vector< std::vector<double> > &injection,
    const double modelingTime
){
    if(1 > injection.size()){
        std::cerr << "Format error in the injection file: matrix is empty." 
            << std::endl;
        
        return false;
    }
    
    for(std::size_t i = 0; i < injection.size(); ++i){
        if(6 != injection[i].size()){
            std::cerr << "Format error in the injection file: matrix size "
                << "must be four." << std::endl;
            
            return false;
        }
    }
    
    for(size_t i = 0; i + 1 < injection.size(); ++i){
        if(injection[i][1] < injection[i + 1][0]){
            injection.insert(
                injection.begin() + i + 1,
                std::vector<double>{
                    injection[i][1],
                    injection[i + 1][0],
                    0.,
                    0.
                }
            );
        }
    }
    
    /// \todo Значения bn и mu подставлять некорректно
    if(modelingTime > injection[injection.size() - 1][1]){
        injection.push_back(
            std::vector<double>{
                injection[injection.size() - 1][1],
                modelingTime,
                0.,
                0.,
                injection[injection.size() - 1][4],
                injection[injection.size() - 1][5],
            }
        );
    }
    
    return true;
}

bool recalculateStressContrast(
    const std::vector< std::vector<double> > &layers,
    std::vector<double> &stress,
    const std::vector<double> &y
){
    std::vector<double> coordinates;
    std::vector<double> coefficients;
    
    double nominalStress = 0.;
    
    for(int i = layers.size() - 1; i >= 0; --i){
        coordinates.push_back(layers[i][0]);
        coordinates.push_back(layers[i][1] - 0.001);
        coefficients.push_back(layers[i][2]);
        coefficients.push_back(layers[i][2]);
        
        if(0. > layers[i][0] * layers[i][1]){
            nominalStress = layers[i][2];
        }
    }
    
    /// \todo Заменить линейный оператор на линейную интерполяцию
    try{
        createLinearOperator(coordinates, coefficients);
    }catch(const std::exception &error){
        std::cerr << "Error during calculations of the linear operator for "
            << "stress barriers." << std::endl;
        
        std::cerr << error.what() << std::endl;
        
        return false;
    }
    
    for(std::size_t i = 0; i < y.size(); ++i){
        stress.push_back(
            (
                calculateValueWithLinearOperator(
                    y[i],
                    coordinates,
                    coefficients
                ) - nominalStress
            )
        );
    }
    
    return true;
}

/// \details Функция сохраняет начальные параметры в файл формата JSON по 
/// установленному шаблону
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
){
    std::ofstream output(filename + ai::string(".json"));

    if(!output.good()){
        throw std::runtime_error(
            ai::string("Exception while saving the data into the file: ") 
            + filename
        );
    }
    
    std::size_t indentLevel = 0;
    
    auto indent = [&indentLevel, &tab]() -> std::string
    {
        std::string output;
        
        for(std::size_t i = 0; i < indentLevel; ++i){
            output += tab;
        }
        
        return output;
    };
    
    output << "{" << std::endl;
    ++indentLevel;
    
    output << indent() << "\"model\": \"Planar3D\"," << std::endl;
    output << indent() << "\"time\": " << modelingTime << "," << std::endl;
    
    output << indent() << "\"mesh\": {" << std::endl;
    ++indentLevel;
        output << indent() << "\"height\": " << meshHeight << "," << std::endl;
        output << indent() << "\"length\": " << meshLength << "," << std::endl;
        output << indent() << "\"cell\": {" << std::endl;
        ++indentLevel;
            output << indent() << "\"height\": " << cellHeight << "," 
                << std::endl;
            output << indent() << "\"length\": " << cellLength
                << std::endl;
        --indentLevel;
        output << indent() << "}" << std::endl;
    --indentLevel;
    output << indent() << "}," << std::endl;
    
    output << indent() << "\"injection\": {" << std::endl;
    ++indentLevel;
        for(size_t i = 0; i < injection.size(); ++i){
            if(injection[i].size() != 6){
                throw std::runtime_error(
                    ai::string("Exception in size of the vector: injection") 
                );
            }
            output << indent() << "\"" << ai::string(i) << "\": {" 
                << std::endl;
            ++indentLevel;
                output << indent() << "\"start time\": " << injection[i][0]
                    << "," << std::endl;
                output << indent() << "\"stop time\": " << injection[i][1]
                    << "," << std::endl;
                output << indent() << "\"pumping rate\": " 
                    << injection[i][2] << "," << std::endl;
                output << indent() << "\"bulk proppant\": " 
                    << injection[i][3] << "," << std::endl;
                output << indent() << "\"rheology index\": " << injection[i][4]
                    << "," << std::endl;
                output << indent() << "\"dynamic viscosity\": " 
                    << injection[i][5] << std::endl;
            --indentLevel;
            output << indent() << "}";
            if(injection.size() > i + 1){
                output << ",";
            }
            output << std::endl;
        }
    --indentLevel;
    output << indent() << "}," << std::endl;
    
    output << indent() << "\"layers\": {" << std::endl;
    ++indentLevel;
        for(size_t i = 0; i < layers.size(); ++i){
            if(layers[i].size() != 6){
                throw std::runtime_error(
                    ai::string("Exception in size of the vector: layers") 
                );
            }
            
            output << indent() << "\"" << ai::string(i) << "\": {" 
                << std::endl;
            ++indentLevel;
                output << indent() << "\"y1\": " << layers[i][1] << ","
                    << std::endl;
                output << indent() << "\"y2\": " << layers[i][0] << ","
                    << std::endl;
                output << indent() << "\"stress\": " << layers[i][2]
                    << "," << std::endl;
                output << indent() << "\"Young's modulus\": "
                    << layers[i][3] << "," << std::endl;
                output << indent() << "\"Poisson's ratio\": " 
                    << layers[i][4] << "," << std::endl;
                output << indent() << "\"Carter's coefficient\": " 
                    << layers[i][5] << std::endl;
            --indentLevel;
            output << indent() << "}";
            if(layers.size() > i + 1){
                output << ",";
            }
            output << std::endl;
        }
    --indentLevel;
    output << indent() << "}" << std::endl;
    
    output << "}";
}
