#pragma once

/// \brief Функция определяет список активных элементов
std::vector< std::vector<size_t> > findActiveElements(
    std::vector< std::vector<Cell> > &mesh,
    double testDistance
);
