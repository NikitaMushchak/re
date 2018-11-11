#pragma once

/// \brief Функция заполняет матрицу коэффициентов влияния
void buildInfluenceMatrix(
    std::vector< std::vector<double> > &influenceMatrix,
    const size_t xSize,
    const size_t ySize
);
