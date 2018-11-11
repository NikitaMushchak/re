#pragma once

/// \brief Функция копирует фрагмент из матрицы коэффициентов влияния
void buildPartialInfluenceMatrix(
    std::vector< std::vector<double> > &newInfluenceMatrix,
    std::vector< std::vector<size_t> > &activeElements,
    std::vector<double> &Wk,
    std::vector<double> &WkNew,
    std::vector< std::vector<double> > &partialInfluenceMatrix,
    std::vector< std::vector<size_t> > &index
);
