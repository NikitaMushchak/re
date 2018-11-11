# Список изменений
Все важные изменения в проекте будут записываться в этот файл.

Используется формат [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) и [Семантическое Версионирование](http://semver.org/spec/v2.0.0.html).

## [3.1.0] - 2018-11-04
### Изменено
- Будьте внимательны! Теперь флаг `--cell-size` задаёт размер ячейки в метрах!
- Будьте внимательны! Теперь флаг `--mesh-size` задаёт размер сетки в ячейках!
- Будьте внимательны! Теперь файл `parameters.json` имеет другую разметку и единицы измерения. Подробнее в файле `README.md`.
- Изменили коды некоторых ошибок.
- Вернули сохранение давления в трещине.
- Теперь в файл `fracture_m.txt` пишется также призабойное раскрытие.
- Теперь в файл `fracture_m.txt` пишется полная длина трещины, а не полудлина
- Исправили расчёт высоты трещины для несимметричных трещин.

## [3.0.0] - 2018-10-26
### Added
- Makefile to build on Windows, macOS and Linux.
- Proppant calculation.
- Parallel execution option using OpenMP. Use flag `--threads` to specify number of parallel threads.
- Use flag `--import` to read parameters from previous calculations.
- Use flag `--layers` to read data about layers.

### Changed
- Be careful! Now flag `--Q=<value>` sets injection in m^3/min!
- Incoming text files format.
- Some function declarations and folder structure.
- Save opening and proppant concentration transposed.
- Don't save pressure except the injection point.

### Deleted
- Flag `--barriers` is now deprecated. Use `--layers`.

## [2.1.1] - 2018-07-10
### Changed 
- Matrix-vector multiplying is faster.
- Compare `<Cell>.type`, `<Cell>.xDerivative` and `<Cell>.yDerivative` with `int` instead of `std::string`.
- Fix calculation of stepToCheck for a leak-off dominated regime.

## [2.1.0] - 2018-06-24
### Added
- Use flag `--injection` to specify injection change with a txt-file.
- Use flag `--viscosity` to set viscosity dominated regime.
- Use flag `--toughness` to set toughness dominated regime.
- Use flag `--leak-off` to set leak-off dominated regime.
- Save information about injection into a txt-file.
- Provide more information about errors while parsing incoming files.

### Changed
- Use `+` as a separator in a build version instead of `~`.
- Use `\t` as a separator in incoming txt-files for injection and barriers instead of ` `.
- Add comments into the txt-file with information about barriers.
- Update AiLibrary to 1.1.0.