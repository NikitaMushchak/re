# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

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