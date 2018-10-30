# Planar3D

Fracture propagation model

## Getting Started

### Prerequisites

Please, install this software to begin with

```
- C++11 Compiler (ISO/IEC 14882:2011)
- Python 2.7
- PyQt 4
- PyInstaller
```

### Building

Building Planar3D with clang++

```
clang++ -std=c++11 -O3 -Wall -pedantic -o ./planar3D.exe ./Sources/Sequential/*.cc
```

Building Planad3D GUI with pyinstaller

```
pyinstaller --windowed ./Sources/Interface/*.pyw
```

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/starobinskii/Planar3D/tags). 

## Authors

* **Alexey Stepanov** - *Initial work*
* **Egor Starobinskii** - *C++ program*