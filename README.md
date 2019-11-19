# Phase-change

## Build Instructions

This project is based on [CMake](https://cmake.org/). Simply generate project, Makefiles, etc. using [CMake](https://cmake.org/) and compile the project with the compiler of your choice. The code was tested with the following configurations:
- Windows 10 64-bit, CMake 3.9.5, superior version Visual Studio 2015
- Ubuntu 16.10 64-bit, CMake 3.5.2, GCC 6.2.0.
- Cuda 9.2, nvvm64_32_0

### build in windows
    mkdir build
    cd build
    MSBuild.exe PhysAni.sln (/property:Configuration=Release)
    ../bin/Demo1

### build in linux
    mkdir build
    cd build
    make PhysAni
    ../bin/Demo1

Note: Please use a 64-bit target on a 64-bit operating system. 32-bit builds on a 64-bit OS are not supported.

## Documentation

## Latest Important Changes


## Features



## Videos


## Screenshots
		

## References

* J. Bender, M. Müller and M. Macklin, "Position-Based Simulation Methods in Computer Graphics", In Tutorial Proceedings of Eurographics, 2015
* J. Bender, M. Müller, M. A. Otaduy, M. Teschner and M. Macklin, "A Survey on Position-Based Simulation Methods in Computer Graphics", Computer Graphics Forum 33, 6, 2014
* M. Macklin, M. Müller, N. Chentanez and T.Y. Kim, "Unified particle physics for real-time applications", ACM Trans. Graph. 33, 4, 2014
* M. Macklin and M. Müller, "Position based fluids", ACM Trans. Graph. 32, 4, 2013
