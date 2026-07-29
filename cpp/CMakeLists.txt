cmake_minimum_required(VERSION 3.18)
project(polybuilder CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# RDKit ships a CMake package config. If CMake cannot find it, point it at
# your install, e.g.:  cmake -DRDKit_DIR=$CONDA_PREFIX/lib/cmake/rdkit ..
find_package(RDKit REQUIRED)
find_package(Boost REQUIRED)

add_executable(polybuilder
  main.cpp
  polymer_builder.cpp
)

target_include_directories(polybuilder PRIVATE
  ${RDKit_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

# Link the RDKit component libraries the engine uses.
target_link_libraries(polybuilder PRIVATE
  RDKitGraphMol
  RDKitSmilesParse
  RDKitSubstructMatch
  RDKitChemTransforms
  RDKitFileParsers
  RDKitRDGeneral
  RDKitDataStructs
)
