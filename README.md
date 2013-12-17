StencilParareal
===============

To build the project use cmake, preferably in a different build directory.
You will have to build the stencil library with "make install" in its build directory.
You will also have to provide the path of the stencil library to cmake by running

cmake -DSTENCIL_LIBRARY_PATH=${path}

where ${path} is the directory containing the "build" directory, the "src" directory and the README file of the stencil library.
