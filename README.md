# Occultation light curves simulator

# Installation
1. Install FFTW from this repository.
```
  cd external/fftw-3.3.10
  configure
  make
  make install
```

2. Install GNU Scientific Library (GSL).

For example, on Arch Linux, you would run
```
sudo pacman -S gsl
```

3. Download the program from GitHub Releases, so you don't need to compile the repository; it just works.

# Contribution
Requirements cmake, a C compiler, GSL and FFTW.

1. Clone the repository
git clone 

2. Cmake
```
cmake -B build && cmake --build build && cd build/bin/
```

