# TITULO

This repository contains the code 
for the XXXXX
for the [paper]() [1]. 

**Authors:** [Mercedes Garcia-Salguero](https://mapir.isa.uma.es/mapirwebsite/?p=1718), 
[Javier Gonzalez-Jimenez](https://mapir.isa.uma.es/mapirwebsite/?p=1536)


**License:** [GPLv3](https://raw.githubusercontent.com/mergarsal/CertHPrior/main/LICENSE)


If you use this code for your research, please cite:

```
xxx
```

The generic certifier can be found in 
```
https://github.com/mergarsal/QCQPIterCertifier
```

## Dependences 
* Eigen 
 ```
        sudo apt install libeigen3-dev
 ```

* Optimization (own fork)
 ```
        git clone https://github.com/mergarsal/Optimization.git
 ```
* Iterative certifier
```
    git clone https://github.com/mergarsal/QCQPIterCertifier.git
```


## Build
```
git clone --recursive https://github.com/mergarsal/CertHPrior.git
cd GNCSO

mkdir build & cd build 

cmake .. 

make -jX

```

Exmaples are found in `bin` folder. For example: 
```
        ./bin/example_basic_plane
```
 







