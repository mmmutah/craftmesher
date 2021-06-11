# CRack And deFecT Mesher
The crack and defect mesher is a tool written in C++ that is able to generate gradated tetrahedral meshes of polycrystals informed by the cracks and defects in the model.

If you use this tool or its key elements in your publication, we kindly ask you cite the following publication:

```
B.R. Phung, J. He, A.D. Spear. "A surface-mesh gradation tool for generating gradated tetrahedral meshes of microstructures with defects" Computational Materials Science XXX (2021) 110622.
```

# INSTALL INSTRUCTIONS

## Prerequisites
* Boost
* Pytrhon
* TetGen

## Linux: autotools
```
aclocal  
autoconf  
automake --add-missing  
./configure  
make  
```

## macOS with Clang (Intel and AS)

### Compile and build boost libraries  
* Download https://www.boost.org/users/download/#live  
* Unpack the archive  
* cd into directory  
* ./bootstrap.sh --prefix=/some/directory/to/install/boost/to  
* ./b2  
* ./b2 install  

### Make
* Set boost inc and lib directories in Makefile_macOS  
* Rename Makefile_macOS to Makefile  
* Run make  

