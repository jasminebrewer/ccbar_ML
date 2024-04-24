
# pythia 8 path
PYTHIA8_INCLUDE=/mnt/users/brewerj/pythia8311/include
PYTHIA8_LIB=/mnt/users/brewerj/pythia8311/lib

# fastjet compilation
FASTJET_INCLUDE=/mnt/users/brewerj/fastjet-install/include
FASTJET_LIB=/mnt/users/brewerj/fastjet-install/lib

g++ ${1}.cc functions_ccbar.cc -o ${1} -w  -I$PYTHIA8_INCLUDE -O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread  -L$PYTHIA8_LIB -Wl,-rpath,$PYTHIA8_LIB -lpythia8 -ldl -I$FASTJET_INCLUDE -L$FASTJET_LIB -Wl,-rpath,$FASTJET_LIB -lfastjet

