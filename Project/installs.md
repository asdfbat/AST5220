* Download tar from ftp://ftp.gnu.org/gnu/gsl/

* Go into the directory and run
./configure --prefix=$HOME/local
make
make install

* Go into the Makefile and edit the INC and LIBS lines to:
INC  = -I$(HOME)/local/include
LIBS = -L$(HOME)/local/lib -lgsl -lgslcblas

* For now, comment out the following lines in the Makefile:
#OPTIONS += -D_COMPLEX_BESSEL
#LIBS += -lgfortran -lcomplex_bessel

* Add the following line to ~/.bashrc
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/lib