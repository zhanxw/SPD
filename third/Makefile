# download all third party library and compile them

LIB = tabix pcre eigen gsl bzip2 samtools
all: $(LIB)

tabix: tabix-0.2.6.tar.bz2
	-mv tabix-0.2.6.tar.bz2* tabix-0.2.6.tar.bz2
	tar jvxf tabix-0.2.6.tar.bz2
	ln -s -f tabix-0.2.6 tabix
	sed -i -e 's:-D_USE_KNETFILE::' tabix/Makefile
	(cd tabix-0.2.6; make)

pcre: pcre-8.21.tar.gz
	-mv pcre-8.21.tar.gz* pcre-8.21.tar.gz
	tar zvxf pcre-8.21.tar.gz
	-(DIR=`pwd`; cd pcre-8.21; ./configure --prefix="$${DIR}"/pcre; make -j; make install)

gsl: gsl-1.15.tar.gz
	tar zvxf gsl-1.15.tar.gz
	-(DIR=`pwd`; cd gsl-1.15; ./configure --prefix="$${DIR}"/gsl; make -j; make install)

eigen: eigen-3.1.4.tar.bz2
	-rm -rf eigen-eigen*
	tar jvxf eigen-3.1.4.tar.bz2
	-mv eigen-eigen-* eigen-3.1.4
	ln -s -f eigen-3.1.4 eigen

bzip2: bzip2-1.0.6.tar.gz
	tar zvxf $<
	ln -s -f bzip2-1.0.6 bzip2
	(cd bzip2-1.0.6; make)

samtools: samtools-0.1.19.tar.bz2
	tar jvxf $<
	ln -s -f samtools-0.1.19 samtools
	(cd samtools; make)

######################################################################
# download sources
eigen-3.1.4.tar.bz2:
	echo "obtain Eigen..."
	wget --no-check-certificate http://bitbucket.org/eigen/eigen/get/3.1.4.tar.bz2 -O $@

gsl-1.15.tar.gz:
	echo "obtain GSL"
	wget http://mirrors.kernel.org/gnu/gsl/gsl-1.15.tar.gz

bzip2-1.0.6.tar.gz:
	echo "obtain Bzip2 library"
	wget http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
clean:
	-rm -rf $(LIB)
deepclean:
	-rm -rf $(LIB)
	-rm -rf tabix-0.2.6
	-rm -rf pcre-8.21
	-rm -rf gsl-1.15
	-rm -rf eigen-3.1.4
	-rm -rf bzip2-1.0.6
	-rm -rf samtools-0.1.19
