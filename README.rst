================================
ecell4
================================

About
=====

E-Cell System, a software platform for modeling, simulation and analysis of complex, heterogeneous and multi-scale systems like the cell.

Install
=======

Do following instructions.

..

  $ ./waf update --files="boost,doxygen" configure --prefix=${PREFIX} build install

Full installation (Read only)

..

  $ cd ${SRCPATH}
  $ git clone git://github.com/ecell/ecell4.git ecell4
  $ git clone git://github.com/headmyshoulder/odeint-v2 odeint-v2
  $ git clone git://github.com/ecell/epdp.git epdp
  $ cd epdp
  $ ./autogen.sh && ./configure && make -j
  $ ln -s _gfrd.so libgfrd.so
  $ cd ../ecell4
  $ PREFIX=/foo/bar LIBRARY_PATH=${SRCPATH}/epdp LD_LIBRARY_PATH=${SRCPATH}/epdp \
  CPLUS_INCLUDE_PATH=${SRCPATH}/epdp:${SRCPATH}/odeint-v2 \
  ./install.sh core bd gillespie ode egfrd
