###############
Installation
###############

* Linux systems:
  Yade can be installed from :ref:`packages <provided-packages>` (pre-compiled binaries) or :ref:`source code <install-from-source-code>`. The choice depends on what you need: if you don't plan to modify Yade itself, package installation is easier. In the contrary case, you must download and install the source code.

* Other Operating Systems:
  Emulating Linux systems including Yade is proposed in this case, through :ref:`docker images <docker-install>` as well `flash-drive or virtual machines <https://yade-dem.org/doc/installation.html#yubuntu>`_ images.

* 64 bit Operating Systems required; no support for 32 bit (i386).

.. toctree::
  :maxdepth: 2

.. _provided-packages:

Packages
----------

**Stable packages**


Since 2011, all Ubuntu (starting from 11.10, Oneiric; and with the exception of Ubuntu 24.04 noble which requires using either daily packages or source code, see below) and Debian (starting from Wheezy) versions
have Yade in their main repositories. There are only stable releases in place.
To install Yade, run the following::

	sudo apt-get install yade

After that you can normally start Yade using the command ``yade`` or ``yade-batch``.

`This image <https://repology.org/badge/vertical-allrepos/yade.svg>`_ shows versions and up to date status of Yade in some repositories.

.. only:: html

	.. image:: https://repology.org/badge/vertical-allrepos/yade.svg
		:alt: Yade versions in different repositories

**Daily packages**

Pre-built packages updated more frequently than the stable versions are provided for all currently supported Debian and Ubuntu
versions and available on `yade-dem.org/packages <http://yade-dem.org/packages/>`_ .

These are "daily" versions of the packages which are being updated regularly and, hence, include
all the newly added features.

To install the daily-version you need to add the repository to your
/etc/apt/sources.list.


- Debian 11 **bullseye**::

	sudo bash -c 'echo "deb http://www.yade-dem.org/packages/ bullseye main" >> /etc/apt/sources.list.d/yadedaily.list'


- Debian 12 **bookworm** also with :ref:`high precision<highPrecisionReal>` ``long double``, ``float128`` and ``mpfr150`` packages::

	sudo bash -c 'echo "deb http://www.yade-dem.org/packages/ bookworm main" >> /etc/apt/sources.list.d/yadedaily.list'


- Debian 13 **trixie** also with :ref:`high precision<highPrecisionReal>` ``long double``, ``float128`` and ``mpfr150`` packages::

	sudo bash -c 'echo "deb http://www.yade-dem.org/packages/ trixie main" >> /etc/apt/sources.list.d/yadedaily.list'


- Ubuntu 18.04 **bionic**::

	sudo bash -c 'echo "deb http://www.yade-dem.org/packages/ bionic main" >> /etc/apt/sources.list.d/yadedaily.list'


- Ubuntu 20.04 **focal**::

	sudo bash -c 'echo "deb http://www.yade-dem.org/packages/ focal main" >> /etc/apt/sources.list.d/yadedaily.list'

- Ubuntu 22.04 **jammy** also with :ref:`high precision<highPrecisionReal>` ``long double``, ``float128`` and ``mpfr150`` packages::

	sudo bash -c 'echo "deb http://www.yade-dem.org/packages/ jammy main" >> /etc/apt/sources.list.d/yadedaily.list'

- Ubuntu 24.04 **noble**::

	sudo bash -c 'echo "deb http://www.yade-dem.org/packages/ noble main" >> /etc/apt/sources.list.d/yadedaily.list'


Add the PGP-key AA915EEB as trusted and install ``yadedaily``::

	wget -O - http://www.yade-dem.org/packages/yadedev_pub.gpg | sudo tee /etc/apt/trusted.gpg.d/yadedaily.asc
	sudo apt-get update
	sudo apt-get install yadedaily


After that you can normally start Yade using the command ``yadedaily`` or ``yadedaily-batch``.
``yadedaily`` on older distributions can have some disabled features due to older library
versions, shipped with particular distribution.


The Git-repository for packaging stuff is available on `GitLab <https://gitlab.com/yade-dev/trunk/tree/master/scripts/ppa_ci>`_.


If you do not need ``yadedaily``-package anymore, just remove the
corresponding line in /etc/apt/sources.list and the package itself::

	sudo apt-get remove yadedaily

To remove our key from keyring, execute the following command::

	sudo apt-key remove AA915EEB

Daily and stable Yade versions can coexist without any conflicts, i.e., you can use ``yade`` and ``yadedaily``
at the same time.

.. _docker-install:

Docker
----------

Yade can be installed using docker images, which are daily built.
Images contain both stable and daily versions of packages, except for Ubuntu 24.04 (see below).
Docker images are based on supported distributions:


- Debian 11 **bullseye**::

	docker run -it registry.gitlab.com/yade-dev/docker-prod:debian-bullseye

- Debian 12 **bookworm**::

	docker run -it registry.gitlab.com/yade-dev/docker-prod:debian-bookworm

- Debian 13 **trixie**::

	docker run -it registry.gitlab.com/yade-dev/docker-prod:debian-trixie

- Ubuntu 18.04 **bionic**::

	docker run -it registry.gitlab.com/yade-dev/docker-prod:ubuntu18.04


- Ubuntu 20.04 **focal**::

	docker run -it registry.gitlab.com/yade-dev/docker-prod:ubuntu20.04

- Ubuntu 22.04 **jammy**::

	docker run -it registry.gitlab.com/yade-dev/docker-prod:ubuntu22.04

- Ubuntu 24.04 **noble**:: (with only the `yadedaily` version)

	docker run -it registry.gitlab.com/yade-dev/docker-prod:ubuntu24.04


After the container is pulled and is running, Yade functionality can be checked::

	yade --test
	yade --check
	yadedaily --test
	yadedaily --check

.. _install-from-source-code:

Source code
------------

Installation from source code is reasonable, when you want to add or
modify constitutive laws, engines, functions etc., or use the recently added features, which are not yet
available in packaged versions.

Doing so, we recommend to separate source-code-folder from a build-place-folder, where Yade will be configured
and where the source code will be compiled. Here is an example for a folder structure::

	myYade/       		## base directory
		trunk/		## folder for git-handled source code, see "Download" section below
		build/		## folder in which the sources will be compiled; build-directory; use cmake here, see "Compilation.." sections below
		install/	## install folder; contains the executables


Download
^^^^^^^^^^

Installing from source, you can adopt either a release
(numbered version, which is frozen) or the current development version
(updated by the developers frequently). You should download the development
version (called ``trunk``) if you want to modify the source code, as you
might encounter problems that will be fixed by the developers. Release
versions will not be updated (except for updates due to critical and
easy-to-fix bugs), but generally they are more stable than the trunk.

#. Releases can be downloaded from the `download page <https://gitlab.com/yade-dev/trunk/-/releases>`_, as compressed archive. Uncompressing the archive gives you a directory with the sources.
#. The development version (``trunk``) can be obtained from the `code repository <https://gitlab.com/yade-dev/>`_ at GitLab.

We use `GIT <http://git-scm.com/>`_ (the ``git`` command) for code
management (install the ``git`` package on your system and create a `GitLab account <https://gitlab.com/users/sign_in>`__). From the top of of the above folder structure::

		git clone git@gitlab.com:yade-dev/trunk.git

will download the whole code repository into ``trunk`` folder. Check out :ref:`yade-gitrepo-label`
for more details on how to collaborate using ``git``.

Alternatively, a read-only checkout is possible via https without a GitLab account (easier if you don't want to modify the trunk version)::

		git clone https://gitlab.com/yade-dev/trunk.git

For those behind a firewall, you can download the sources from our `GitLab <https://gitlab.com/yade-dev>`__ repository as compressed archive.

Release and trunk sources are compiled in exactly the same way.

.. _prerequisites:

Prerequisites
^^^^^^^^^^^^^

Yade compilation and execution rely on a number of mandatory and optional external softwares; they are checked before the compilation starts.
Following dependencies are for instance mandatory:

* `cmake <http://www.cmake.org/>`_ build system
* `gcc <https://gcc.gnu.org/>`_ compiler (``g++``); other compilers will not work; you need g++>=4.2 for openMP support
* `boost <http://www.boost.org/>`_ 1.47 or later
* `Qt <http://www.qt.io/>`_ library
* `freeglut3 <http://freeglut.sourceforge.net>`_
* `libQGLViewer <http://www.libqglviewer.com>`_
* `python <http://www.python.org>`_, `numpy <https://www.numpy.org/>`_, `ipython <https://ipython.org/>`_, `sphinx <https://www.sphinx-doc.org/en/master/>`_
* `matplotlib <http://matplotlib.sf.net>`_
* `eigen <http://eigen.tuxfamily.org>`_ algebra library (minimal required version 3.2.1)
* `gdb <http://www.gnu.org/software/gdb>`_ debugger
* `sqlite3 <http://www.sqlite.org>`_ database engine

They can be installed from the command line of your Linux distribution, assuming you have root privileges.

**For Ubuntu 20.04, 18.04**, **Debian 9, 10, 11** and their derivatives, just copy&paste to the terminal the following code block for installing all mandatory and optional dependencies::

		sudo apt install cmake git freeglut3-dev libboost-all-dev fakeroot \
		dpkg-dev build-essential g++ python3-dev python3-ipython python3-matplotlib \
		libsqlite3-dev python3-numpy python3-tk gnuplot libgts-dev python3-pygraphviz \
		libvtk6-dev libeigen3-dev python3-xlib python3-pyqt5 pyqt5-dev-tools python3-mpi4py \
		python3-pyqt5.qtwebkit gtk2-engines-pixbuf python3-pyqt5.qtsvg libqglviewer-dev-qt5 \
		python3-pil libjs-jquery python3-sphinx python3-git libxmu-dev libxi-dev libcgal-dev \
		help2man libbz2-dev zlib1g-dev libopenblas-dev libsuitesparse-dev \
		libmetis-dev python3-bibtexparser python3-future coinor-clp coinor-libclp-dev \
		python3-mpmath libmpfr-dev libmpfrc++-dev libmpc-dev texlive-xetex python3-pickleshare python3-ipython-genutils

Note: on Ubuntu 22.04 and newer, the VTK library should be ``libvtk9-dev`` instead of ``libvtk6-dev``.

Most of the list above is very likely already packaged for your distribution. In case you are still confronted
with some errors concerning not available packages (e.g., package ``libmetis-dev`` is not available) it may be necessary
to add yade external ppa from https://launchpad.net/~yade-users/+archive/external (see below) as well as http://www.yade-dem.org/packages (see the top of this page)::

	sudo add-apt-repository ppa:yade-users/external
	sudo apt-get update

If you are using other distributions than Debian or its derivatives you should
install by yourself the software packages listed above. Their names in other distributions can differ from the
names of the Debian-packages.

Some of the above packages are only required for some choice of Yade compilation options, for desired Yade features, in the subsequent ``cmake`` configuration of compilation. If a required package is eventually not installed the related features will be disabled automatically with a message appearing during ``cmake`` output (at the end, in particular). Generally speaking, it is advised to watch for notes and warnings/errors, which are shown by ``cmake`` in the following.

Compilation configuration
^^^^^^^^^^^^^^^^^^^^^^^^^

Then, inside the build-directory of the above folder structure, you should call ``cmake`` to configure the compilation process, passing a path to install folder (as an option) and the path to sources::

	cmake -DCMAKE_INSTALL_PREFIX=../install ../trunk

In the above, note the ``cmake -DOPTION1=VALUE1 -DOPTION2=VALUE2`` syntax which is here applied to the lone ``CMAKE_INSTALL_PREFIX`` option, being part of a first group of ``cmake`` options that control the compilation process in itself or just slightly modify the behavior of the executable:

	* CMAKE_INSTALL_PREFIX: path where Yade should be installed (/usr/local by default)
	* CMAKE_VERBOSE_MAKEFILE: output additional information during compiling (OFF by default)
	* CHOLMOD_GPU link Yade to custom SuiteSparse installation and activate GPU accelerated PFV, see :ref:`GPUacceleration` (OFF by default)
	* DEBUG: compile in debug-mode, enabling a more convenient debugging or profiling by the user and leading to a much (1 or 2 orders of magnitude) slower executable (OFF by default)
	* DISABLE_ALL: for switching off all available boolean options, before possibly enabling explicitely just some of them, e.g. ``cmake -DDISABLE_ALL=ON -DENABLE_VTK=ON`` (OFF by default)
	* DISABLE_PKGS: comma-separated list of disabled packages i.e. names of source subdirectories under `pkg`, `preprocessing` or `postprocessing`, e.g. ``cmake -DDISABLE_PKGS=fem,pfv,image``. If empty all packages will be built. The packages `common` and `dem` are required to run, but the project can be compiled without them. (EMPTY by default)
	* ENABLE_ASAN: AddressSanitizer allows detection of memory errors, memory leaks, heap corruption errors and out-of-bounds accesses but it is slow (OFF by default)
	* ENABLE_FAST_NATIVE: use maximum optimization compiler flags including ``-Ofast`` and ``-mtune=native``. Note: ``native`` means that code will **only** run on the same processor type on which it was compiled. Observed speedup was 2% (below standard deviation measurement error) and above 5% if clang compiler was used. (OFF by default)
	* ENABLE_OAR: generate a script for oar-based task scheduler, as discussed :ref:`here<OARsection>` (OFF by default)
	* ENABLE_USEFUL_ERRORS: enable useful compiler errors which help a lot in error-free development (ON by default)
	* LIBRARY_OUTPUT_PATH: path to install libraries (lib by default)
	* MAX_LOG_LEVEL: :ref:`set maximum level <maximum-log-level>` for LOG_* macros compiled with below ``ENABLE_LOGGER``, (default is 5)
	* NOSUFFIX: do not add a suffix after binary-name, see also ``SUFFIX`` option (OFF by default)
	* PYTHON_VERSION: force Python version to the given one, e.g. ``-DPYTHON_VERSION=3.5``. Set to -1 to automatically use the last version on the system (-1 by default)
	* REAL_PRECISION_BITS, REAL_DECIMAL_PLACES: specify either of them to use a custom calculation precision of ``Real`` type. By default ``double`` (64 bits, 15 decimal places) precision is used as the ``Real`` type. See :ref:`high precision documentation<highPrecisionReal>` for additional details.
	* runtimePREFIX: used for packaging, when install directory is not the same as runtime directory (/usr/local by default)
	* SUFFIX: suffix, added after binary-names, see also ``NOSUFFIX`` option (version number by default)
	* SUITESPARSEPATH: define this variable with the path to a custom suitesparse install
	* USE_QT5: use QT5 for GUI. It is actually the only choice when GUI is requested through ``ENABLE_GUI`` option below, since libQGLViewer of version 2.6.3 and higher are compiled against Qt5 on Debian/Ubuntu operating systems (ON by default)
	* VECTORIZE: enables vectorization and alignment in Eigen3 library, experimental (OFF by default)
	* YADE_VERSION: explicitly set version number (is defined from git-directory by default)

As a more precise alternative to the above ``DISABLE_*`` options, other ``cmake`` options will select or unselect specific Yade classes for compilation, enabling or disabling additional Yade features while possibly requiring additional dependencies in form of external packages. They obey a ``ENABLE_OPTION=ON`` or ``OFF`` syntax as follows (see also the `source code <https://gitlab.com/yade-dev/trunk/blob/master/CMakeLists.txt>`_ for a most up-to-date list):

	* ENABLE_CGAL: enables a number of code sections using the `CGAL <http://www.cgal.org/>`_ library, requires ``libcgal-dev`` package (ON by default)
	* ENABLE_COMPLEX_MP: use `boost multiprecision complex <https://www.boost.org/doc/libs/1_77_0/libs/multiprecision/doc/html/boost_multiprecision/tut/complex.html>`__ and `mpc <http://www.multiprecision.org/mpc/>`_ (as an extension to MPFR, see ``ENABLE_MPFR``) for ``ComplexHP<N>``, otherwise use ``std::complex<RealHP<N>>``. See :ref:`high precision documentation<highPrecisionReal>` for additional details. Requires ``libmpc-dev`` (ON by default if possible: requires boost >= 1.71)
	* ENABLE_DEFORM: enable the constant volume deformation approach for bodies [Haustein2017]_ (OFF by default)
	* ENABLE_FEMLIKE: enable FEM-like meshed solids (ON by default)
	* ENABLE_GL2PS: enable GL2PS-option (ON by default)
	* ENABLE_GTS: enable GTS-option (ON by default)
	* ENABLE_GUI: enable a Qt5 GUI. Requires ``python-pyqt5 pyqt5-dev-tools`` (ON by default)
	* ENABLE_LBMFLOW: enable LBM computations, e.g. the use of :yref:`HydrodynamicsLawLBM` (ON by default)
	* ENABLE_LS_DEM: enable a :yref:`LevelSet` shape description (ON by default)
	* ENABLE_LINSOLV: enable the use of optimized algebra libraries `SuiteSparse <http://www.suitesparse.com>`_ (sparse algebra, requires eigen>=3.1), `OpenBLAS <http://www.openblas.net/>`_ (optimized and parallelized alternative to the standard blas+lapack) and `Metis <http://glaros.dtc.umn.edu/gkhome/metis/metis/overview/>`_ (matrix preconditioning) for the optional fluid coupling :yref:`FlowEngine`, see ``ENABLE_PFVFLOW`` below. Requires ``libopenblas-dev libsuitesparse-dev libmetis-dev`` packages (ON by default)
	* ENABLE_LIQMIGRATION: enable LIQMIGRATION-option, see [Mani2013]_ for details (OFF by default)
	* ENABLE_LOGGER: provides :ref:`logging<logging>` possibilities for each class thanks to `boost::log <https://www.boost.org/doc/libs/release/libs/log/>`_ library. See also ``MAX_LOG_LEVEL`` in the above (ON by default)
	* ENABLE_MASK_ARBITRARY: enable arbitrary precision of bitmask variables (only ``Body::groupMask`` yet implemented) (experimental). If ON, use -DMASK_ARBITRARY_SIZE=int to set number of used bits (256 by default) (OFF by default)
	* ENABLE_MPFR: use `mpfr <https://www.mpfr.org/>`_ in ``C++`` and `mpmath <http://mpmath.org/>`_ in ``python`` for higher precision ``Real`` or for CGAL exact predicates, see :ref:`high precision documentation<highPrecisionReal>` for more details. Requires ``python3-mpmath libmpfr-dev libmpfrc++-dev`` packages (OFF by default)
	* ENABLE_MPI: enable MPI environment and communication thanks to `OpenMPI <https://www.open-mpi.org/software/>`_ and `python3-mpi4py <https://bitbucket.org/mpi4py/>`_ (see also `there <https://mpi4py.readthedocs.io/en/stable/>`_), for parallel distributed computing (distributed memory) and Yade-OpenFOAM coupling. Requires ``python3-mpi4py`` (ON by default)
	* ENABLE_MULTI_REAL_HP: allow using twice, quadruple or higher precisions of ``Real`` as ``RealHP<2>``, ``RealHP<4>`` or ``RealHP<N>`` in computationally demanding sections of ``C++`` code. See :ref:`high precision documentation<highPrecisionReal>` for additional details (ON by default).
	* ENABLE_OPENMP: enable OpenMP-parallelizing of Yade execution (ON by default)
	* ENABLE_PARTIALSAT : enable the partially saturated clay engine :yref:`PartialSatClayEngine`, under construction (ON by default)
	* ENABLE_PFVFLOW: enable PFV :yref:`FlowEngine` (ON by default)
	* ENABLE_POTENTIAL_BLOCKS: enable :yref:`PotentialBlock` shape description thanks for instance to the `COIN-OR <https://github.com/coin-or/Clp>`_ Linear Programming Solver, requires ``coinor-clp coinor-libclp-dev libopenblas-dev`` (ON by default)
	* ENABLE_POTENTIAL_PARTICLES: enable :yref:`PotentialParticle` shape description, requires ``libopenblas-dev`` (ON by default)
	* ENABLE_PROFILING: enable profiling, e.g., shows some more metrics, which can define bottlenecks of the code (OFF by default)
	* ENABLE_SPH: enable Smoothed Particle Hydrodynamics (OFF by default)
	* ENABLE_THERMAL : enable :yref:`ThermalEngine` (ON by default, experimental)
	* ENABLE_TWOPHASEFLOW: enable :yref:`TwoPhaseFlowEngine` (ON by default)
	* ENABLE_VTK: enable exports of data using the `VTK <http://www.vtk.org/>`_ library, e.g. :yref:`VTKRecorder` engine, requires a ``libvtk*-dev`` (e.g., ``libvtk9-dev`` on Ubuntu 22.04) package (ON by default)

Maintaining a consistent choice for options values, in addition to using the same version of source code, is often necessary for successfully reloading previous Yade saves, see :yref:`O.load<Omega.load>`.

For using more extended parameters of cmake, please follow the corresponding
documentation on `https://cmake.org/documentation <https://cmake.org/documentation/>`_.

.. warning:: If you have Ubuntu 14.04 Trusty, you need to add -DCMAKE_CXX_FLAGS=-frounding-math
 during the configuration step of compilation (see below) or to install libcgal-dev
 from our `external PPA <https://launchpad.net/~yade-users/+archive/external/>`_.
 Otherwise the following error occurs on AMD64 architectures::

    terminate called after throwing an instance of 'CGAL::Assertion_exception'
    what():  CGAL ERROR: assertion violation!
    Expr: -CGAL_IA_MUL(-1.1, 10.1) != CGAL_IA_MUL(1.1, 10.1)
    File: /usr/include/CGAL/Interval_nt.h
    Line: 209
    Explanation: Wrong rounding: did you forget the  -frounding-math  option if you use GCC (or  -fp-model strict  for Intel)?
    Aborted


.. _yadeCompilation:

Compilation and usage
^^^^^^^^^^^

If cmake finishes without errors, you will see all enabled
and disabled options at the end. Then start the actual compilation process with::

	make

The compilation process can take a considerable amount of time, be patient.
If you are using a multi-core system you can use the parameter ``-j`` to speed-up the compilation
and split the compilation onto many cores. For example, on 4-core machines
it would be reasonable to set the parameter ``-j4``. Note, Yade requires
approximately 3GB RAM per core for compilation, otherwise the swap-file will be used
and compilation time dramatically increases.

The installation is performed with the following command::

	make install

The ``install`` command will in fact also recompile if source files have been modified.
Hence there is no absolute need to type the two commands separately. You may receive make errors if you don't have permission to write into the target folder.
These errors are not critical but without writing permissions Yade won't be installed in /usr/local/bin/.

After the compilation finished successfully,
the new built can be started by navigating to /path/to/installfolder/bin and calling yade via (based on version yade-2014-02-20.git-a7048f4)::

    cd /path/to/installfolder/bin
    ./yade-2014-02-20.git-a7048f4

.. comment: is it possible to invoke python yade.config.revision and put it above as a text in the doc?

For building the documentation you should at first execute the command ``make install``
and then ``make doc`` to build it, provided that package ``texlive-xetex`` is present. On some multi-language systems an error ``Building format(s) --all. This may take some time... fmtutil failed.`` may occur, in that case a package ``locales-all`` is required.

The generated files will be stored in your current
install directory /path/to/installfolder/share/doc/yade-your-version. Once again writing permissions are necessary for installing into /usr/local/share/doc/. To open your local documentation go into the folder html and open the file index.html with a browser.

``make manpage`` command generates and moves manpages in a standard place.
``make check`` command executes standard test to check the functionality of the compiled program.

Yade can be compiled not only by GCC-compiler, but also by `CLANG <http://clang.llvm.org/>`_
front-end for the LLVM compiler. For that you set the environment variables CC and CXX
upon detecting the C and C++ compiler to use::

	export CC=/usr/bin/clang
	export CXX=/usr/bin/clang++
	cmake -DOPTION1=VALUE1 -DOPTION2=VALUE2

Clang does not support OpenMP-parallelizing for the moment, that is why the
feature will be disabled.

Supported linux releases
^^^^^^^^^^^^^^^^^^^^^^^^

`Currently supported <https://gitlab.com/yade-dev/trunk/pipelines?scope=branches>`_ [#buildLog]_ linux releases and their respective `docker <https://docs.docker.com/>`_ `files <https://docs.docker.com/engine/reference/builder/>`_ are:

* `Ubuntu 18.04 bionic <https://gitlab.com/yade-dev/docker-yade/blob/ubuntu18.04/Dockerfile>`_
* `openSUSE 15 <https://gitlab.com/yade-dev/docker-yade/blob/suse15/Dockerfile>`_

These are the bash commands used to prepare the linux distribution and environment for installing and testing yade.
These instructions are automatically performed using the `gitlab continuous integration <https://docs.gitlab.com/ee/ci/quick_start/>`_ service after
each merge to master. This makes sure that yade always works correctly on these linux distributions.
In fact yade can be installed manually by following step by step these instructions in following order:

1. Bash commands in the respective Dockerfile to install necessary packages,

2. do ``git clone https://gitlab.com/yade-dev/trunk.git``,

3. then the ``cmake_*`` commands in the `.gitlab-ci.yml file <https://gitlab.com/yade-dev/trunk/blob/master/.gitlab-ci.yml>`_ for respective distribution,

4. then the ``make_*`` commands to compile yade,

5. and finally the ``--check`` and ``--test`` commands.

6. Optionally documentation can be built with ``make doc`` command, however currently it is not guaranteed to work on all linux distributions due to frequent interface changes in `sphinx <http://www.sphinx-doc.org/en/master/>`_.

These instructions use ``ccache`` and ``ld.gold`` to :ref:`speed-up compilation <speed-up>` as described below.

.. [#buildLog] To see details of the latest build log click on the *master* branch.

Python 2 backward compatibility
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Following the end of Python 2 support (beginning of 2020), Yade compilation on a Python 2 ecosystem is no longer garanteed since the `6e097e95 <https://gitlab.com/yade-dev/trunk/-/tree/6e097e95368a9c63ce169a040f418d30c7ba307c>`_ trunk version. Python 2-compilation of the latter is still possible using the above ``PYTHON_VERSION`` cmake option, requiring Python 2 version of prerequisites packages whose list can be found in the corresponding paragraph (Python 2 backward compatibility) of the `historical doc <https://gitlab.com/yade-dev/trunk/-/blob/6e097e95368a9c63ce169a040f418d30c7ba307c/doc/sphinx/installation.rst>`_.

Ongoing development of Yade now assumes a Python 3 environment, and you may refer to some notes about :ref:`converting Python 2 scripts into Python 3<convert-python2-to3>` if needed.

.. _speed-up:

Speed-up compilation
---------------------

Compile with ccache
^^^^^^^^^^^^^^^^^^^^^

Caching previous compilations with `ccache <https://ccache.samba.org/>`_ can significantly speed up re-compilation::

	cmake -DCMAKE_CXX_COMPILER_LAUNCHER=ccache [options as usual]

Additionally one can check current ccache status with command ``ccache --show-stats`` (``ccache -s`` for short) or change the default `cache size <https://wiki.archlinux.org/index.php/ccache#Set_maximum_cache_size>`_ stored in file ``~/.ccache/ccache.conf``.

Compile with distcc
^^^^^^^^^^^^^^^^^^^^^

When spliting the compilation on many cores (``make -jN``), ``N`` is limited by the available cores and memory. It is possible to use more cores if remote computers are available, distributing the compilation with `distcc <https://wiki.archlinux.org/index.php/Distcc>`_  (see distcc documentation for configuring slaves and master)::

	export CC="distcc gcc"
	export CXX="distcc g++"
	cmake [options as usual]
	make -jN

The two tools can be combined, adding to the above exports::

	export CCACHE_PREFIX="distcc"

Compile with cmake UNITY_BUILD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This option concatenates source files in batches containing several ``*.cpp`` each, in order to share the overhead of include directives (since most source files include the same boost headers, typically). It accelerates full compilation from scratch (quite significantly). It is activated by adding the following to cmake command, ``CMAKE_UNITY_BUILD_BATCH_SIZE`` defines the maximum number of files to be concatenated together (the higher the better, main limitation might be available RAM)::

	-DCMAKE_UNITY_BUILD=ON -DCMAKE_UNITY_BUILD_BATCH_SIZE=18

This method is helpless for incremental re-compilation and might even be detrimental since a full batch has to be recompiled each time a single file is modified. If it is anticipated that specific files will need incremental compilation they can be excluded from the unity build by assigning their full path to cmake flag ``NO_UNITY`` (a single file or a comma-separated list)::

	-DCMAKE_UNITY_BUILD=ON -DCMAKE_UNITY_BUILD_BATCH_SIZE=18 -DNO_UNITY=../trunk/pkg/dem/CohesiveFrictionalContactLaw.cpp

Link time
^^^^^^^^^^^^^^^^^^^^^

The link time can be reduced by changing the default linker from ``ld`` to ``ld.gold``. They are both in the same package ``binutils`` (on opensuse15 it is package ``binutils-gold``). To perform the switch execute these commands as root::

	ld --version
	update-alternatives --install "/usr/bin/ld" "ld" "/usr/bin/ld.gold" 20
	update-alternatives --install "/usr/bin/ld" "ld" "/usr/bin/ld.bfd" 10
	ld --version

To switch back run the commands above with reversed priorities ``10`` ↔ ``20``. Alternatively a manual selection can be performed by command: ``update-alternatives --config ld``.

Note: ``ld.gold`` is incompatible with the compiler wrapper ``mpicxx`` in some distributions, which is manifested as an error in the ``cmake`` stage.
We do not use ``mpicxx`` for our gitlab builds currently. If you want to use it then disable ``ld.gold``. Cmake ``MPI``-related failures have also been reported without the ``mpicxx`` compiler, if it happens then the only solution is to disable either ``ld.gold`` or the ``MPI`` feature.

Cloud Computing
----------------

It is possible to exploit cloud computing services to run Yade. The combo Yade/Amazon Web Service has been found to work well, namely. Detailed instructions for migrating to amazon can be found in the section :ref:`CloudComputing`.

GPU Acceleration
----------------

The FlowEngine can be accelerated with CHOLMOD's GPU accelerated solver. The specific hardware and software requirements are outlined in the section :ref:`GPUacceleration`.

Special builds
--------------

The software can be compiled by a special way to find some specific bugs and problems in it: memory corruptions, data races, undefined behaviour etc.


The listed sanitizers are runtime-detectors. They can only find the problems in the code, if the particular part of the code
is executed. If you have written a new C++ class (constitutive law, engine etc.) try to run your Python script with
the sanitized software to check, whether the problem in your code exist.

.. _address-sanitizer:

AddressSanitizer
^^^^^^^^^^^^^^^^

`AddressSanitizer <https://clang.llvm.org/docs/AddressSanitizer.html>`_ is a memory error detector, which helps to find heap corruptions,
out-of-bounds errors and many other memory errors, leading to crashes and even wrong results.

To compile Yade with this type of sanitizer, use ENABLE_ASAN option::

	cmake -DENABLE_ASAN=1

The compilation time, memory consumption during build and the size of build-files are much higher than during the normall build.
Monitor RAM and disk usage during compilation to prevent out-of-RAM problems.

To find the proper libasan library in your particular distribution, use ``locate`` or ``find /usr -iname "libasan*so"`` command. Then, launch your yade executable in connection with that libasan library, e.g.::

    LD_PRELOAD=/some/path/to/libasan.so yade

By default the leak detector is enabled in the asan build. Yade is producing a lot of leak warnings at the moment.
To mute those warnings and concentrate on other memory errors, one can use detect_leaks=0 option. Accounting for the latter, the full command
to run tests with the AddressSanitized-Yade on Debian 10 Buster is::

	ASAN_OPTIONS=detect_leaks=0:verify_asan_link_order=false yade --test

If you add a new check script, it is being run automatically through the AddressSanitizer in the CI-pipeline.

Yubuntu
-------

If you are not running a Linux system there is a way to create an Ubuntu `live-usb <http://en.wikipedia.org/wiki/Live_USB>`_ on any usb mass-storage device (minimum size 10GB). It is a way to boot the computer on a linux system with Yadedaily pre-installed without affecting the original system. More informations about this alternative are available `here <http://people.3sr-grenoble.fr/users/bchareyre/pubs/yubuntu/>`_ (see the README file first). Note that the images there date back from 2018 and use ubuntu16.04, for newer versions of yade see below.

Alternatively, images of a linux virtual machine can be downloaded `here (ubuntu20.04) <https://yade-dem.org/publi/virtual/>`_ , or for older (ubuntu16.04) versions `here <http://people.3sr-grenoble.fr/users/bchareyre/pubs/yubuntu/>`_. They should run on any system with a virtualization software (tested with VirtualBox and VMWare).


