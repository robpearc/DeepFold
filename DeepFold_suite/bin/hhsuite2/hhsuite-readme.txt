This is a fork for the old hhsuite 2.x, downloaded from 
http://gwdu111.gwdg.de/~compbiol/data/hhsuite/releases/all/hhsuite-2.0.16.tar.gz
with patch from 
https://gist.github.com/milot-mirdita/fd4b193a2423cc8e71868a3a68ef940f

In addition to hhsuite 2.0.16, this package also include, under bin/ folder, 
binary executables from:
kClust 1.0
Clustal Omega 1.2.4
PSIPRED 4.0.1
legacy BLAST 2.2.26
solvpred (from MetaPSICOV 1.04)

It also contain, under bin/ folder, custom binary executables from:
qhmmer (https://github.com/kad-ecoli/qhmmer)
MSAParser (https://github.com/kad-ecoli/MSAParser)

The following is the original README file for hhsuite 2.0.15 and 2.0.16.

*****************************************************************************
*****************************************************************************
   HH-suite for sensitive sequence searching version 2.0.15 (June 2012)

 (C) Johannes Soeding, Michael Remmert, Andreas Hauser, Andreas Biegert 2012

This README only serves as a quick start guide.
For full documentation see the user guide in hhsuite-userguide.pdf

########################################################################
The hhsuite contains in file hhprefilter.C code adapted from Michael 
Farrar (http://sites.google.com/site/farrarmichael/smith-waterman). 
His code is marked in the file hhprefilter.C. For the copy right of that 
code, please see the LICENSE file that comes with HHsuite.
Reference: Farrar M. Striped Smith-Waterman speeds database searches six 
times over other SIMD implementations. Bioinformatics. 2007, 23:156-61.
Many posthumous thanks to Michael Farrar for his great code!
########################################################################


*****************************************************************************
 Installation
*****************************************************************************

1. Downloading

Download the sources from ftp://toolkit.genzentrum.lmu.de/HH-suite/, e.g.

$ mkdir ~/programs/hh/
$ cd ~/programs/hh/
$ wget ftp://toolkit.genzentrum.lmu.de/HH-suite/hhsuite-latest.tar.gz

Unzip and untar

$ tar -xzvf hhsuite-latest.tar.gz

This will unpack the sources to hhsuite-<VERSION>.

2. Compilation 

The following dependencies exist: perl, libpng, libz

$ cd hhsuite-<VERSION>/
$ make

Compilation produces by default static binaries. If you encounter missing
library errors, also make sure you have the static versions installed, e.g.
glibc-static, zlib-static and libpng-static. Package names and whether they are
automatically installed or varies with distribution.
On RHEL/Fedora/SL/CentOS etc. eg.:
$ yum install glibc-static libpng-static zlib-static

A dynamically linked version of the programs can be compiled with:

$ make all

On Mac OS X only dynamic linking is supported.

Notes on SSE support.

To compile code for CPUs without SSE3 you can compile with NO_SSE3 set:
$ make NO_SSE3=1
To see whether your CPU supports SSE3 you can check the output of:
$ grep pni /prc/cpuinfo

If your CPU supports SSE4.1 you can enable with with WITH_SSE41=1:
$ make WITH_SSE41=1
To see whether your CPU supports SSE4.1 you can check the output of:
$ grep sse4_1 /prc/cpuinfo

If you only use the programs on your local machine or of the exact type
you can further enable the use of all other instructions supported by
your CPU with WITH_NATIVE:
$ make WITH_NATIVE=1
This is independent of the SSE flags.


3. Installation

Either install in current directory:

$ make install

Or set INSTALL_DIR to the base directory (absolute path), e.g. /usr/local,
where you want to install:

$ make install INSTALL_DIR=/usr/local


4. Set paths 

4.1. In your shell set environment variable HHLIB to $INSTALL_DIR/lib/hh, 
e.g (for bash, zsh, ksh):

$ export HHLIB=/usr/local/lib/hh

HHsearch and HHblits look for the column state library file cs219.lib
and the context library file context_data.lib in $HHLIB/data/. The hh-suite
perl scripts also read HHLIB to locate 

4.2. Specify paths to BLAST, PSIPRED, PDB, DSSP etc. in 
$INSTALL_DIR/lib/hh/scripts/HHPaths.pm
where they are read by the perl scripst in hh-suite.

4.3 For dynamic builds only

To use a dynamic build, in the lib subdirectory must either be in the system
library path, e.g. the ones configured in /etc/ld.so.conf on linux, or a variable
must be set for the run-time linker to find the libraries.

Assuming INSTALL_DIR was set to /opt/hhsuite, on linux LD_LIBRARY_PATH must
be set, e.g. in bourne shell format (bash, zsh etc):

$ export LD_LIBRARY_PATH=/opt/hhsuite/lib

on Mac OSX DYLD_LIBRARY_PATH must be set, e.g.:

$ export DYLD_LIBRARY_PATH=/opt/hhsuite/lib


5. Download Databases

Download current databases from our FTP-server:
ftp://toolkit.genzentrum.lmu.de/HH-suite/databases/hhsuite_dbs/
ftp://toolkit.genzentrum.lmu.de/HH-suite/databases/hhsearch_dbs/

The latter directory contains old-style HHsearch databases which will in the longer
run be replaced by the HH-suite formatted versions.

To build up multiple sequences alignments using HHblits, either one of uniprot20 
or nr20 dbs is sufficient.


*****************************************************************************
* Usage
*****************************************************************************


For performing a single search iteration of HHblits, run HHblits with the 
following command:

$ hhblits -i <input-file> -o <result-file> -n 1 -d <database-basename>

For generating an alignment of homologous sequences:

$ hhblits -i <input-file> -o <result-file> -oa3m <result-alignment> -d <database-basename>

You can get a detailed list of options for HHblits by running HHblits with the "-h" option.

Example (if HHLIB is set):

$ hhblits -d $HHLIB/hhblits_dbs/uniprot20_02Sep11 -i $HHLIB/data/query.a3m


*****************************************************************************
* License
*****************************************************************************

The HHsearch/HHblits software package is distributed under Gnu Public Licence, Version 3.
This means that the HH-suite is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License for more details.

See the copy of the GNU General Public License in the LICENSE file. 
If you do not have this file, see http://www.gnu.org/licenses/


-- 
For full documentation see the user guide in hhsuite-userguide.pdf


We are very grateful for bug reports! 
Please contact us at soeding@genzentrum.lmu.de
