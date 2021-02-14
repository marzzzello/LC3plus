/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

=====================================
LC3plus ETSI Testvector script V0.0.1
=====================================


=====================
WHAT IS THIS SCRIPT?
=====================

This is the LC3plus ETSI testvector script. It checks whether the output of
your LC3plus build produces the expected results and matches the precalculated
MD5 hashes for a number of operating points. It is meant to be used for the
fixed-point version of LC3plus only.

The following configurations are tested:

 Samplingrate [Hz] | Bitrate [bps]  | EP Mode [0 = off, 4 = highest protection]
-------------------+----------------+------------------------------------------
 8000              | 32000          | 0, 4
 16000             | 32000          | 0, 4 
 24000             | 48000          | 0, 4
 32000             | 48000          | 0, 4
 44100             | 64000          | 0, 4
 48000             | 64000          | 0, 4
 
For all test cases, PLC is also triggered using the error pattern file 'plc_fer.dat'.
The script compares the MD5 hashes of bitstreams and output WAVEs.

IMPORTANT NOTICE: This is not the conformance tool and tests only a limited number
of operating points of the LC3plus codec as a pre-test. The full conformance tool
can be found in 'PACKAGE/conformance/lc3_conformance.py'.

=======================
TO DO BEFORE FIRST USE
=======================

1) Build the LC3plus executable in the 'src/fixed_point' folder. The script expects an
executable called 'LC3plus' in the '../src/fixed_point' folder. Alternatively adjust
the following line in the script:

<<<my $EXE_TST = './../src/fixed_point/LC3plus';>>>

to match the path to your LC3plus test executable.

======
USAGE
======

perl testvectorCheck.pl [-create] [-clean] [-log] [-verbose]

Command line options:

-create  : Only create bitstreams and decoded files with test executable
-clean   : Cleanup after usage -> delete all temporary files after script has finished
-log     : Log all commands into logfile
-verbose : Do not hide system call output
-h       : Print help screen

Default usage: Create vectors and test MD5 hash against reference values.


=======
OUTPUT
=======

The pass/fail result is printed on the console. If [-log] was also selected, the pass/fail 
result will also be printed into the logfile. Unless you use [-clean], the script will save 
all files in a temporary folder labeled lc3plus_testvec_MM_DD_YYYY_HH_MM_SS.

