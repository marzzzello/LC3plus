#! /bin/bash
# set -x

# /******************************************************************************
# *                        ETSI TS 103 634 V1.1.1                               *
# *              Low Complexity Communication Codec Plus (LC3plus)              *
# *                                                                             *
# * Copyright licence is solely granted through ETSI Intellectual Property      *
# * Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
# * estoppel or otherwise.                                                      *
# ******************************************************************************/

cat << EOF
(C) Copyright Ericsson AB and Fraunhofer Gesellschaft zur Foerderung der
    angewandten Forschung e.V. for its Fraunhofer IIS 2019

Example script file how to simulate a streaming session with redundancy LC3plus frames.    

EOF

if [ "$#" -ne 5 ]; then
cat << EOF
usage: lc3plus_redundancy_example.sh input.wav primary_rate secondary_rate error_pattern offset

where
  input.wav         input wav audio file
  primary_rate      bit rate of primary LC3plus frames
  secondary_rate    bit rate of secondary LC3plus frames
  error_apttern     error pattern file, 16bit or g192 format
  offset            temporal offset of secondary frames compared to primary frames in 20ms
    
The decoder wav is available in the file out.dec.wav

EOF
exit
fi

INPUT=$1                 # input wav file
OUTPUT=out.dec.wav       # output wac file
MAINRATE=$2              # primary frame rate in kbbs
HELPRATE=$3              # secondary rate in kbps
HELPBANDWIDTH=4000       # bandwidth of helper payload
ERRORPATTERN=$4          # error pattern file
OFFSET=$5                # helper payload offset in packets
BFI=3                    # signal of helper payload

if [ $MAINRATE == $HELPRATE ]; then
    BFI=0                # if helper paylaod rate equal to main payload, signal regular payload
fi


EXE=../src/fixed_point/LC3plus

echo "Coding primary payload"
$EXE -E -formatG192 $INPUT full_bitstream_file.g192 $MAINRATE
echo "Coding secondary payload"
$EXE -E -formatG192 -bandwidth $HELPBANDWIDTH $INPUT help_bitstream_file.g192 $HELPRATE

echo "Combining payload"
./lc3plus_redundancy_simulator.pl full_bitstream_file.g192 help_bitstream_file.g192 $ERRORPATTERN output_bitsream.g192 $OFFSET 3

$EXE -D -formatG192 output_bitsream.g192 out.dec.wav 64000





