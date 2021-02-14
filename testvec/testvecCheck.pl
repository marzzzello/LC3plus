#!/usr/bin/perl -w

# ======================================
# LC3plus ETSI Testvectors script V0.0.1
# ======================================
#
# /******************************************************************************
# *                        ETSI TS 103 634 V1.1.1                               *
# *              Low Complexity Communication Codec Plus (LC3plus)              *
# *                                                                             *
# * Copyright licence is solely granted through ETSI Intellectual Property      *
# * Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
# * estoppel or otherwise.                                                      *
# ******************************************************************************/

##################################################
my $EXE_TST = './../src/fixed_point/LC3plus'; # Path to test LC3plus executable
my $md5_bin = './md5_bin.txt';
my $md5_dec = './md5_dec.txt';
my $MY_MD5  = 'md5sum'; # System dependent MD5 call
my $epf = '-epf ./plc_fer.dat';
##################################################

use strict;
use File::Basename;
use File::Path 'rmtree';
use POSIX;

my $VERSION = 'V0.0.1';
my $timestamp = strftime "%m_%d_%Y_%H_%M_%S", localtime;
my $tmp_folder = "lc3plus_testvectors".$timestamp;
my $inputFile = "./input/thetest";
my $output_folder_stream_tst = $tmp_folder."/bitstream_tst";
my $output_folder_decoded_tst = $tmp_folder."/decoded_tst";
my $report = "lc3plus_testvectors_report_".$timestamp.".txt";
my $fh;
my $quiet = '>/dev/null 2>&1';

my $testvectors_fail = 0;

my @SR = (8000, 16000, 24000, 32000, 44100, 48000);
my @BR_8 = (32000);
my @BR_16 = (32000);
my @BR_24 = (48000);
my @BR_32 = (48000);
my @BR_441 = (64000);
my @BR_48 = (64000);

my @rates = ([@BR_8], [@BR_16], [@BR_24], [@BR_32], [@BR_441], [@BR_48]);
my $srLen = @SR;
my @EP_MODES = (0, 4);

# Set MD5 command according to OS
getOS();

# Get command line arguments
my ($create, $test, $clean, $log) = (0, 0, 0, 0);
getArgs(\$create, \$clean, \$log, \$quiet, \@ARGV);

$test = 1; # Run complete test by default
if ($create)
{
    $test = 0; # Don't run test if only create option selected
}

# Check if both executables are here
checkExe($EXE_TST, $md5_bin, $md5_dec);

# Check if MD5 command exists
checkMD5($MY_MD5, $quiet);

# Create all directories
createDirs($tmp_folder, $output_folder_stream_tst, $output_folder_decoded_tst);

# Check if input files are present
checkDir($inputFile);

# Create report file if user wants logging
if ($log)
{
    open($fh, '>', $report) or die "Could not open file '$report' $!";
}

open(my $md5stream, '<', $md5_bin) or die "Could not open file '$md5_bin' $!";
open(my $md5decoded, '<', $md5_dec) or die "Could not open file '$md5_dec' $!";

print("Script started... \n");
if ($create || $test)
{
    print("Creating files with test executable...");
    
    for (my $srIndex = 0; $srIndex < $srLen; $srIndex++)
    {
        my $sr = $SR[$srIndex];
        my $input = "$inputFile".(floor($sr/1000)).".wav";
        my $base = basename($input);
        my $filename = $base;
        $base =~ s/\.[^.]+$//;

        my $ref = $rates[$srIndex];

        foreach my $br (@$ref)
        {
            foreach my $ep (@EP_MODES)
            {
                my $output_stream = $output_folder_stream_tst."/".$base."_".$br."_EP".$ep.".lc3plus";
                system("$EXE_TST $epf -epmode $ep -E $input $output_stream $br $quiet");

                if ($log)
                {
                    print $fh "$EXE_TST $epf -epmode $ep -E $input $output_stream $br $quiet\n";
                }

                my $output_decoded = $output_folder_decoded_tst."/".$base."_".$br."_EP".$ep.".wav";
                system("$EXE_TST $epf -epmode $ep $input $output_decoded $br $quiet");

                if ($log)
                {
                    print $fh "$EXE_TST $epf -epmode $ep $input $output_decoded $br $quiet\n";
                }

                if ($test)
                {
                    compare($output_stream, $output_decoded, $log, $fh, $md5stream, $md5decoded, \$testvectors_fail);
                }
            }
        }
    }
    
    print("...done!\n");
}

if ($clean)
{
    cleanup($tmp_folder)
}

if ($test)
{
    my $result = "passed";
    if ($testvectors_fail)
    {
        $result = "NOT passed";
    }

    print("\nThe testvector check test was $result!\n");

    # Close text file for logging
    if ($log)
    {
        print $fh "\n The testvectors check test was $result!\n";
        close $fh;
    }
}

if ($log)
{
    print("Please see logfile for more information: $report\n");
}



# Functions
sub cleanup
{
   my ($tmp_dir) = @_;
   rmtree($tmp_dir);
}

sub checkMD5
{
    my ($cmd, $quiet) = @_;
    
    my $ret = system("$cmd $0 $quiet");
    
    if ($ret != 0)
    {
        print("Error: cannot find md5 command: $cmd\n");
        print("Please adjust md5 command in script!\n");
        exit(0);
    }
}

sub compare
{
    my ($stream_tst, $dec_tst, $log, $filehandle_log, $md5stream, $md5dec, $testvectors_fail) = @_;
    my $ret;
    
    my $md5_ref_stream = <$md5stream>;
    my $md5_ref_dec = <$md5dec>;
    my $base_stream = basename($stream_tst);
    my $base_dec = basename($dec_tst);
    chomp($md5_ref_stream); chomp($md5_ref_dec);
    
    my @stream=split(/:/,$md5_ref_stream);
    my @dec=split(/:/,$md5_ref_dec);
    
    if ($stream[0] ne $base_stream || $dec[0] ne $base_dec)
    {
        print("Error: strings do not match!\n");
        exit(1);
    }
    
    if ($log)
    {
        print $filehandle_log "Checking MD5 of files: $stream_tst and $dec_tst\n";
        print $filehandle_log "Reference MD5: stream=$stream[1], dec=$dec[1]\n";
    }
    
    my $md5 = qx($MY_MD5 $stream_tst);
    my ($hash_stream) = $md5 =~ /([a-z0-9]{32})/i;
    $hash_stream =~ s/^.//;
    
    $md5 = qx($MY_MD5 $dec_tst);
    my ($hash_dec) = $md5 =~ /([a-z0-9]{32})/i;
    $hash_dec =~ s/^.//;
    
    if ($stream[1] ne $hash_stream)
    {
        ${$testvectors_fail} = 1;
    }
    
    if ($dec[1] ne $hash_dec)
    {
        ${$testvectors_fail} = 1;
    }
    
    if ($log)
    {
        print $filehandle_log "Test MD5: stream=$hash_stream, dec=$hash_dec\n";
        if ($stream[1] ne $hash_stream || $dec[1] ne $hash_dec)
        {
            print $filehandle_log "The MD5 sums do NOT match!\n";
        } else {
            print $filehandle_log "The MD5 sums do match!\n";
        }
    }
}

sub checkDir
{
    my ($fileBase) = @_;
    
    my $test8 = $fileBase."8.wav";
    my $test16 = $fileBase."16.wav";
    my $test24 = $fileBase."24.wav";
    my $test32 = $fileBase."32.wav";
    my $test44 = $fileBase."44.wav";
    my $test48 = $fileBase."48.wav";
    
    if (! -e $test8 || ! -e $test16 || ! -e $test24 || ! -e $test32 || ! -e $test44 || ! -e $test48)
    {
        print("Cannot find input files! Please provide them in ./input/\n");
        printUsage();
    }
}

sub createDirs
{
    my ($tmp, $tmp_tst_stream, $tmp_tst_dec) = @_;
    
    mkdir($tmp);
    mkdir($tmp_tst_stream);
    mkdir($tmp_tst_dec);
}

sub getArgs 
{
    my ($create, $clean, $log, $quiet, $args) = @_;
    
    my @arg = @{$args};
    my $size = @arg;
    my $help = 0;
    
    foreach my $argument (@arg)
    {
        if ($argument eq '-create')
        {
            ${$create} = 1;
        } elsif ($argument eq '-clean')
        {
            ${$clean} = 1;
        } elsif ($argument eq '-log')
        {
            ${$log} = 1;
        } elsif ($argument eq '-verbose')
        {
            ${$quiet} = '';
        } elsif ($argument eq '-h' || $argument eq '--h' || $argument eq '-help' || $argument eq '--help')
        {
            $help = 1;
        } else {
            print("Unknown parameter detected: $argument!\n");
            printUsage();
        }
    }
    
    if ($help)
    {
        printUsage();
    }
}

sub checkExe
{
    my ($exe_tst, $md5bin, $md5ref) = @_;
    
    if (! -e $exe_tst)
    {
        print("Cannot find test executable! Please place LC3plus in src/fixed_point/ or adjust the variable EXE_TST in this script.\n");
        printUsage();
    }
    
    if (! -e $md5bin)
    {
        print("Cannot find MD5 reference file: $md5bin!\n");
        printUsage();
    }
    
    if (! -e $md5ref)
    {
        print("Cannot find MD5 reference file: $md5ref!\n");
        printUsage();
    }
}

sub getOS
{
    my $osname = $^O;
    
    if ($osname eq "darwin")
    {
        $MY_MD5 = "md5";
    } else {
        $MY_MD5 = "md5sum";
    }
}

sub printUsage
{
    print("\nLC3plus ETSI Testvectors script $VERSION\n");
    print("Default: Create files and compare MD5 hash with reference values.\n");
    print("Usage: lc3plus_testvectors.pl [-create] [-clean] [-log] [-verbose]\n");
    print("[-create] : Only create bitstreams and decoded files for testvectors check using test executable.\n");
    print("[-clean] : Delete all temporary files after usage.\n");
    print("[-log] : Log all commands.\n");
    print("[-verbose] : Show all system output (supressed by default).\n");
    print("[-h] : Display this information.\n");
    exit(0);
}
