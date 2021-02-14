#!/usr/bin/perl -w

# /******************************************************************************
# *                        ETSI TS 103 634 V1.1.1                               *
# *              Low Complexity Communication Codec Plus (LC3plus)              *
# *                                                                             *
# * Copyright licence is solely granted through ETSI Intellectual Property      *
# * Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
# * estoppel or otherwise.                                                      *
# ******************************************************************************/

use strict;
use File::Copy;

our %config;

our $ep_offset = 0;    # offset to ep file
our $framing   = 20;   # 10 or 20 ms frame size

if ($#ARGV < 5) {
    printHelp();
}

&parseCommandline();

# copy main config to output config
$config{mainCfg} = $config{mainBs}.".cfg";
$config{outCfg} = $config{output}.".cfg";
copy($config{mainCfg},$config{outCfg}) or die "Copy failed: $config{mainCfg} -> $config{outCfg}";

# parse input files
my $mainBs = &parseG192 ($config{mainBs});
my $helpBs = &parseG192 ($config{helpBs});
my $epf    = &parseEPF ($config{epf});

# extend epf pattern if required
my $epf_initial;
@$epf_initial = @$epf;
while (@$epf < @$mainBs) {
        push @$epf, @$epf_initial;
}

my $burst_index = 0;
my @newBs;
my %stat;
$stat{main} = 0;
$stat{help} = 0;
$stat{nodata} = 0;

my @tmp;
# assemble new bit stream
for (my $f=0; $f<@$mainBs; $f++) {
    my $newFrame;
    my $bfi = $epf->[$f + $ep_offset];
    if ($bfi==0) {                       # no packet loss, take main payload
        $burst_index=0;
        push @newBs, $mainBs->[$f];
        $stat{main}++;
        push @tmp, 0;
    } else {
        $burst_index++;
        if ($burst_index <= $config{offset}) {           # packet loss, but secondary payload available
            $helpBs->[$f]->{sync} += ($config{signal}==3);
            push @newBs, $helpBs->[$f];
            $stat{help}++;
            push @tmp, $config{signal};
        } else {
            $mainBs->[$f]->{sync} -= 1;                  # packet loss, no secondary payload available
            $mainBs->[$f]->{data} = "";
            for (my $i=0; $i<$mainBs->[$f]->{length}; $i++) {
                $mainBs->[$f]->{data} .= (pack "v", 0);
            }
            push @newBs, $mainBs->[$f];
            $stat{nodata}++;
            push @tmp, 1;
        }
    }
}

print "Stats: main $stat{main}, help $stat{help}, nodata $stat{nodata}\n";

# write new stream
&writeG192(\@newBs,$config{output});

### END


sub writeG192 {
    my $bs = shift;
    my $file = shift;

    open FILE, ">$file" or die $!;
    binmode FILE;

    foreach my $frame (@$bs) {
        print FILE (pack "v", $frame->{sync});
        print FILE (pack "v", $frame->{length});
        print FILE $frame->{data};
    }

    close FILE;
}


sub parseCommandline {

    my $i=0;
    $config{mainBs} = $ARGV[$i++];
    $config{helpBs} = $ARGV[$i++];
    $config{epf}    = $ARGV[$i++];
    $config{output} = $ARGV[$i++];
    $config{offset} = $ARGV[$i++] * $framing/10;
    $config{signal} = $ARGV[$i++];
}

sub printHelp {
  print <<FOO;

/******************************************************************************
*                        ETSI TS 103 634 V1.1.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

Assembles new bit stream out of primary and secondary LC3plus frames based
on an error pattern file. Secondary frames are inserted, if no primary data
available.\n
Usage: $0 primary_bitstream.g192 secondary_bitstream.g192 error_pattern output_bs.192 offset signal\n
  where
     primary_bitstream.g192     primary LC3plus bitstream; g192 format
     secondary_bitstream.g192   secondary LC3plus bitstream; g192 format
     error_pattern              16bit or ITU format error error_pattern
     output_bs                  resulting bit bitstream
     offset                     temporal offset of the secondary frame compared to
                                prinary frame in 20ms packets
     signal                     singaling for secondary frames, should be 3

FOO
    exit;
}

sub parseG192 {
        my $input = shift;
        my @bs_data;

        open FILE, $input or die $!;
        binmode FILE;
        my $eof= 0;
        my ($sync, $tsync, $length, $tlen, $data);
        do {
            $eof = read (FILE, $sync, 2);
            if ($eof > 0) {
                my %frame_data;
                $tsync = (unpack "v", $sync);
                $frame_data{sync} = $tsync;
                if ($tsync == int (0x6b21)) {
                    $eof = read (FILE, $length, 2);
                    $tlen = unpack "v", $length;
                    $frame_data{length} = $tlen;
                    read (FILE, $data, $tlen*2);
                    $frame_data{data} = $data;
                    push @bs_data, \%frame_data;
                }
            }
        }
        while ( $eof );
        close (FILE);


        return \@bs_data;
}

sub parseEPF {
    my $input = shift;
    my @ep_data;

    # 16 bit epf file
    if (-B $input) {
        open FILE, $input or die $!;
        binmode FILE;
        my $eof= 0;
        my ($v, $t);
        do {
            $eof = read (FILE, $v, 2);
            if ($eof > 0) {
                $t = unpack "v", $v;
                push @ep_data, $t;
            }

        }
        while ( $eof );
        close FILE;
    } else {
        # .g192 file
        open FILE, $input or die $!;
        my $tmp = <FILE>;
        my @data = split("k", $tmp);
        foreach (@data) {
            if ($_ eq "!") {
                push @ep_data, 0;
            } else {
                push @ep_data, 1;
            }
        }
        close FILE;
    }

    return \@ep_data;
}
