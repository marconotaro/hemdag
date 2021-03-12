#!/usr/bin/perl

use strict;
use warnings;

print "#!/bin/sh\n\n";
print "tot_start=\$(date +%s)\n\n";

my @orgs=  qw(7227.drome); 
my @flats= qw(svmlinear);  
my @algs=  qw(gpav isotprTF isotprW isodescensTF isodescensW isodescensTAU);
my @onts=  qw(mf); ## qw(bp mf cc)

my $m=12; ## block of parallel evaluations to be computed before going to the next block
my $k=0;  ## cpu number on which binding a task 
foreach my $org (@orgs){
    foreach my $flat (@flats){
        foreach my $alg (@algs){
            foreach my $ont (@onts){
                $k++;
                my $cpu= $k-1;
                print "taskset -c $cpu Rscript hemdag-perf-eval.R -o $org -d $ont -f $flat -a $alg > $org.$flat.$alg.go.$ont.perfmeas.out 2> /dev/null &\n";
                if($k % $m ==0){
                    print "\n";
                    print "wait\n";
                    print "\n";
                    $k=0;
                } 
            }
        }
    }
}

print "\n";
print "tot_end=\$(date +%s)\n";
print "tot_elapsed_s=\$((tot_end-tot_start))\n";
print "tot_elapsed_m=\$((tot_elapsed_s/60))\n";
print "tot_elapsed_h=\$((tot_elapsed_m/60))\n";
print "printf \"grand total elapsed time:\t\$((tot_elapsed_s)) SECONDS\t\$((tot_elapsed_m)) MINUTES\t\$((tot_elapsed_h)) HOURS\"\n";
print "echo\n\n";

print "echo \"compute performance done\"\n\n";


