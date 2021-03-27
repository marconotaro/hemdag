#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long 'HelpMessage';

GetOptions(
    'exp=s' => \(my $exp='ho'),
    'chk=i' => \(my $chk='12'),
    'help'  => sub {HelpMessage(0)},
) or HelpMessage(1);

print "#!/bin/sh\n\n";
print "tot_start=\$(date +%s)\n\n";

my @orgs=  qw(7227_drome); ## organism(s) list
my @flats= qw(svmlinear);  ## flat classifier(s) list
my @algs=  qw(gpav isotprTF isotprW isodescensTF isodescensW isodescensTAU);  ## HEMDAG algorithm(s) list
my @onts=  qw(mf);  ## GO domain(s) list (bp mf cc)

my $k=0;  ## cpu number on which binding a task
foreach my $org (@orgs){
    foreach my $flat (@flats){
        foreach my $alg (@algs){
            foreach my $ont (@onts){
                $k++;
                my $cpu= $k-1;
                print "taskset -c $cpu Rscript hemdag-eval.R -o $org -d $ont -e $exp -f $flat -a $alg > $org"."_go_"."$ont"."_"."$flat"."_"."$alg"."_perfmeas.out 2> /dev/null \n";
                if($k % $chk ==0){
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

__END__

=pod

=head1 NAME

call-hemdag-per-eval - generate HEMDAG evaluation calls chunks

=head1 SYNOPSIS

  --exp,-e   type of dataset on which evaluate HEMDAG. It can be: ho or cv (by def. ho)
  --blk,-b   number of parallel evaluations to be computed before going to the next block (by def 12)
  --help,-h  print this help

=cut

