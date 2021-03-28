#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long 'HelpMessage';

GetOptions(
    "org=s"  => \(my @orgs=""),
    "flat=s" => \(my @flats=""),
    "alg=s"  => \(my @algs=""),
    "dmn=s"  => \(my @onts=""),
    'exp=s'  => \(my $exp='ho'),
    'chk=i'  => \(my $chk='12'),
    'help'   => sub {HelpMessage(0)},
) or HelpMessage(1);

@orgs=  split(/,/,join(',',@orgs));
@flats= split(/,/,join(',',@flats));
@algs=  split(/,/,join(',',@algs));
@onts=  split(/,/,join(',',@onts));

print "#!/bin/sh\n\n";
print "tot_start=\$(date +%s)\n\n";

my $k=0;  ## cpu number on which binding a task
foreach my $org (@orgs){
    foreach my $flat (@flats){
        foreach my $alg (@algs){
            foreach my $ont (@onts){
                $k++;
                my $cpu= $k-1;
                print "taskset -c $cpu Rscript hemdag-eval.R -o $org -d $ont -e $exp -f $flat -a $alg > $org"."_go_"."$ont"."_"."$flat"."_"."$alg"."_perfmeas.out 2> /dev/null &\n";
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

call-hemdag-chunk - generate HEMDAG evaluation calls chunks

=head1 SYNOPSIS

    --org,-o    organism list. The list can have one or more elements. The elements must be comma separated (7227_drome,9031_chick,..)
    --flt,-f    flat classifier list. The list can have one or more elements. The elements must be comma separated (svmlinear,ranger,..)
    --alg,-a    hierarchical correction algorithm list. The list can have one or more elements. The elements must be comma separated (gpav,isodescensTF,isodescensTAU,..)
    --dmn,-d    gene ontology domain list. The list can have one or more elements. The elements must be comma separated (bp,mf,cc)
    --exp,-e    type of dataset on which evaluate HEMDAG. It can be: ho or cv (def. ho)
    --chk,-c    number of parallel evaluations to be computed before going to the next block (def 12)
    --help,-h   print this help

=cut

