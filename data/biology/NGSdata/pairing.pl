#! /usr/bin/perl -w

###################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html)
#
# Intellectual property belongs to IRD and SouthGreen developpement plateform
#
# Written by Francois Sabot 
#
#################### 


use strict;
use Getopt::Long;
use lib qw(/home/sabotf/bin/);
use alldays qw(all);


my $courriel="francois.sabot-at-ird.fr";
my ($nomprog) = $0 =~/([^\/]+)$/;
my $MessAbruti ="\nUsage:
\t$nomprog -f forward_seq -r reverse_seq 
or
\t$nomprog -forward forward_seq -reverse reverse_seq

DESCRIPTION

From a forward and a reverse fastq files unpaired, it will generate three new files forward.fastq, reverse.fastq and single.fastq in the current folder
The forwad.fastq and reverse.fastq files will contain the paired sequences, and the single.fastq the remaining non-paired ones.
The sequences will be renamed using

	contact: $courriel\n\n";
	

unless (@ARGV) 
	{
	print "\nType --help for more informations\n\n";
	exit;
	}

my ($forward, $reverse,$help);

GetOptions("help|?|h" => \$help,	
	   "f|forward=s"=>\$forward,
           "r|reverse=s"=>\$reverse);
			
if ($help){print $MessAbruti; exit;}

#Checking rights for reading and availability of files
alldays::read_rights($forward);
alldays::read_rights($reverse);

#Opening inputs and outputs
open IN, $forward;
open IN2,$reverse;

open MATEF, ">forward.fastq";
open MATER, ">reverse.fastq";
open SINGLE, ">single.fastq";

#Picking up seq ID from the forward file
my %hashforward;
while (<IN>)
	{
	my $line = $_;
	chomp $line;
        next if ($line =~ m/^$/);
        my $next = $line."\n";
	$line =~ s/\/\d$//;
        $line =~ s/\s\d:\w:\d:\w{1,10}$//;
	$next .= <IN>;
	$next .= <IN>;
	$next .= <IN>;
	$hashforward{$line}=$next;
	}
close IN;

#Comparing with the reverse seq ID
while (<IN2>)
	{
	my $line = $_;
	chomp $line;
        next if ($line =~ m/^$/);
        my $next = $line."\n";
	$line =~ s/\/\d$//;
        $line =~ s/\s\d:\w:\d:\w{1,10}$//;
	$next .= <IN2>;
	$next .= <IN2>;
	$next .= <IN2>; 
	#printing outputs for paired files
	if (exists $hashforward{$line})
		{
		my $out = $hashforward{$line};
		print MATEF $out;
		my $out2 = $next;
		print MATER $out2;
		delete $hashforward{$line}; # To save memory and to conserve only the singles
		}
        #printing output for singles from reverse file
	else
		{
		my $out2 = $next;
		print SINGLE $out2;
		}
	}

#printing output for singles from forward file
foreach my $remain (keys %hashforward)
	{
	my $out = $hashforward{$remain};
	print SINGLE $out;
	}

#Closing files
close IN2;
close SINGLE;
close MATEF;
close MATER;

exit;

__END__