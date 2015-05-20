#! usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Data::Dumper;
$|++; # command buffering

# COMMAND LINE ARGUMENTS

# CREATE a LOG FILE
# as it is a typeglob use \*FILEHANDLE to pass it through functions
# also it could be just LOG not a new variable $LOG
print "\033[0;33m[==========]\033[0m CREATING LOG FILE ...\n\n";
open(my $LOG, "> log.txt") || die "### ERROR ### Cannot open file: log.txt\n";
print "\033[0;32m[       OK ]\033[0m DONE\n\n";


# VARIABLE DEFINITIONS
print "\033[0;33m[==========]\033[0m DEFINING VARIABLES ...\n\n";
my %FILTERED_VCF      = ();
print "\033[0;32m[       OK ]\033[0m DONE\n\n";

# FUNCTION CALLS
print "\033[0;33m[==========]\033[0m CALLING FUNCTIONS ...\n\n";
&filter_create_hash_bcf_files(\%FILTERED_VCF, 1, \*$LOG);
exit(0);

# FUNCTION DEFINITION
sub filter_create_hash_bcf_files($$$$){
  my($filteredVcf, $DEBUG, $LOG) = @_;

  print "\t\t\033[0;33m[==========]\033[0m FUNCTION filter_create_hash_bcf_files WORKING ...\n\n";

  # Open directory
  # If you put this script inside the directory where you have all the .bcf files
  # use "./"
  my $dir = "/home/sergimash/NGS_sequencing/samtools-1.1_WORKING_WELL/bcf_final/";
  print "\t\t\033[0;33m[==========]\033[0m OPENING DIRECTORY $dir ...\n\n";
  opendir(BCF, $dir) || die "### ERROR ### Cannot open directory: $dir\n\n";
  print "\t\t\033[0;32m[       OK ]\033[0m DONE\n\n";

  # Read each file in directory
  print "\t\t\033[0;33m[==========]\033[0m READING FILES ...\n\n";
  my $files = 0; # start file counter
  READDIR: while(my $bcf = readdir(BCF)){
    # We are using .bcf file (which is a binary version of .vcf)
    next READDIR if $bcf !~ /\.bcf$/; # skip files that DO NOT end in .bcf
    # Uncomment the below line for .vcf filtering
    # next READDIR if $bcf !~ /\.vcf$/;

<<OPTIONS;
         -Q INT    minimum RMS mapping quality for SNPs [\$opts{Q}]
         -d INT    minimum read depth [\$opts{d}]
         -D INT    maximum read depth [\$opts{D}]
         -a INT    minimum number of alternate bases [\$opts{a}]
         -w INT    SNP within INT bp around a gap to be filtered [\$opts{w}]
         -W INT    window size for filtering adjacent gaps [\$opts{W}]
         -1 FLOAT  min P-value for strand bias (given PV4) [\$opts{1}]
         -2 FLOAT  min P-value for baseQ bias [\$opts{2}]
         -3 FLOAT  min P-value for mapQ bias [\$opts{3}]
         -4 FLOAT  min P-value for end distance bias [\$opts{4}]
         -e FLOAT  min P-value for HWE (plus F<0) [\$opts{e}]
         -p        print filtered variants

         DEFAULT VALUES
         -Q 10
         -d 2
         -D 10000000
         -a 2
         -w 3
         -W 10
         -1 1e-4
         -2 1e-100
         -3 0
         -4 1e-4
         -e 1e-4
         -p (either put '-p' on the command line or don't put anything)
OPTIONS

    # Here is using the DEFAULT values for each filter, but you can change them accordingly
    my $cmd = "bcftools view $dir$bcf | /home/sergimash/NGS_sequencing/bcftools-1.2/vcfutils.pl varFilter -Q 10 -d 2 -D 10000000 -a 2 -w 3 -W 10 -1 1e-4 -2 1e-100 -3 0 -4 1e-4 -e 1e-4 > $bcf.FILTERED.vcf";
    system($cmd);

    # Open each .bcf file
    open(vcfFIL, "< $bcf.FILTERED.vcf") || die "### ERROR ### Cannot open file: $bcf.FILTERED.vcf\n";

    $files++; # increment file counter by 1 wicth each file read

    # READ each .bcf file from the directory
    READvcf: while(<vcfFIL>){
      # FIELDS in .bcf file
      # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	GENOTYPE
      next READvcf if /^\#/o; # skip line that are comments (start with #)
      next READvcf if /^\s*$/o; # skip blank lines
      chomp; # remove \n character. In windows/files that are creates there
             # you might want to do $_ =~ s/\n|\n\r|\r\n|\r//o;

      my($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $genotype);

      ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $genotype) = split /\t/o,$_,10;

      exists($filteredVcf->{$pos}) || do {
        $filteredVcf->{$pos} = {
          'CHROM'    => $chrom,
          'POS'      => $pos, # again I know
          'ID'       => $id,
          'FILES'    => {}
        };
      };

      # IF each .bcf file is a SAMPLE/PATIENT/INDIVIDUAL is useful that
      # foreach position to get how many SAMPLE/PATIENT/INDIVIDUAL have that position as a
      # variant
      exists($filteredVcf->{$pos}{'FILES'}{$bcf}) || do {
        $filteredVcf->{$pos}{'FILES'}{$bcf} = {
          'REF'      => $ref,
          'ALT'      => $alt,
          'QUAL'     => $qual,
          'INFO'     => $info,
          'FORMAT'   => $format,
          'GENOTYPE' => $genotype
        };
      };

    }; # while READvcf

    # Cleanup. Remove the generated files after creating the hash
    system "rm $bcf.FILTERED.vcf";

    # Here you can do all sort of stats or operations will the data obtained in the hash.
    # For example, for each position how many 'FILES' have it as a variant


    # exists()
  }; # while READDIR
  if($files == 0){
    print "\t\t\033[0;31m[   ERROR  ]\033[0m NO .bcf/.vcf FILES ...\n\n";
  };

  print "\t\t\033[0;32m[       OK ]\033[0m DONE\n\n";
  print "\t\t\033[0;33m[==========]\033[0m CLEANING UP ...\n\n";
  print "\t\t\033[0;32m[       OK ]\033[0m DONE\n\n";
  print "\t\t\033[0;33m[==========]\033[0m PERFORMING CALCULATIONS ...\n\n";
  foreach my $pos (keys %{$filteredVcf}){
    $filteredVcf->{$pos}{'Nº_files'} = scalar keys $filteredVcf->{$pos}{'FILES'};
  };
  print "\t\t\033[0;32m[       OK ]\033[0m DONE\n\n";
  print "\t\t\033[0;33m[==========]\033[0m READ $files files from $dir ###\n\n";
  print "\t\t\033[0;33m[==========]\033[0m DUMPING DATA INTO LOG ...\n\n";
  print $LOG Data::Dumper->Dump([$filteredVcf],[qw(*FILTERED_VCFs)]), "\n" if $DEBUG;
  print "\t\t\033[0;32m[       OK ]\033[0m DONE\n\n";
};
