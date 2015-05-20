This is a Perl script useful for reading multiple files (.bcf/.vcf) in order to
filter them using the __bcftools__ and __vcfutils.pl__ from the same tool.

Get [bcftools](https://github.com/samtools/bcftools).

### INSTRUCTIONS:

1. install bcftools and change directory 'cd' to the install directory
2. write 'pwd' on the command line. Copy the result
3. replace that directory in the variable **__'$cmd'__** of the script
after the pipe '|' like: '/home/user/bcftools/vcfutils.pl'  

### USAGE:

perl fil_hash_bcfs.pl

### OUTPUT

log.txt (dumping of the hash)
