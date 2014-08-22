#!/usr/local/bin/perl
#
# This script takes two argument:
#	1- list of path 
#	2- list of library names
#
# It is assumed that they are ordered and matching. 
# All the script does is generate a list of fully qualified
# library names.

($dir,$lib) = @ARGV;

# Separate directories from libraries
@dirlist = split(/ /, $dir);
@liblist = split(/ /, $lib);

# Invert library list 
while ($l = pop @liblist) {
	push @libstack, $l;
}

# Match them
foreach $d (@dirlist) {
	$l = pop @libstack;
	push @result, "$d/$l";
}

print "@result\n";

