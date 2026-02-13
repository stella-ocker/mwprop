#!/usr/bin/perl
# 15 May 2002
# General version: command-line input includes field wanted 
$nargs = ($#ARGV + 1);
if($nargs < 5) {
   print "Usage:\n"; 
   print "	run_NE2001      l    b       DM/D         ndir    field \n"; 
   print "	               deg  deg    pc-cm^{-3}     1,-1    D etc \n";
   print "	                           or kpc \n";
   print "	Possible Fields (case insensitive):\n";
   print "	Dist, DM, SM, EM, TAU, SBW, SCINTIME, THETA_G, THETA_X, NU_T, ALL\n";
   exit 0;
}
if($nargs == 5) {
   $l = $ARGV[0];
   $b = $ARGV[1];
   $DMD = $ARGV[2];
   $ndir = $ARGV[3]; 
   $field = $ARGV[4];
   $FIELD = uc($field);
}
#print "$field $FIELD\n";
$bindir = "./run_NE2001/";
open(TEMP02, "$bindir/NE2001 $l $b  $DMD  $ndir | ") || die "Couldn't run NE2001: $!\n";
while (<TEMP02>) {
	if ($field eq "ALL" || $field eq "all") {print};
	if (/\s+$field\s+/ || /\s+$FIELD\s+/ ) {
		($value, $name, $units, $description)=split(' ',$_, 4); 
		print "$name = $value $units         $description";
	}
}
close(TEMP02);
exit 0;

