#! /usr/bin/perl
#
# This file is part of the SYMPLER package.
# https://github.com/kauzlari/sympler
#
# Copyright 2002-2017, 
# David Kauzlaric <david.kauzlaric@frias.uni-freiburg.de>,
# and others authors stated in the AUTHORS file in the top-level 
# source directory.
#
# SYMPLER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SYMPLER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SYMPLER. If not, see <http://www.gnu.org/licenses/>.
#
# Please cite the research papers on SYMPLER in your own publications. 
# Check out the PUBLICATIONS file in the top-level source directory.
#
# You are very welcome to contribute extensions to the code. Please do 
# so by making a pull request on https://github.com/kauzlari/sympler
# 
#
use POSIX;
use Getopt::Long;
use Cwd;

#An automated test script to compile, run and debug the Sympler code	

$pwd=`pwd`;
use FindBin '$Bin';

if($pwd !~ $Bin)
{
	print "\nPlease make sure that you are in '$Bin/' directory while running this test script\n\n";
	exit;
}

my $builddir= "/tmp/BUILDDIR2";	#default directory for build 
GetOptions(
          "builddir=s"           =>      \$builddir,	#here takes a directory as input for build	
          );

$builddir= Cwd::abs_path($builddir); 	#to get absolute path of directory
my @builddir = split('/', $builddir); 	#creates an array of Dir names
pop @builddir; 	#removes the currnet directory 
$build_par = join( '/', @builddir);	 #Joins names back into parent path
chdir "$build_par";                #change to parent directory
$findsvn = system("find .svn >/tmp/findres.txt 2>&1"); #check wheather it contains .svn directoy, $find returns 0 if it contains
$findgit = system("find .git >/tmp/findres.txt 2>&1"); #check wheather it contains .svn directoy, $find returns 0 if it contains
if($findsvn eq 0)
{
	print "\nCannot create directory here because under subversion control! Choose another\n\n";
	exit;
}
if($findgit eq 0)
{
	print "\nCannot create directory here because under git control! Choose another\n\n";
	exit;
}

# create a build directory 
if (-d $builddir) 
{

	chdir "$builddir";
}
else 
{
	mkdir "$builddir" or die "$!\n"; 
	chdir "$builddir";
}
#Compile the code
print"CONFIGURATION PROGRESSING...\n";
$com= `$Bin/../configure --with-tnt --with-superlu CXXFLAGS="-O3 -g" > configout.txt 2>&1`; 
if($com == 0)
{
	print"CONFIGURATION DONE\n";
}
else
{
	print"CONFIGURATION FAILED!!!\n";
	exit;
}
print"COMPILATION PROGRESSING...\n";
$num_cores=`grep 'core id' /proc/cpuinfo | sort -u | wc -l`; #count #of cores on cpu
$numcores  = floor($num_cores); 
$N=$numcores*3;
$ret=system("make -j$N > complout.txt 2>&1");

#Checking whether compilation successful
if($ret == 0) 
{ 
	print "COMPILATION DONE SUCCESSFULLY\n"; 
}
else 
{ 
	print "COMPILATION FAILED!!!\n";
	print "For more information read '$builddir/complout.txt'  \n";
	exit;
}
system "pwd";
# create a temporary directory for Simulation
if(-d "/tmp/SIMULATION")
{
     	system('rm','-r',"/tmp/SIMULATION");
     	mkdir "/tmp/SIMULATION" or die "$!\n"; 	
	chdir "/tmp/SIMULATION";
}
else
{
	mkdir "/tmp/SIMULATION" or die "$!\n"; 	
	chdir "/tmp/SIMULATION";
}
system "pwd";	
# running all test input files against sympler 
opendir(TSTIN_DIR, "$Bin/TSTIN/");
@dirs = grep { $_ ne '.' && $_ ne '..' && $_ ne '.svn'} readdir TSTIN_DIR;
closedir TSTIN_DIR;
foreach $dir (@dirs) 
{
	system("cp $Bin/TSTIN/$dir/* /tmp/SIMULATION/");
	eval {
 	  local %SIG;
 	  $SIG{ALRM}=
 	  sub{ die "Timeout, test input '$dir.xml' took more than 2 minutes to run!\n"; };
	  alarm 120;
  	  $test=system("$builddir/src/sympler $Bin/TSTIN/$dir/test.xml > /tmp/test_$dir.txt 2>&1");
	  alarm 0;
       	 };

	alarm 0;

	if($@) { print "Error: $@\n"; }

 	if($test)
	{  
		print "Error Occured while running input file $Bin/TSTIN/$dir/test.xml !!! Check '/tmp/test_$dir.txt' \n";
	} 
	
	#compare produced result in '/tmp/SIMULATION/' directory with backuped one 'testsuite/RESULTS/' directory
	opendir RES, "$Bin/RESULTS/" or die "$!";
	my @res= grep { $_ ne '.' && $_ ne '..' && $_ ne '.svn' && $_ ne '.git'} readdir RES;
	closedir RES;

	opendir SIM, "/tmp/SIMULATION/" or die "$!";
	my @sim= grep { $_ ne '.' && $_ ne '..' } readdir SIM;
	closedir SIM;
	

	foreach my $res(@res)
	{
		foreach my $sim(@sim)
		{
			system("perl -pi -e 's/-0.000000/0.000000/g' $sim");
			if($res eq $sim)
			{	
				$ret=system("diff -E -b -t $Bin/RESULTS/$res /tmp/SIMULATION/$sim >/tmp/diff_$sim");
				if ($ret){print "RESULT CHANGED FOR $res and $sim AND CHANGED RESULTS ARE STORED IN '/tmp/diff_$sim \n";}
				else {print "RESULT IS SAME FOR $res and $sim\n";} 
			}
			
		}
	}

#delete contents of SIMULATION
	system ("rm -rf /tmp/SIMULATION/*");

#end foreach
}

#delete SIMULATION after Simulation 
system('rm','-r',"/tmp/SIMULATION");
print("---------END OF SIMULATION----------\n");

#####END OF SIMULATION#####
