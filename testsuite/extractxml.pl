#!/usr/bin/perl
#
# This file is part of the SYMPLER package.
# https://github.com/kauzlari/sympler
#
# Copyright 2002-2013, 
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

#This script extracts the tag names (modules) from sympler input files and stores the entries in a spreadsheet.
#Then it checks which module is tested in which input file and displays them in spreadsheet. 


use Spreadsheet::ParseExcel;
use Spreadsheet::WriteExcel;
use Spreadsheet::ParseExcel::SaveParser;
use FindBin '$Bin';
use Getopt::Long;
use Cwd;

my $sheetdir= "$Bin";	#default directory for spreadsheet 
GetOptions(
          "sheetdir=s"           =>      \$sheetdir,	#here takes a directory as input for storing spreadsheet
          );

$sheetdir= Cwd::abs_path($sheetdir); 	#to get absolute path of directory

# create a build directory 
if (-d $sheetdir) 
{
	chdir "$sheetdir";
}
else 
{
	mkdir "$sheetdir" or die "$!\n"; 
	chdir "$sheetdir";
}
      
#Extracting Tag values from xml file	
opendir(TSTIN_DIR, "$Bin/TSTIN/");
@dirs = grep { $_ ne '.' && $_ ne '..' && $_ ne '.svn' && $_ ne '.git'} readdir TSTIN_DIR;
closedir TSTIN_DIR;
foreach $dir (@dirs) 
{
 	open(XMLIN, "$Bin/TSTIN/$dir/test.xml");
 	$line = <XMLIN>;
 	open(OUT, ">>/tmp/out.txt");
 	open(OUTDIR, ">>/tmp/out_$dir.txt");
 	while($line ne "") 
	{
		if ($line=~ m/<\s*(\w+)/)
		{
			print OUT "$1\n" ;
			print OUTDIR "$1\n" ;
		}
	$line = <XMLIN>;
	}
	close OUT;
	close OUTDIR;
}
    
system("sort -d /tmp/out.txt > /tmp/sorted.txt");

#Creating desired spreadsheet
 
# Create a new spreadsheet workbook called XMUNLCHECK.xls for temporary usage
my $workbook = Spreadsheet::WriteExcel->new('/tmp/XMLUNCHECK.xls');

#Add worksheet TAGS
my $worksheet = $workbook->add_worksheet("TAGS");
my $row=1;
#Add format for Entries      
$format = $workbook->add_format();
$format->set_size(10);
$format->set_color('blue');
$format->set_align('justify');
$format->set_align('vjustify');        
        
$format1 = $workbook->add_format();
$format1->set_bold();
$format1->set_size(10);
$format1->set_color('purple');
$format1->set_align('justify');
$format1->set_align('vjustify');        

#Set column width
$worksheet->set_column(0, 0, 30);
$worksheet->set_column(1, 3, 20);

#Adding coloumn headers
$worksheet->write(0, 0, "Module", $format1);
opendir(TSTIN_DIR, "$Bin/TSTIN/");
@dirs = grep { $_ ne '.' && $_ ne '..' && $_ ne '.svn' && $_ ne '.git'} readdir TSTIN_DIR;
closedir TSTIN_DIR;
my $col=1;
foreach $dir (@dirs) 
{
	my $row =0;
        $worksheet->write($row, $col, "$dir.xml", $format1);
        $col++;
}
         
# Adding row headers
my $file = '/tmp/sorted.txt';
my %seen = ();
{
	local @ARGV = ($file);
        while(<>)
	{
	        $seen{$_}++;
	        next if $seen{$_} > 1;
        	my $col =0;
                #Adding Entries to the sheet
           	$worksheet->write($row, $col, $_, $format);
	   	$row++;
        }
}
        
system("rm /tmp/out.txt");
system("rm /tmp/sorted.txt");
$workbook->close();	

#Marking checked modules

my $FILE = "/tmp/XMLUNCHECK.xls";
my $SHEETNAME = "TAGS";
my $parser   = Spreadsheet::ParseExcel::SaveParser->new();
my $excel = $parser->Parse($FILE);
my $sheet = $excel->Worksheet($SHEETNAME);

$format2->{AlignH};


opendir(TSTIN_DIR, "$Bin/TSTIN/");
@dirs = grep { $_ ne '.' && $_ ne '..' && $_ ne '.svn' && $_ ne '.git'} readdir TSTIN_DIR;
closedir TSTIN_DIR;
foreach $dir (@dirs) 
{
	open (TESTF,"/tmp/out_$dir.txt"); 
	my $line = <TESTF>;
	while($line ne "")
	{
  		foreach my $row ($sheet->{MinRow} .. $sheet->{MaxRow}) 
		{
	                my $cell = $sheet->{Cells}[$row][0];
                        if ($cell ->{Val} eq $line) 
			{
	                        foreach my $col ($sheet->{MinCol} ..  $sheet->{MaxCol})
				{    
	                        	my $cell = $sheet->{Cells}[0][$col];
	                        	if ($cell ->{Val} eq "$dir.xml")
	                        	{ 
	                        		$sheet->AddCell( $row, $col ,'X', $format2 );	
	                        		$excel->SaveAs("$sheetdir/XMLCHECKED.ods"); #Here saves the final spreadsheet(in ods format) in desired directory
	                        		$col++;
	                    		}
	                	}
	        	}
		}
		$line = <TESTF>;
	}
	close TESTF;
	system("rm /tmp/out_$dir.txt");
}


