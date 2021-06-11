#!/usr/bin/env perl
use warnings;
#This script built a figure for polishing results for the SLAG paper.
open (RSC, "> $ARGV[0]"); #adhoc06082020a.R
$pdffile = $ARGV[1]; #polishingresults06082020a.svg
$pdfheight = $ARGV[2]; #3
$pdfwidth = $ARGV[3]; #4
$ubound = $ARGV[4]; #3550
$lbound = $ARGV[5]; #3090
$highbound = $ARGV[6]; #98
$lowbound = $ARGV[7]; #90
print RSC "x <- c(0, 1, 2, 3, 4, 5)\n";
#Convert to percentages.
print RSC "meanidentities <- c(90.1137, 97.8266, 96.9667, 96.5518, 96.3501, 96.2017)\n";
print RSC "meanhitlengths <- c(3543.66, 3402.79, 3222.35, 3183.72, 3124.77, 3101.04)\n";
print RSC "meancontiglengths <- c(3523.49, 3498.29, 3365.68, 3350.77, 3333.07, 3326.28)\n";
if ($pdffile =~ m/pdf$/) {print RSC ")\npdf(\"$pdffile\", height = $pdfheight, width = $pdfwidth)\n";}
elsif ($pdffile =~ m/svg$/) {print RSC "svg(\"$pdffile\", height = $pdfheight, width = $pdfwidth)\n";}
else {print RSC ")\npng(\"$pdffile\", height = $pdfheight, width = $pdfwidth, units = \"in\")\n";}
print RSC "par(mar = c(5, 5, 3, 5))\n";
print RSC "plot(x, meancontiglengths, ylim = c($lbound, $ubound), type = \"l\", xlab = \"Times polished\", ylab = \"Length\", lty = 1, col = \"black\")\n";
print RSC "lines(x, meanhitlengths, lty = 2, col = \"black\")\n";
print RSC "par(new = TRUE)\n";
print RSC "plot(x, meanidentities, ylim = c($lowbound, $highbound),type = \"l\", xlab = \"\", ylab = \"\", xaxt = \"n\", yaxt = \"n\", lty = 3, col = \"black\")\n";
print RSC "axis(side = 4)\nmtext(\"Percent identity\", side = 4, line = 3)\ndev.off()\n";
close (RSC);

#print RSC "plot(x, rallcontigs, ylim = c($lowbound, $highbound), type = \"l\", xlab = \"\", ylab = \"\", xaxt = \"n\", yaxt = \"n\", lty = 3, col = \"black\")\n";
#print RSC "plot(x, rmatchinglongests, ylim = c(0, $ubound), type = \"l\", xlab = \"Cycle\", ylab = \"Length\", lty = 1, col = \"blue\")\n";
#print RSC "lines(x, rmeanmatchinglengths, lty = 2, col = \"blue\")\n";
#print RSC "lines(x, rmeanalllengths, lty = 2, col = \"black\")\n";
