#!/bin/bash

#
# server_info.sh - display server hardware info
#
# 2008 - Mike Golvach - eggi@comcast.net
#
# Creative Commons Attribution-Noncommercial-Share Alike 3.0 United States License
# 

hwinfo="/usr/sbin/hwinfo --short"
hostname=`hostname`
separator="----------------------------------------"
echo $separator
echo "System Information For $hostname"
echo $separator
echo $separator
echo SERVER - MEMORY
echo $separator
/usr/sbin/hwinfo --bios|egrep 'OEM id:|Product id:|CPUs|Product:|Serial:|Physical Memory Array:|Max. Size:|Memory Device:|Location:|Size:|Speed:|Location:'|sed -e 's/"//g' -e '/^ *Speed: */s/Memory Device:/\n  Memory Device:/' -e 's/\(Max. Speed:\)/CPU \1 MHz/' -e 's/\(Current Speed\)/CPU \1 MHz/'
echo $separator
echo SMP
echo $separator
$hwinfo --smp
echo $separator
echo CPU
echo $separator
$hwinfo --cpu
echo $separator
echo CD_ROM
echo $separator
/usr/sbin/hwinfo --cdrom|egrep '24:|Device File:|Driver:'|awk -F":" '{ if ( $1 ~ /[0-9][0-9]*/ ) print $0; else print "  " $2}'|sed -e 's/^.*[0-9] //' -e 's/ //' -e 's/"//g'
echo $separator
echo DISK
echo $separator
$hwinfo --disk
echo $separator
echo PARTITION
echo $separator
$hwinfo --partition
echo $separator
echo NETWORK
echo $separator
$hwinfo --network
echo $separator
echo NETCARD
echo $separator
$hwinfo --netcard
echo $separator
