my $rc = -1;
my $cmd;
#$rc = qx(ping por017lin.hre.hydro.com);
#$rc = qx(ls);
#$rc = `ls > /dev/null`;

$cmd = 'ping -c2 por017lin.hre.hydro.com > /dev/null';
#$cmd = 'ls > /dev/null';
$rc = system($cmd);

print "cmd= $cmd rc= $rc\n";
