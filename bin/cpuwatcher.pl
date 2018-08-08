#!/usr/bin/perl -w
#
#
sub usage
{
   return <<THE_END;

 usage $0 user min_number_of_cpus mailadress dt

 where dt is how often (in seconds) we will check the top list.

 for some stupid reason it does not work on linux...

THE_END
}

sub message
{
   print STDERR "$0: ",@_, "\n";
}
use constant SAFETY_LIMIT => 2;

my $mailprog = 'mail';
my $subject   = "Subject: cpuwarning\n";

my $user     = shift || die usage();
my $min_cpus = shift || die usage();
my $mailadr  = shift || die usage();
my $dt       = shift || die usage();



my $count = 0;
while (1) {
   my $cmd = sprintf("top 16 | grep $user");
   my @toplist = `$cmd`;
#print "cmd= $cmd\n";
#print @toplist;
   message("number of entries on top-list for user $user = " . @toplist);

   if (@toplist < $min_cpus) {
      ++$count;
      message("used cpus below level ($min_cpus) for the $count time");
   }

   # will only send mail when we are sure there is a problem ...
   if ($count > SAFETY_LIMIT) {
      message("sending mail to $mailadr");

      open(MAIL,"|$mailprog $mailadr") || die "cannot open mailprogram $mailprog";
      print MAIL $subject;
      print MAIL "\nNumber of cpus used by $user is now " . @toplist;
      close MAIL;

      # stop looping. dont want a lot of mails generated.
      last;
   }
   sleep($dt);
}
