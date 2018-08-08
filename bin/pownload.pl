#!/usr/bin/perl --
#
# pownload.pl   -- a tool for downloading video from popular videoservices
#
use strict;
use warnings;
use LWP::UserAgent;
 
$| = 1;
 
my %services = (
   'http://cs.videosaver.ru/xurl/'
                => 'href="([^"]*)"\s+title="FLV',
   'http://keepvid.com/'
                => '<a href="([^"]*)"[^>]*>[^:]*\.flv',
   'http://vidirect.ru/backend/main_backend.php'
                => '^(.*)$',
   'http://0download.ru/'
                => q{copytoclipboard\('([^']*)'\)},
);
 
die "Usage: pownload.pl VideoURL [outputfile]\n" unless @ARGV;
 
my ($videourl, $outflv) = @ARGV;
 
my $i = 0;
$outflv = sprintf('video_%04d.flv', $i++)
   while !$outflv or -e $outflv;
 
my $ua = LWP::UserAgent->new;
 
foreach my $service (keys %services) {
    print "Trying $service...\n";
    my $response = $ua->get($service . '?url=' . $videourl);
    next unless $response->is_success
                && $response->content =~ /$services{$service}/;
 
    print "Saving video to $outflv...\n";
 
    my $received_size = 0;
    open my $fh, '>' , $outflv or die $!;
    binmode $fh;
    $ua->get($1, ':content_cb' => sub {
       my ($data, $response) = @_;
       $received_size += length $data;
       print {$fh} $data;
       printf "\r" . [qw(- \\ | /)]->[++$i % 4] . ' %d%%         ',
         100 * $received_size / ($response->header('Content-Length') || 1);
     });
    print "\nDone.\n";
    exit;
}
 
print "Sorry.\n";
