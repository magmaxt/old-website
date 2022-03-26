# $id$
# Get a picture from the wiki and compress it

use strict;
use LWP::UserAgent;
require HTTP::Request;

die "Usage \"perl getpic.pl <name>\"" if ($#ARGV < 0);

my $name = $ARGV[0];
my $ua = new LWP::UserAgent;


my $request = HTTP::Request->new(GET => "https://sourceforge.net/apps/mediawiki/dealii/index.php?title=Image:$name");

my $r = $ua->request($request);

die "HTTP request failed: ", $r->message unless ($r->is_success);
die "No HTML: ", $r->content_type unless ($r->content_type =~ /^text\/html/);

${$r->content_ref} =~ m/=\"(\/apps/mediawiki[^\":]*$name)\">/;
my $path = $1;

print $path, "\n";

my $imrequest = new HTTP::Request(GET => "http://sourceforge.net$path");
$r = $ua->request($imrequest);

die "Not an image: ", $r->content_type unless ($r->content_type =~ /^image/);

print "Content type: ", $r->content_type, "\n";

open FH, ">raw$name";
syswrite FH, ${$r->content_ref};
close FH;

system ('convert', '-size', '100x150', '-resize', '100x150', "raw$name", $name);
