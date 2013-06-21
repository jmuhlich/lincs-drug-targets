#!/usr/bin/perl -w
use strict;
use LWP::Simple;

my @cols = qw(
               gd_hgnc_id
               gd_app_sym
               gd_app_name
               gd_prev_sym
               gd_aliases
               gd_name_aliases
               gd_pub_acc_ids
               gd_gene_fam_name
               md_prot_id
             );

my %params = (
               status => 'Approved',
               status_opt => '2',
               where => '',
               order_by => 'gd_hgnc_id',
               format => 'text',
               limit => '',
               hgnc_dbtag => 'on',
               submit => 'submit'
             );

my $url_base = 'http://www.genenames.org/cgi-bin/hgnc_downloads?';

my $cols = join('&', map "col=$_", @cols);
my $params = join('&', map "$_=$params{$_}", keys %params);
my $url = "$url_base$cols&$params\n";

print get( $url );
