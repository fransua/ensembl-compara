use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use Bio::EnsEMBL::Compara::Production::Projection::RunnableDB::ProjectOntologyXref;
use Bio::EnsEMBL::Compara::Production::Projection::RunnableDB::RunnableLogger;  
use Bio::EnsEMBL::Registry;

my $log_config = <<LOGCFG;
log4perl.logger=DEBUG, Screen
log4perl.appender.Screen=Log::Log4perl::Appender::Screen
log4perl.appender.Screen.stderr=0
log4perl.appender.Screen.Threshold=DEBUG
log4perl.appender.Screen.layout=Log::Log4perl::Layout::PatternLayout
log4perl.appender.Screen.layout.ConversionPattern=%d %p> %F{1}:%L %M - %m%n
LOGCFG

my @options = qw( 
  source=s 
  target=s 
  engine=s 
  compara=s 
  write_to_db 
  file=s 
  registry=s 
  verbose help man 
);

#The only thing we run in global
run();
#End of the script

sub run {
  my $opts = _get_opts();
  _initalise_log();
  my $runnable = _build_runnable($opts);
  $runnable->run_without_hive();
  return;
}

sub _get_opts {
  my $opts = {};
  GetOptions($opts, @options) or pod2usage(1);
  pod2usage( -exitstatus => 0, -verbose => 1 ) if $opts->{help};
	pod2usage( -exitstatus => 0, -verbose => 2 ) if $opts->{man};
	
	#Source & target check
	_exit('No -source option given', 1, 1) if ! $opts->{source};
	_exit('No -target option given', 1, 1) if ! $opts->{target};
	_exit('No -compara option given', 1, 1) if ! $opts->{compara};
	
	#Registry work
	my $reg = $opts->{registry};
	_exit('No -registry option given', 2, 1) if ! $reg && ! -f $reg;
	my @args = ($reg);
	push @args, 1 if $opts->{verbose};
	Bio::EnsEMBL::Registry->load_all(@args);
	
	#Engine work
	$opts->{engine} = 'Bio::EnsEMBL::Compara::Production::Projection::GOAProjectionEngine' if ! $opts->{engine};
	if(! $opts->{write_to_db} && ! $opts->{file}) {
	  _exit('-write_to_db and -file were not specified. We need one', 3, 1);
	}
	
  return $opts;
}

sub _build_runnable {
  my ($opts) = @_;
  my %args = (
    -PROJECTION_ENGINE => _build_engine($opts),
    -TARGET_GENOME_DB => _get_genome_db($opts, $opts->{target}),
    -DEBUG => $opts->{verbose}
  );
  $args{-FILE} = $opts->{file} if $opts->{file};
  $args{-WRITE_DBA} = _get_adaptor($opts->{target}, 'core') if $opts->{write_to_db};
  return Bio::EnsEMBL::Compara::Production::Projection::RunnableDB::ProjectOntologyXref->new_without_hive(%args);
}

sub _build_engine {
  my ($opts) = @_;
  my $mod = $opts->{engine};
  _runtime_import($mod, 1);
  return $mod->new(
    -GENOME_DB => _get_genome_db($opts, $opts->{source}),
    -DBA => _get_adaptor($opts->{compara}, 'compara'),
    _log()
  );
}

sub _get_genome_db {
  my ($opts, $name) = @_;
  my $compara_dba = _get_adaptor($opts->{compara}, 'compara');
  my $core_dba = _get_adaptor($name, 'core');
  my $gdb_a = $compara_dba->get_GenomeDBAdaptor();
  my $gdb = $gdb_a->fetch_by_core_DBAdaptor($core_dba);
  return $gdb;
}

sub _get_adaptor {
  my ($name, $group) = @_;
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($name, $group);
  if(! defined $dba) {
    _exit("No adaptor for ${name} and ${group}. Check your registry and try again", 5, 1);
  }
  return $dba;
}

sub _exit {
  my ($msg, $status, $verbose) = @_;
  print STDERR $msg, "\n";
  pod2usage( -exitstatus => $status, -verbose => $verbose);
}

my $log4perl_available = 0;

sub _initalise_log {
  if(_runtime_import('Log::Log4perl')) {
    Log::Log4perl->init(\$log_config);
    $log4perl_available = 1;
  }
}

#If log4perl was available let the module get it's own logger otherwise we 
#build our own
sub _log {
  my ($opts) = @_;
  if($log4perl_available) {
    return;
  }
  my $log = Bio::EnsEMBL::Compara::Production::Projection::RunnableDB::RunnableLogger->new(
    -DEBUG => $opts->{verbose}
  );
  return ( -LOG => $log );
}

sub _runtime_import {
  my ($mod, $die) = @_;
  eval "require ${mod}";
  _exit "Cannot import ${mod}: $@", 5, 1 if $die && $@;
  return ($@) ? 0 : 1;
}

__END__
=pod

=head1 NAME

project_dbentry.pl

=head1 SYNOPSIS

  ./project_dbentry.pl -registry REG -source SRC -target TRG -compara COM [-engine ENG] [-write_to_db] [-file FILE] [-verbose] [-help | -man]

=head1 DESCRIPTION

This script is a thin-wrapper around the RunnableDB instance and is used
for the ad-hoc testing & running of the Xref projection engine. At the moment
this is configured for projecting GO terms from one species to another
however it will operate on any Xref so long as you can provide the correct
projection engine implementation.

The script can also add data back into a database but to do so we must
assume that a core DBAdaptor for the target species is linked to 
a read/write account. Otherwise you will not be able to perform the 
linkage.

For a flavor of what the pipeline can do pass the script a file name which
will produce a CSV of what I<would> have been written back to the DB.

=head1 OPTIONS

=over 8

=item B<--registry> - The registry to use

=item B<--source> - The source species (species with GOs)

=item B<--target> - The target species (species without GOs)

=item B<--compara> - The compara database to use

=item B<--engine> - The engine to use; defaults to GOAProjectionEngine (must be a fully qualified package)

=item B<--write_to_db> - Indicates we want Xrefs going back to the core DB. If used we assume the registry's core DBAdaptor is writable

=item B<--file> - Location to write output to. Can be a directory (so an automatically generated name will be given) or a full path

=item B<--verbose> - Start emitting more messages

=item B<--help> - Basic help with options

=item B<--man> - Manual version of the help. More complete 

=back

=head1 REQUIREMENTS

=over 8

=item EnsEMBL core (v60+)

=item EnsEMBL compara

=item EnsEMBL hive

=item Log::Log4perl - if not present on PERL5LIB or @INC messages go to STDOUT

=item Text::CSV (for file writing)

=item Data::Predicate 

=back

=cut