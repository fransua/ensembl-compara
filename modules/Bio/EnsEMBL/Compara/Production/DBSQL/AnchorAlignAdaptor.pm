=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 NAME

Bio::EnsEMBL::Compara::Production::DBSQL::AnchorAlignAdaptor

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 APPENDIX

=cut

use strict;
use warnings;


package Bio::EnsEMBL::Compara::Production::DBSQL::AnchorAlignAdaptor;

use Data::Dumper;

use Bio::EnsEMBL::Compara::Production::EPOanchors::AnchorAlign;

use Bio::EnsEMBL::Utils::Exception qw(throw);

use base qw(Bio::EnsEMBL::Compara::DBSQL::BaseAdaptor);


#############################
#
# store methods
#
#############################

=head2 store

  Arg[1]     : one or many DnaFragChunk objects
  Example    : $adaptor->store($chunk);
  Description: stores DnaFragChunk objects into compara database
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub store {
  my ($self, $anchor_align)  = @_;

  throw() unless($anchor_align);
  throw() unless(UNIVERSAL::isa($anchor_align, 'Bio::EnsEMBL::Compara::Production::EPOanchors::AnchorAlign'));

  my $query = qq{
  INSERT INTO anchor_align
    (method_link_species_set_id, anchor_id, dnafrag_id, dnafrag_start,
    dnafrag_end, dnafrag_strand, score, num_of_organisms, num_of_sequences, evalue)
  VALUES (?,?,?,?,?,?,?,?,?,?)};

  my $sth = $self->prepare($query);
  my $insertCount =
     $sth->execute($anchor_align->method_link_species_set_id,
        $anchor_align->anchor_id,
        $anchor_align->dnafrag_id,
        $anchor_align->dnafrag_start,
        $anchor_align->dnafrag_end,
        $anchor_align->dnafrag_strand,
        $anchor_align->score,
        $anchor_align->num_of_organisms,
        $anchor_align->num_of_sequences,
	$anchor_align->evalue,
        );
  if($insertCount>0) {
    #sucessful insert
    $anchor_align->dbID( $self->dbc->db_handle->last_insert_id(undef, undef, 'anchor_align', 'anchor_align_id') );
    $sth->finish;
  }

  $anchor_align->adaptor($self);

  return $anchor_align;
}

sub store_mapping_hits {
	my $self = shift;
	my $batch_records = shift;
	my $out_put_mlssid = shift;
	throw() unless($batch_records);
	
        # FIXME: disconnect_when_inactive(): why do we need a LOCK here ?
	my $dcs = $self->dbc->disconnect_when_inactive();
	$self->dbc->disconnect_when_inactive(0);
	$self->dbc->do("LOCK TABLE anchor_align WRITE");

	my $query = qq{
	INSERT INTO anchor_align (method_link_species_set_id, anchor_id, dnafrag_id, dnafrag_start,	
	dnafrag_end, dnafrag_strand, score, num_of_organisms, num_of_sequences, evalue, anchor_status)
	VALUES (?,?,?,?,?,?,?,?,?,?,?)};

	my $sth = $self->prepare($query);
	foreach my $anchor_hits( @$batch_records ) {
		$sth->execute( @{ $anchor_hits } );
	}	
	$sth->finish;
	$self->dbc->do("UNLOCK TABLES");
	$self->dbc->disconnect_when_inactive($dcs);
	return 1;
}

sub store_exonerate_hits {
        my $self = shift;
        my $batch_records = shift;
        my $out_put_mlssid = shift;
        throw() unless($batch_records);
    
        my $query = qq{ 
        INSERT INTO anchor_align (method_link_species_set_id, anchor_id, dnafrag_id, dnafrag_start,     
        dnafrag_end, dnafrag_strand, score, num_of_organisms, num_of_sequences)
        VALUES (?,?,?,?,?,?,?,?,?)};

        my $sth = $self->prepare($query);
        foreach my $row(@$batch_records) {
                $sth->execute( split(":", $row) );
        }    
        $sth->finish;
        return 1;
}


###############################################################################
#
# fetch methods
#
###############################################################################


sub fetch_dnafrag_id {
	my $self = shift;
	my($coord_sys, $dnafrag_name, $target_genome_db_id) = @_;
	unless (defined($coord_sys) and defined($dnafrag_name) and defined($target_genome_db_id)) {
		throw("fetch_dnafrag_id must have a coord_sys, dnafrag_name and target_genome_db_id");
	}
	my $query = qq{
		SELECT dnafrag_id FROM dnafrag WHERE name = ? AND 
		coord_system_name = ? AND genome_db_id = ?};
	my $sth = $self->prepare($query);
	$sth->execute($dnafrag_name, $coord_sys, $target_genome_db_id) or die $self->errstr;
	while (my$row = $sth->fetchrow_arrayref) {
		return $row->[0];
	}
}


##########################


=head2 fetch_all_by_anchor_id_and_mlss_id

  Arg[1]     : anchor_id, string
  Arg[2]     : method_link_species_set_id, string
  Example    : my $anchor = $anchor_align_adaptor->fetch_all_by_anchor_id_and_mlss_id($self->input_anchor_id,$self->method_link_species_set_id);
  Description: returns hashref of cols. from anchor_align table using anchor_align_id as unique hash key
  Returntype : hashref 
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_by_anchor_id_and_mlss_id {
	my ($self, $anchor_id, $method_link_species_set_id) = @_;
	unless (defined $anchor_id && defined $method_link_species_set_id) {
		throw("fetch_all_by_anchor_id_and_mlss_id must have an anchor_id and a method_link_species_set_id");
	}

	my $query = qq{
		SELECT anchor_align_id, method_link_species_set_id, anchor_id, 
		dnafrag_id, dnafrag_start, dnafrag_end, dnafrag_strand, score, 
		num_of_organisms, num_of_sequences FROM anchor_align WHERE 
		anchor_id = ? AND method_link_species_set_id = ? AND anchor_status IS NULL};
	my $sth = $self->prepare($query);
	$sth->execute($anchor_id, $method_link_species_set_id) or die $self->errstr;
	return $sth->fetchall_hashref("anchor_align_id");
}


=head2 update_anchor_status

  Arg[1]     : anchor_id, hashref 
  Arg[2]     : integer: new "anchor_status" value
  Example    : $anchor_align_adaptor->update_anchor_status($self->input_anchor_id, 3333);
  Description: updates anchor_status field
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub update_anchor_status {
	my($self, $failed_anchor_hash_ref, $new_status, $mlssid) = @_;
	unless (defined $failed_anchor_hash_ref ){
		throw( "No failed_anchor_id : update_anchor_status failed");
	} 
	unless (defined $mlssid) {
		throw("No mlssid : update_anchor_status failed");
	}

	my $update = qq{
		UPDATE anchor_align SET anchor_status = ? WHERE anchor_id = ? AND method_link_species_set_id = ?};
	my $sth = $self->prepare($update);
	foreach my $failed_anchor(%{$failed_anchor_hash_ref}) {
		$sth->execute($new_status, $failed_anchor, $mlssid) or die $self->errstr;
	}
	return 1;
}


# FIXME: should in theory return the genome_dbs of mlss' species_set
sub fetch_all_genome_db_ids_for_mlssid {
	my($self, $mlssid) = @_;
	my $dnafrag_query = qq{
		SELECT DISTINCT(df.genome_db_id) FROM anchor_align aa
		INNER JOIN dnafrag df on aa.dnafrag_id = df.dnafrag_id 
		WHERE aa.method_link_species_set_id = ?};
        return $self->dbc->db_handle->selectcol_arrayref($dnafrag_query, undef, $mlssid);
}

=head2 fetch_all_anchors_by_genome_db_id_and_mlssid 

  Arg[0]     : genome_db_id, string
  Arg[1]     : mlssid, string
  Example    : 
  Description: 
  Returntype : arrayref 
  Exceptions : none
  Caller     : general

=cut

#HACK

sub fetch_all_anchors_by_genome_db_id_and_mlssid {
	my($self, $genome_db_id, $mlssid) = @_;
	unless (defined $genome_db_id && defined $mlssid) {
		throw("fetch_all_anchors_by_genome_db_id_and_mlssid must have a genome_db_id and a mlssid");
	}
	my $dnafrag_query = qq{
		SELECT aa.dnafrag_id, aa.anchor_align_id, aa.anchor_id, aa.dnafrag_start, aa.dnafrag_end 
		FROM anchor_align aa
		INNER JOIN dnafrag df ON df.dnafrag_id = aa.dnafrag_id 
		WHERE df.genome_db_id = ? AND aa.method_link_species_set_id = ? AND anchor_status 
		IS NULL ORDER BY dnafrag_start, dnafrag_end};
	my $sth = $self->prepare($dnafrag_query);
	$sth->execute($genome_db_id, $mlssid) or die $self->errstr;
	return $sth->fetchall_arrayref();
}


############################
#
# INTERNAL METHODS
#
############################

#internal method used in multiple calls above to build objects from table data

sub _tables {
  my $self = shift;

  return (['anchor_align', 'aa'] );
}

sub _columns {
  my $self = shift;

  return qw (aa.anchor_align_id
             aa.method_link_species_set_id
             aa.anchor_id
             aa.dnafrag_id
             aa.dnafrag_start
             aa.dnafrag_end
             aa.dnafrag_strand
             aa.score
             aa.num_of_organisms
             aa.num_of_sequences
             aa.anchor_status
            );
}



sub _objs_from_sth {
  my ($self, $sth) = @_;

  my @anchor_aligns = ();

  while( my $row_hashref = $sth->fetchrow_hashref()) {
    my $this_anchor_align = Bio::EnsEMBL::Compara::Production::EPOanchors::AnchorAlign->new();

    $this_anchor_align->adaptor($self);
    $this_anchor_align->dbID($row_hashref->{'anchor_align_id'});
    $this_anchor_align->method_link_species_set_id($row_hashref->{'method_link_species_set_id'});
    $this_anchor_align->anchor_id($row_hashref->{'anchor_id'});
    $this_anchor_align->dnafrag_id($row_hashref->{'dnafrag_id'});
    $this_anchor_align->dnafrag_start($row_hashref->{'dnafrag_start'});
    $this_anchor_align->dnafrag_end($row_hashref->{'dnafrag_end'});
    $this_anchor_align->dnafrag_strand($row_hashref->{'dnafrag_strand'});
    $this_anchor_align->score($row_hashref->{'score'});
    $this_anchor_align->num_of_organisms($row_hashref->{'num_of_organisms'});
    $this_anchor_align->num_of_sequences($row_hashref->{'num_of_sequences'});
    $this_anchor_align->anchor_status($row_hashref->{'anchor_status'});

    push @anchor_aligns, $this_anchor_align;
  }
  $sth->finish;

  return \@anchor_aligns;
}


1;
