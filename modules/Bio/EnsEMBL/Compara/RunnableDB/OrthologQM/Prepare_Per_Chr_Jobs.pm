=pod
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

        Bio::EnsEMBL::Compara::RunnableDB::OrthologQM::Prepare_Per_Chr_Jobs

=head1 SYNOPSIS

=head1 DESCRIPTION
  Takes as input an hash of reference and non reference species DBIDs and  dnafrag DBIDs as keys and the list of homolog DBIDs as values.
  Unpacks the hash into seperate hashes each containing a single dnafrag DBID as the key to a list of ordered homologs.
  if the goc_reuse_db parameter is set,
   theis runnable will reuse goc scores from the previous db and only calculate goc scores for new homologs and the 2 homologies closest to them on each side.

    Example run

  standaloneJob.pl Bio::EnsEMBL::Compara::RunnableDB::OrthologQM::Prepare_Per_Chr_Jobs

=cut

package Bio::EnsEMBL::Compara::RunnableDB::OrthologQM::Prepare_Per_Chr_Jobs;

use strict;
use warnings;
use Data::Dumper;
use base 'Bio::EnsEMBL::Compara::RunnableDB::OrthologQM::Compare_orthologs';

use Bio::EnsEMBL::Registry;

sub param_defaults {
  my $self = shift;
  return {
    %{ $self->SUPER::param_defaults() },
#      'ortholog_info_hashref' => {'13200473' => {
#                          '153474' => '77514854',
#                          '70578' => '69816616',
#                          '225827' => '102025881',
#                          '153125' => '89554512',
#                          '171749' => '82542192',
#                          '40769' => '115084016'
#                        },

#    },
#    'ref_species_dbid'    => 122,
#    'non_ref_species_dbid'    => 144,
#    'goc_reuse_db'  =>  'mysql://ensro@compara2/mp14_protein_trees_85',
#    'compara_db'  => 'mysql://ensadmin:'.$ENV{ENSADMIN_PSW}.'@compara5/cc21_protein_trees_no_reuse_86',
#    'goc_mlss_id'  => '50062',
    };
}


=head2 fetch_input

    Description : Implements fetch_input() interface method of Bio::EnsEMBL::Hive::Process that is used to read in parameters and load data.
    Here I need to retrieve the ortholog hash that was data flowed here by Prepare_Othologs.pm

=cut

sub fetch_input {
  my $self = shift;
  $self->debug(4);
  print "Bio::EnsEMBL::Compara::RunnableDB::OrthologQM::Prepare_Per_Chr_Jobs --------------------------------START  \n  " if ( $self->debug );
  my $mlss_id = $self->param_required('goc_mlss_id');
  $self->param('homolog_adaptor', $self->compara_dba->get_HomologyAdaptor);
  $self->param('mlss_adaptor', $self->compara_dba->get_MethodLinkSpeciesSetAdaptor);

  my $ref_species_dbid = $self->param('ref_species_dbid');
  $self->param('mlss_check', 0); # will be use to check if the mlss exist the reuse db ortherwise the goc for this mlss will be recalculated
  #this variables name was change to 'goc_reuse_db' as it was getting mixed up with another variable. The old name was left here to make it compatible with older runs of the runnable.

  if ($self->param('goc_reuse_db') ) {
    my $mlss = $self->param('mlss_adaptor')->fetch_by_dbID($mlss_id);
    my $prev_release_mlss_id = $mlss->get_tagvalue('prev_release_mlss_id');
    print "\n\n  curr_release_mlss_id   :  $mlss_id  , and , prev_release_mlss_id :  $prev_release_mlss_id\n\n" if ( $self->debug >3 );

    if (defined $prev_release_mlss_id) {
      $self->param('mlss_check', 1);

  #preload all gene members to make quering the homologs faster later on
      my $homologs = $self->param('homolog_adaptor')->fetch_all_by_MethodLinkSpeciesSet($mlss);
      my $sms = Bio::EnsEMBL::Compara::Utils::Preloader::expand_Homologies($homologs->[0]->adaptor->db->get_AlignedMemberAdaptor, $homologs);
      Bio::EnsEMBL::Compara::Utils::Preloader::load_all_GeneMembers($homologs->[0]->adaptor->db->get_GeneMemberAdaptor, $sms);
      my $preloaded_homologs_hashref;

  #loop through the preloaded homologies to create a hash table homology id => homology object. this will serve as a look up table 
      while ( my $ortholog = shift( @{ $homologs} ) ) {
        $preloaded_homologs_hashref->{$ortholog->dbID()} = $ortholog;
      }
      $self->param('preloaded_homologs', $preloaded_homologs_hashref);

  # now we create a hash object of homology id mapping for this mlss id 
      my $query = "SELECT curr_release_homology_id, prev_release_homology_id FROM homology_id_mapping where mlss_id = $mlss_id";
      my $hID_map = $self->compara_dba->dbc->db_handle->selectall_arrayref($query);
      my %homologyID_map = map { $_->[0] => $_->[1] } @{$hID_map} ;
      $self->param('homologyID_map', \%homologyID_map);
    
    #now we will query the ortholog_goc_metric table uploaded from from the previous db using the prev mlss id that maps to the new mlss id
      my $sql = "select homology_id, goc_score, left1, left2, right1, right2 from prev_rel_goc_metric ogm
              join prev_rel_gene_member gm using (gene_member_id) where gm.genome_db_id = $ref_species_dbid and method_link_species_set_id = $prev_release_mlss_id";
      my $prev_goc_hashref = $self->compara_dba->dbc->db_handle->selectall_hashref($sql, "homology_id");
      print Dumper($prev_goc_hashref) if ( $self->debug >5 );
      $self->param('prev_goc_hashref', $prev_goc_hashref);
    }
  }
}

sub run {
  my $self = shift;
  print "the runnable Prepare_Per_Chr_Jobs ----------start \n mlss_id ---->   ", $self->param('goc_mlss_id')   if ( $self->debug >3);
  if ($self->param('mlss_check')) {
    $self->_reusable_species();
    $self->_insert_goc_scores();
  }
  else {
    $self->_non_reusable_species();
  }

  print "the runnable Prepare_Per_Chr_Jobs ----------END \n mlss_id \n", $self->param('goc_mlss_id') , if ( $self->debug >3);
}

sub _reusable_species {
  my $self = shift;
  print "\n Starting ------  _reusable_species \n" if ( $self->debug );
  my @goc_score_array;
  $self->param('goc_score_arrayref', \@goc_score_array);
  my $ortholog_hashref = $self->param_required('ortholog_info_hashref');
  my $count_homologs = 0;
  my $count_recal_homologs = 0;
  my $count_new_homologs =0;

  while (my ($ref_dnafragID, $chr_orth_hashref) = each(%$ortholog_hashref) ) {
    my $orth_sorted = $self->_order_chr_homologs($chr_orth_hashref); # will contain the orthologs ordered by the dnafrag start position #sorting the orthologs by dnafrag start position
    for (my $homolog_index = 0; $homolog_index < scalar @{$orth_sorted}; $homolog_index++ ) {
      $count_homologs ++;
      my $curr_homology_id = $orth_sorted->[$homolog_index];
      my $prev_homology_id = $self->param('homologyID_map')->{$curr_homology_id};


      if ( defined $prev_homology_id) {
        my $curr_homology_obj = $self->param('preloaded_homologs')->{$curr_homology_id};
        my $curr_ref_gmem_dbID = $curr_homology_obj->get_all_GeneMembers($self->param('ref_species_dbid'))->[0]->dbID;

        #now we will query the prev_goc_hashref gotten from the previous db using the prev homology id that maps to the new homology id
        my $homology_goc_score = $self->param('prev_goc_hashref')->{$prev_homology_id};

        if (! defined $homology_goc_score->{'goc_score'} ) {
          $count_new_homologs +=1;
          #this only happens if an homolog is a 1 to many ortholog and that 1 also happens to be the only gene on the chr. there will still be goc scores for the homologmany genes on the
          if ($homolog_index == 0 ) {
            $self->_new_homolog_at_0($orth_sorted,$homolog_index,$ref_dnafragID);
            $homolog_index +=2;
            $count_homologs +=2;
            $count_recal_homologs +=3;
          }
          elsif ($homolog_index == 1) {
            $self->_new_homolog_at_1($orth_sorted,$homolog_index,$ref_dnafragID);
            $homolog_index +=2;
            $count_homologs +=2;
            $count_recal_homologs +=4;
          }
          else {
            $self->_new_homolog($orth_sorted,$homolog_index,$ref_dnafragID);
            $homolog_index +=2;
            $count_homologs +=2;
            $count_recal_homologs +=5;
          }
            next;
        }

        $homology_goc_score->{'gene_member_id'} = $curr_ref_gmem_dbID;
        $homology_goc_score->{'homology_id'} = $curr_homology_id;
        $homology_goc_score->{'dnafrag_id'} = $ref_dnafragID;
        $homology_goc_score->{'method_link_species_set_id'} = $self->param_required('goc_mlss_id');
        push (@{$self->param('goc_score_arrayref')}, $homology_goc_score) ;
      }
      else {
        $count_new_homologs +=1;
        if ($homolog_index == 0 ) {
          $self->_new_homolog_at_0($orth_sorted,$homolog_index,$ref_dnafragID);
          $homolog_index +=2;
          $count_homologs +=2;
          $count_recal_homologs +=3;
        }
        elsif ($homolog_index == 1) {
          $self->_new_homolog_at_1($orth_sorted,$homolog_index,$ref_dnafragID);
          $homolog_index +=2;
          $count_homologs +=2;
          $count_recal_homologs +=4;
        }
        else {
          $self->_new_homolog($orth_sorted,$homolog_index,$ref_dnafragID);
          $homolog_index +=2;
          $count_homologs +=2;
          $count_recal_homologs +=5;
        }
      }
    }
  }

  print "\n\n This is how many homology ids were in the input : $count_homologs  \n this is how many were recalculated  : $count_recal_homologs  \n this is the actual count of new homologs :  $count_new_homologs \n" if ( $self->debug );

}

sub _new_homolog_at_0 {
  my $self = shift;
  print "\n Starting ------  _new_homolog_at_0 \n" if ( $self->debug );
  my ($sorted_homologs,$h_index,$curr_ref_dnafragID) = @_;
  my ($left1, $left2, $right1, $right2, $query) = (undef,undef,undef,undef);
      #recalculate the goc score for the query (new) homolog
  ($right1, $right2) = ($sorted_homologs->[$h_index +1], $sorted_homologs->[$h_index +2]);
  $query = $sorted_homologs->[$h_index];
  my $goc_score_hashref = $self->_compute_ortholog_score($left1, $left2, $query, $right1, $right2, $curr_ref_dnafragID );
  push @{$self->param('goc_score_arrayref')}, $goc_score_hashref ;

  #recalculate the goc score for the homolog at right1 of the query homolog (query position +1 is now the query)
  if (defined $sorted_homologs->[$h_index +1]) {
  ($left2, $left1, $query, $right1, $right2) = ($sorted_homologs->[$h_index -1], $sorted_homologs->[$h_index], $sorted_homologs->[$h_index +1], $sorted_homologs->[$h_index +2], $sorted_homologs->[$h_index +3]);
    $goc_score_hashref = $self->_compute_ortholog_score($left1, $left2, $query, $right1, $right2, $curr_ref_dnafragID );
    push @{$self->param('goc_score_arrayref')}, $goc_score_hashref ;

    if (defined $sorted_homologs->[$h_index +2]) {
    #recalculate the goc score for the homolog at right2 of the query homolog (query position +2 is now the query)
      ($left2, $left1, $query, $right1, $right2) = ($sorted_homologs->[$h_index], $sorted_homologs->[$h_index +1], $sorted_homologs->[$h_index +2], $sorted_homologs->[$h_index +3], $sorted_homologs->[$h_index +4]);
      $goc_score_hashref = $self->_compute_ortholog_score($left1, $left2, $query, $right1, $right2, $curr_ref_dnafragID );
      push @{$self->param('goc_score_arrayref')}, $goc_score_hashref ;
    }
  }
}

sub _new_homolog_at_1 {
  my $self = shift;
  print "\n Starting ------  _new_homolog_at_1 \n" if ( $self->debug );
  my ($sorted_homologs,$h_index,$curr_ref_dnafragID) = @_;
  my ($left1, $left2, $right1, $right2, $query) = (undef,undef,undef,undef);
  $self->_delete_from_goc_score_array();
  #recalculate the goc score for the homolog at left1 of the query homolog (query position -1 is now the query)
  ($query, $right1, $right2) = ($sorted_homologs->[$h_index -1], $sorted_homologs->[$h_index], $sorted_homologs->[$h_index +1]);
  my $goc_score_hashref = $self->_compute_ortholog_score($left1, $left2, $query, $right1, $right2, $curr_ref_dnafragID );
  push @{$self->param('goc_score_arrayref')}, $goc_score_hashref ;
  #recalculate the goc score for the query (new) homolog
  ($left1, $query, $right1, $right2) = ($sorted_homologs->[$h_index -1],$sorted_homologs->[$h_index],$sorted_homologs->[$h_index +1],$sorted_homologs->[$h_index +2]);
  $goc_score_hashref = $self->_compute_ortholog_score($left1, $left2, $query, $right1, $right2, $curr_ref_dnafragID );
  push @{$self->param('goc_score_arrayref')}, $goc_score_hashref ;
  #recalculate the goc score for the homolog at right1 of the query homolog (query position +1 is now the query)
  if (defined $sorted_homologs->[$h_index +1]) {
  ($left2, $left1, $query, $right1, $right2) = ($sorted_homologs->[$h_index -1], $sorted_homologs->[$h_index], $sorted_homologs->[$h_index +1], $sorted_homologs->[$h_index +2], $sorted_homologs->[$h_index +3]);
    $goc_score_hashref = $self->_compute_ortholog_score($left1, $left2, $query, $right1, $right2, $curr_ref_dnafragID );
    push @{$self->param('goc_score_arrayref')}, $goc_score_hashref ;

    if (defined $sorted_homologs->[$h_index +2]) {
    #recalculate the goc score for the homolog at right2 of the query homolog (query position +2 is now the query)
      ($left2, $left1, $query, $right1, $right2) = ($sorted_homologs->[$h_index], $sorted_homologs->[$h_index +1], $sorted_homologs->[$h_index +2], $sorted_homologs->[$h_index +3], $sorted_homologs->[$h_index +4]);
      $goc_score_hashref = $self->_compute_ortholog_score($left1, $left2, $query, $right1, $right2, $curr_ref_dnafragID );
      push @{$self->param('goc_score_arrayref')}, $goc_score_hashref ;
    }
  }
}

sub _new_homolog {
  my $self = shift;
  print "\n Starting ------  _new_homolog \n" if ( $self->debug );
  my ($sorted_homologs,$h_index,$curr_ref_dnafragID) = @_;
  my ($left1, $left2, $right1, $right2, $query) = (undef,undef,undef,undef);
  $self->_delete_from_goc_score_array();
  $self->_delete_from_goc_score_array();
  #recalculate the goc score for the homolog at left2 of the query homolog (query position -2 is now the query)
  ($left2, $left1, $query, $right1, $right2) = ($sorted_homologs->[$h_index -4], $sorted_homologs->[$h_index -3], $sorted_homologs->[$h_index -2], $sorted_homologs->[$h_index -1], $sorted_homologs->[$h_index]);
  my $goc_score_hashref = $self->_compute_ortholog_score($left1, $left2, $query, $right1, $right2, $curr_ref_dnafragID );
  push @{$self->param('goc_score_arrayref')}, $goc_score_hashref ;
  #recalculate the goc score for the homolog at left1 of the query homolog (query position -1 is now the query)
  ($left2, $left1, $query, $right1, $right2) = ($sorted_homologs->[$h_index -3], $sorted_homologs->[$h_index -2], $sorted_homologs->[$h_index -1], $sorted_homologs->[$h_index], $sorted_homologs->[$h_index +1]);
  $goc_score_hashref = $self->_compute_ortholog_score($left1, $left2, $query, $right1, $right2, $curr_ref_dnafragID );
  push @{$self->param('goc_score_arrayref')}, $goc_score_hashref ;
  #recalculate the goc score for the query (new) homolog
  ($left2, $left1, $query, $right1, $right2) = ($sorted_homologs->[$h_index -2], $sorted_homologs->[$h_index -1], $sorted_homologs->[$h_index], $sorted_homologs->[$h_index +1], $sorted_homologs->[$h_index +2]);
  $goc_score_hashref = $self->_compute_ortholog_score($left1, $left2, $query, $right1, $right2, $curr_ref_dnafragID );
  push @{$self->param('goc_score_arrayref')}, $goc_score_hashref ;
  #recalculate the goc score for the homolog at right1 of the query homolog (query position +1 is now the query)
  if (defined $sorted_homologs->[$h_index +1]) {
  ($left2, $left1, $query, $right1, $right2) = ($sorted_homologs->[$h_index -1], $sorted_homologs->[$h_index], $sorted_homologs->[$h_index +1], $sorted_homologs->[$h_index +2], $sorted_homologs->[$h_index +3]);
    $goc_score_hashref = $self->_compute_ortholog_score($left1, $left2, $query, $right1, $right2, $curr_ref_dnafragID );
    push @{$self->param('goc_score_arrayref')}, $goc_score_hashref ;

    if (defined $sorted_homologs->[$h_index +2]) {
    #recalculate the goc score for the homolog at right2 of the query homolog (query position +2 is now the query)
      ($left2, $left1, $query, $right1, $right2) = ($sorted_homologs->[$h_index], $sorted_homologs->[$h_index +1], $sorted_homologs->[$h_index +2], $sorted_homologs->[$h_index +3], $sorted_homologs->[$h_index +4]);
      $goc_score_hashref = $self->_compute_ortholog_score($left1, $left2, $query, $right1, $right2, $curr_ref_dnafragID );
      push (@{$self->param('goc_score_arrayref')}, $goc_score_hashref ) ;
    }
  }
}

sub _delete_from_goc_score_array {
  my $self = shift;
  print "\n Starting ------  _delete_from_goc_score_array \n" if ( $self->debug );
  my $trash = pop @{$self->param('goc_score_arrayref')};
  print "this score has been infected so must recalculated " if ( $self->debug >2);
  print Dumper($trash) if ( $self->debug >3);

}

#this method will create one sql insert statement out of the goc_scores array and insert it into the ortholog_goc_metric table
sub _insert_goc_scores {
  my $self = shift;
  print "\n Starting ------  _insert_goc_scores \n" if ( $self->debug );

  my $sql_insert_query = 'INSERT INTO ortholog_goc_metric (method_link_species_set_id,homology_id,gene_member_id,dnafrag_id,goc_score,left1,left2,right1,right2) VALUES ' ;
  my $size_goc_score_arrayref = scalar @{$self->param('goc_score_arrayref')};
  print "\n\n This is the how many homology ids we have goc scores for and inserting into the ortholog goc table :  \n  $size_goc_score_arrayref \n\n " if ( $self->debug );
  foreach ( @{$self->param('goc_score_arrayref')} ) {
    my $l1 = defined $_->{'left1'} ? $_->{'left1'} : 'NULL';
    my $l2 = defined $_->{'left2'} ? $_->{'left2'} : 'NULL';
    my $r1 = defined $_->{'right1'} ? $_->{'right1'} : 'NULL';
    my $r2 = defined $_->{'right2'} ? $_->{'right2'} : 'NULL';
    defined $_->{'method_link_species_set_id'} or die " the method_link_species_set_id can not be empty. homology_id : $_->{'homology_id'} ";
    defined $_->{'goc_score'} or die " the goc score should atleast be 0. mlss id :  $_->{'method_link_species_set_id'} , homology_id : $_->{'homology_id'}, goc score : $_->{'goc_score'} , l1 : $_->{'left1'} , l2 : $_->{'left2'} , r1 : $_->{'right1'}, r2 : $_->{'right2'}  ";
    my $tmp = ' (' . $_->{'method_link_species_set_id'} . ',' . $_->{'homology_id'} . ',' . $_->{'gene_member_id'} . ',' . $_->{'dnafrag_id'} .
              ',' . $_->{'goc_score'} . ',' . $l1 . ',' . $l2 . ',' . $r1 . ',' . $r2 . '),' ;

    $sql_insert_query .= $tmp ;
  }

  my $insert_query = substr($sql_insert_query, 0, -1);
  my $goc_insert =  $self->compara_dba->dbc->db_handle->prepare($insert_query) or die "Couldn't prepare query: ".$self->compara_dba->dbc->db_handle->errstr;

  $goc_insert->execute() or die "Couldn't execute query: ".$self->compara_dba->dbc->db_handle->errstr;
}

#this method creates a per chromosome array of ordered homologies which it dataflows to the compare_orthologs.pm
sub _non_reusable_species {
  my $self = shift;
  print "\n Starting ------  _non_reusable_species \n" if ( $self->debug );
  my $ortholog_hashref = $self->param_required('ortholog_info_hashref');
  while (my ($ref_dnafragID, $chr_orth_hashref) = each(%$ortholog_hashref) ) {
    my $orth_sorted = $self->_order_chr_homologs($chr_orth_hashref); # will contain the orthologs ordered by the dnafrag start position #sorting the orthologs by dnafrag start position
    my $chr_job = {};
    $chr_job->{$ref_dnafragID} = $orth_sorted;
    $self->dataflow_output_id( {'chr_job' => $chr_job, 'ref_species_dbid' => $self->param('ref_species_dbid'), 'non_ref_species_dbid' => $self->param('non_ref_species_dbid') }, 2 );
  }
}

#this method sorts an hash of homology id and dnafrag_starts into an ordered array of homology_ids based on the dnafrag starts
sub _order_chr_homologs {
  my $self = shift;
  print "\n Starting ------  _order_chr_homologs \n" if ( $self->debug );
  my ($unsorted_chr_orth_hashref) = @_;
  my @sorted_orth;

  foreach my $name (sort { ($unsorted_chr_orth_hashref->{$a} <=> $unsorted_chr_orth_hashref->{$b}) || ($a <=> $b) } keys %$unsorted_chr_orth_hashref ) {
    printf "%-8s %s \n", $name, $unsorted_chr_orth_hashref->{$name} if ( $self->debug >3);
    push @sorted_orth, $name;

  }

  return \@sorted_orth;
}

1;
