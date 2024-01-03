#!/usr/bin/perl

use warnings;
use strict;
use Cwd;

my $filename = $ARGV[0];

my $seq;
open( my $fh, '<', $filename ) or die "cannot open file $filename";
{
    local $/;
    $seq = <$fh>;
}
close($fh);

my $sgRNA_size_CRISPR = $ARGV[1];
my $blaster           = $ARGV[2];
my $lower_gc_CRISPR   = $ARGV[3];
my $higher_gc_CRISPR  = $ARGV[4];
my $sgRNA_PAM_CRISPR  = $ARGV[5];

my $seed_l   = $ARGV[6];
my $seed_mis = $ARGV[7];
my $tot_mis  = $ARGV[8];

my $algo = $ARGV[9];

my $id = $ARGV[10];


my %IUB = (
    'G' => 'G',
    'A' => 'A',
    'T' => 'T',
    'C' => 'C',
    'R' => '[AG]',
    'Y' => '[CT]',
    'M' => '[AC]',
    'K' => '[GT]',
    'S' => '[CG]',
    'W' => '[AT]',
    'H' => '[ACT]',
    'B' => '[CGT]',
    'V' => '[ACG]',
    'D' => '[AGT]',
    'N' => '[ACGT]'
);
my %IUB_com = (
    'G' => 'C',
    'A' => 'T',
    'T' => 'A',
    'C' => 'G',
    'R' => '[TC]',
    'Y' => '[GA]',
    'M' => '[TG]',
    'K' => '[CA]',
    'S' => '[GC]',
    'W' => '[TA]',
    'H' => '[TAG]',
    'B' => '[GCA]',
    'V' => '[TGC]',
    'D' => '[TAC]',
    'N' => '[ACGT]'
);

my %PAMHASH = (
    "T" => 0,
    "A" => 1,
    "C" => 1,
    "G" => 1
);
######################

####Global Strings
my $outfile;

#temp blast file
my $blast_rna;

#holder for sgRNA file
my @lines;

my $gRNAlength       = $sgRNA_size_CRISPR;
my $targetPAM        = $sgRNA_PAM_CRISPR;
my $targetPAM_length = length($targetPAM);
my $pam1_com         = reverse $targetPAM;
my @array_pam1       = split( //, $pam1_com );
my $pam1             = '';
for ( my $i = 0 ; $i <= $#array_pam1 ; $i++ ) {
    $pam1 = $pam1 . $IUB_com{ $array_pam1[$i] };
}
my @array_pam2 = split( //, $targetPAM );
my $pam2 = '';
for ( my $i = 0 ; $i <= $#array_pam2 ; $i++ ) {
    $pam2 = $pam2 . $IUB{ $array_pam2[$i] };
}
my $seq_len = length($seq);
my @array = split( //, $seq );

my $count = 0;

#blast file of sgRNAs
$blast_rna = "_blast_rna.temp";
open my $temp_out, ">>", $blast_rna or die;
my $output = "_sgRNA_.temp";
open my $out, ">>", $output or die;

print $out "Number\tgRNA\tstrand\tstart\tend\tGC%\tPAMPN Score\tGC%-6\tCRISPRater score\tprimers\toff-target hits\n";

for ( my $i = 0 ; $i <= $#array - $targetPAM_length + 1 ; $i++ ) {
    my $seq_pam = '';
    for ( my $num = $i ; $num <= $i + $targetPAM_length - 1 ; $num++ ) {
        $seq_pam = $seq_pam . $array[$num];
    }
    if ( $seq_pam =~ /^$pam1$/ ) {
        if ( $i + $gRNAlength + $targetPAM_length - 1 <= $#array ) {
            my $gRNA_rev = '';
            for ( my $p = $i ;
                $p < $i + $targetPAM_length + $gRNAlength ; $p++ )
            {
                $gRNA_rev = $gRNA_rev . $array[$p];
            }
            my $gRNA_com  = reverse $gRNA_rev;
            my @com_array = split( //, $gRNA_com );
            my $gRNA      = '';
            for ( my $i = 0 ; $i <= $#com_array ; $i++ ) {
                $gRNA = $gRNA . $IUB_com{ $com_array[$i] };
            }
            my $gRNA_end = $seq_len - $i - 1 + 1;
            my $gRNA_start =
              $seq_len - $i - $targetPAM_length - $gRNAlength + 1;
              my $tail   = substr $gRNA, 0, $gRNAlength;
            if (   ( gc_calc($tail) <= $higher_gc_CRISPR )
                && ( gc_calc($tail) >= $lower_gc_CRISPR ) )
            {
                $count++;
                my $PAMPN = 0;
                my $p20 = substr $tail, -1;  
                my $p16_20 = substr $tail, -4, 3;
                if (($p20 eq "G")||($p20 eq "g")){
                $PAMPN += 2;
                }
                my $PAM3 = calcPAMP($p16_20);
                my $totPAMPN = ($PAM3 + $PAMPN)/5; 
                my $gc_pc = gc_calc($tail);
                my $gRNA_tail   = substr $tail, -6;
                my $gc_pc6 = gc_calc($gRNA_tail);
                my $prim3 = primer3($seq, $gRNA_start, $gRNA_end);
                my $crisprater = getScore($gRNA);
                print $out "$count\t";
                print $out "$gRNA\t-\t$gRNA_start\t$gRNA_end\t$gc_pc\t$totPAMPN\t$gc_pc6\t$crisprater\t$prim3\n";
                print $temp_out ">_" . "$count\n$gRNA\n";
            }
        }
    }
    if ( $seq_pam =~ /^$pam2$/ ) {
        if ( $i - $gRNAlength >= 0 ) {
            my $gRNA = '';
            for (
                my $q = $i - $gRNAlength ;
                $q <= $i + $targetPAM_length - 1 ;
                $q++
              )
            {
                $gRNA = $gRNA . $array[$q];
            }
            my $gRNA_start = $i - $gRNAlength + 1;
            my $gRNA_end   = $i + $targetPAM_length - 1 + 1;
            my $tail   = substr $gRNA, 0, $gRNAlength;
            if (   ( gc_calc($gRNA) <= $higher_gc_CRISPR )
                && ( gc_calc($gRNA) >= $lower_gc_CRISPR ) )
            {
                $count++;
                my $PAMPN = 0;
                my $p20 = substr $tail, -1;  
                my $p16_20 = substr $tail, -4, 3; 
                if (($p20 eq "G")||($p20 eq "g")){
                $PAMPN += 2;
                }
                my $PAM3 = calcPAMP($p16_20);
                my $totPAMPN = ($PAM3 + $PAMPN)/5; 
                my $gc_pc = gc_calc($tail);
                my $gRNA_tail   = substr $tail, -6;
                my $gc_pc6 = gc_calc($gRNA_tail);
                my $crisprater = getScore($gRNA);
                my $prim3 = primer3($seq, $gRNA_start, $gRNA_end);
                print $out "$count\t";
                print $out "$gRNA\t+\t$gRNA_start\t$gRNA_end\t$gc_pc\t$totPAMPN\t$gc_pc6\t$crisprater\t$prim3\n";
                print $temp_out ">_" . "$count\n$gRNA\n";
            }
        }

    }
}

my $count_l = `wc -l < $blast_rna`;
die "wc failed: $?" if $?;
chomp($count_l);
my $num_guides = $count_l / 2;

#print $num_guides;

#blast the sgRNA temp file
my $blastresult = "_blast_sgRNA_results.temp";
open my $blast_fh_out, ">>", $blastresult
  or print "unable to open file";

my $tot_factor = ( $tot_mis * 30 ) + 30;

#offtarg arrays
my @off_target_array;
my @off_target_array_;


    my $command =
        "./bowtie -a -n "
      . $seed_mis . " -l "
      . $seed_l . " -e "
      . $tot_factor
      . " -y --quiet ./index/"
      . $blaster . " -f "
      . $blast_rna
      . " -S ./"
      . $blastresult;
    system($command);

    close $out;
    close $temp_out;

    #close $blast_fh_out;
    open $blast_fh_out, "<", $blastresult or print "unable to open file";
    my $gRNA_count = 1;
    my @file_in    = <$blast_fh_out>;
    close $blast_fh_out;
    ##parse blast output and update off-target scores
    for ( my $i = 0 ; $i < $num_guides ; $i++ ) {
        my $regex = "_" . $gRNA_count;
        my $has_match = grep ( /$regex\t/, @file_in );
        push @off_target_array, ( $has_match - 1 );
        $gRNA_count++;
    }


# turn negative numbers into zero
@off_target_array_ = map { $_ < 0 ? 0 : $_ } @off_target_array;

#slurp temp file into array
open( my $fout, "<", "_sgRNA_.temp" )
  or die "Failed to open file: $!\n";
while (<$fout>) {
    chomp;
    push @lines, $_;
}
close $fout;
my @newLines = grep( s/\t/","/g, @lines );

my $output_final = "_sgRNA_candidates".$id.".js";
open my $out_final, ">", $output_final or die;

print $out_final "let dataTableInput = [";

#print Dumper @lines;
for my $i ( 1 .. $#newLines ) {
    my $temp = shift @off_target_array_;
    print $out_final "[\"" . $newLines[$i] . "\",\"" . $temp . "\"],";
}

print $out_final "];";

#empty array
undef @off_target_array_;
undef @off_target_array;

#close file handle
close $out_final;

#remove tmp files
clean_up();

##########################

sub gc_calc {
    my $dna = shift;

    #Count Gs and Cs
    my $countGC = 0;
    while ( $dna =~ /G|C/ig ) {
        $countGC++;
    }

    #Calculate percent GC
    my $percentGC = 100 * $countGC / length($dna);
    return sprintf( "%0.1f", $percentGC );
}

########################################
sub calcPAMP {
    my $sequence = shift;
    my $dH       = 0;
    my $seq_len  = length($sequence);
    my $i;

    # Compute dH
    for ( $i = 0 ; $i <= $seq_len ; $i++ ) {
        my $pair = substr( $sequence, $i, 1 );
        $dH += $PAMHASH{$pair};
    }


    return $dH;
}

########################################
sub gc_freq {
    my $dna     = shift;
    my $countGC = 0;
    while ( $dna =~ /G|C/ig ) {
        $countGC++;
    }
    my $freqGC = $countGC / length($dna);
    return sprintf( "%0.1f", $freqGC );
}

########################################
sub calcFeatures {
    my $dna    = shift;
    my @dnaArr = split //, $dna;
    my @feat   = [];
    my $seq    = substr $dna, 3, 13;
    $feat[0] = gc_freq($seq);
    if ( $dnaArr[19] eq "G" ) {
        push @feat, 1;
    }
    else {
        push @feat, 0;
    }
    if ( $dnaArr[2] eq "T" or $dnaArr[2] eq "A" ) {
        push @feat, 1;
    }
    else {
        push @feat, 0;
    }
    if ( $dnaArr[11] eq "G" or $dnaArr[11] eq "A" ) {
        push @feat, 1;
    }
    else {
        push @feat, 0;
    }
    if ( $dnaArr[5] eq "G" ) {
        push @feat, 1;
    }
    else {
        push @feat, 0;
    }
    if ( $dnaArr[3] eq "T" or $dnaArr[3] eq "A" ) {
        push @feat, 1;
    }
    else {
        push @feat, 0;
    }
    if ( $dnaArr[17] eq "G" or $dnaArr[17] eq "A" ) {
        push @feat, 1;
    }
    else {
        push @feat, 0;
    }
    if ( $dnaArr[4] eq "C" or $dnaArr[4] eq "A" ) {
        push @feat, 1;
    }
    else {
        push @feat, 0;
    }
    if ( $dnaArr[13] eq "G" ) {
        push @feat, 1;
    }
    else {
        push @feat, 0;
    }
    if ( $dnaArr[14] eq "A" ) {
        push @feat, 1;
    }
    else {
        push @feat, 0;
    }
    return (@feat);
}

########################################
sub getScore {
    my $dna          = shift;
    my @model_weight = (
        0.14177385,  0.06966514,  0.04216254,  0.03303432,
        0.02355430,  -0.04746424, -0.04878001, -0.06981921,
        -0.07087756, -0.08160700
    );
    my $model_offset = 0.6505037;
    my @features     = calcFeatures($dna);
    my $len          = scalar(@features);
    my $score        = 0;
    for my $i ( 0 .. $#features ) {
        $score += $features[$i] * $model_weight[$i];
    }
    $score = $score + $model_offset;
    return (sprintf("%.3f", $score));
}


########################################
sub primer3 {
    my $sequence = $_[0];
    my $start = $_[1];
    my $finish = $_[2];
    my $middle = ($start + $finish) / 2 ;
    my $p3       = "primer3_output";
    open my $fh, ">>", $p3 or die;
    close $fh;

    my $new_sequence =
        "SEQUENCE_TEMPLATE="
      . $sequence . "\n"
      . "SEQUENCE_TARGET="
      . $middle . ",8" . "\n"
      . "PRIMER_TASK=generic" . "\n"
      . "PRIMER_PICK_LEFT_PRIMER=1" . "\n"
      . "PRIMER_PICK_RIGHT_PRIMER=1" . "\n"
      . "PRIMER_PRODUCT_SIZE_RANGE="
      . 50 . "-"
      . 1000 . "\n"
      . "PRIMER_SALT_MONOVALENT="
      . 50 . "\n"
      . "PRIMER_GC_CLAMP="
      . 0 . "\n"
      . "PRIMER_MAX_TM="
      . 85 . "\n"
      . "PRIMER_MIN_TM="
      . 30 . "\n"
      . "PRIMER_MIN_SIZE="
      . 18 . "\n"
      . "PRIMER_MAX_SIZE="
      . 35 . "\n"
      . "PRIMER_MAX_GC="
      . 100 . "\n"
      . "PRIMER_MIN_GC="
      . 0 . "\n"
      . "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=./html/cgi-bin/primer3_config/" . "\n"
      . "PRIMER_MAX_NS_ACCEPTED=1" . "\n"
      . "PRIMER_EXPLAIN_FLAG=1" . "\n" . "=";

    my $temp_seq = "primer3_temp.txt";

    open my $p3_fh, ">", $temp_seq or die;
    print $p3_fh($new_sequence);

    my $command = "./html/cgi-bin/primer3_core -format_output < $temp_seq >> $p3";
    system($command);

    my $resultF;
    my $resultR;

    open $fh, "<", $p3 or die;
    my $count = 0;
    while ( my $row = <$fh> ) {

        #print $row;
        if ( $row =~
/LEFT PRIMER\s+\d*\s+\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s(.*)\n/
          )
        {
            $resultF =  $1;
        }
        elsif ( $row =~
/RIGHT PRIMER\s+\d*\s+\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s+\d*\.\d+\s(.*)\n/
          )
        {
            $resultR =  $1;

        }
    }
    my $concatPrimer = "F:".$resultF." R:".$resultR;
    return($concatPrimer);
}

########################################
sub clean_up {
    my $dir = cwd();
    unlink glob "$dir/*.temp";

    #print "\nDirectory has been cleaned.\n";
}

sub slurp {
    my $file = shift;
    open my $fh, '<', $file or die;
    local $/ = undef;
    my $cont = <$fh>;
    close $fh;
    return $cont;
}