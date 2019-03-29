#!/usr/bin/perl -w    
use Carp;

my $MAN = <<TEXT;
test_stage1.pl (v. 1.0) script calculates:

-TTC value                          ("Tvalue_TTC" in the output file), 
-Default GGg-based component of TC1 ("Tvalue_GGg" in the output file),
-Full Gelfand-Ghosh statistic       ("Tvalue_GGf" in the output file)
for the test of absolute model-data fit (Goremykin (2019)).

command line synopsis:

-i name of a file containing observed alignment(s) in fasta or sequential phylip format (1 sequence per line)
-r name of a file containing parametric replicates in fasta or sequential phylip format (1 sequence per line)
-t list of filenames of model tree files, in newick format, one tree file name per line
-o output file name. If subsequent analyzes with the test_stage2.pl script should be conducted, please make 
   the output file names compatible with the format that the test_stage2.pl script requires (described when 
   the test_stage2.pl script is called without arguments)
-s (optional) alphabet type
   "-s dna":  DNA alphabet
   "-s prot": Protein alphabet
   "-s ry":   alphabet which has R and Y characters only
   Note: The performance of the test has been tested with DNA data alphabet only. Other options are added for
   exploratory purposes. If the alphabet type is not specified, the script assumes the DNA alphabet by default.
-p (optional) sets a P factor value (default 10,000)
-g (optional) enables calculation of full Gelfand and Ghosh statistic

   
TEXT
$infile_defined = 0;
$repfile_defined = 0;
$treefile_defined = 0;
$outfile_defined = 0;
$seqtype_defined = 0;
$pval_defined = 0;
$fullGG = 0; 

for $i (0..$#ARGV)
{
   if ($ARGV[$i] eq "-i")
   {
      $infile = $ARGV[$i + 1];$infile_defined = 1;
      open (INPUTOBS, "<$infile") or (print "ERROR: cannot open file $infile for reading\n\n" and  print ($MAN) and die);      
   }
   if ($ARGV[$i] eq "-r")
   {
      $bootlist = $ARGV[$i + 1];$repfile_defined = 1;
      open (INPUTREP, "<$bootlist") or (print "ERROR: cannot open file $bootlist for reading\n\n" and  print ($MAN) and die);     
   }
   if ($ARGV[$i] eq "-t")
   {
      $treelist = $ARGV[$i + 1];$treefile_defined = 1;
      open (TREELIST, "<$treelist") or (print "ERROR: cannot open file $treelist for reading\n\n" and  print ($MAN) and die); 
   }
   if ($ARGV[$i] eq "-o")
   {
      $outex = $ARGV[$i + 1];$outfile_defined = 1;
      open (VIVOD, ">$outex") or (print "ERROR: cannot open file $outex for writing\n\n" and  print ($MAN) and die);   
   }
   if ($ARGV[$i] eq "-s")
   {
      $seqtype = $ARGV[$i + 1];
      if ($seqtype eq "dna" or $seqtype eq "prot" or $seqtype eq "ry")
      {
         $seqtype_defined = 1;
      }
      else
      {
         print "ERROR: the sequence type specified ($seqtype) is not implemented\n\n" and print ($MAN) and die;	 
      }
   }
   if ($ARGV[$i] eq "-p")
   {
      if ($ARGV[$i + 1] > 1)
      { 
         $PVAL = $ARGV[$i + 1];
         $pval_defined = 1;
      }
      else
      {
         print "ERROR: the P values <= 1 should not be used\n\n" and print ($MAN) and die;	 
      }
      
   }
   if ($ARGV[$i] eq "-g")
   {
      $fullGG = 1;      
   }



}
if ($pval_defined == 0)
{
   $PVAL = 10000;
}
if ($seqtype_defined == 0)
{
   $seqtype = "dna";
}


unless ( $infile_defined == 1 and $repfile_defined == 1 and $treefile_defined == 1 and $outfile_defined == 1)
{
   $infile_defined == 0   and print "ERROR: file contailing observed alignment(s) has not been specified\n\n" and print ($MAN) and die;
   $repfile_defined == 0  and print "ERROR: file contailing parametric replicates has not been specified\n\n" and print ($MAN) and die;
   $treefile_defined == 0 and print "ERROR: file listing model trees has not been specified\n\n" and print ($MAN) and die;
   $outfile_defined == 0  and print "ERROR: output file name has not been specified\n\n" and print ($MAN) and  die;
}



#######################################################################
#reading observed data (linear phylip or fasta formats)
#######################################################################
while (defined($boot = <INPUTOBS>))
{
   if ($boot =~ /^\s*\d+\s+\d+/)
   {
      $list_mode = 0; $list_format = "phylip";
   }
   elsif ($boot =~ /^\>/)
   {
      $list_mode = 0; $list_format = "fasta";
   }
   else
   {
      $list_mode = 1;
   }
   last;
}
close  INPUTOBS;


$countreplicas1 = 0;
$count_species = 0;
if ($list_mode == 1)
{
   open (INPUTOBS, "<$infile"); 	
   while (defined($boot = <INPUTOBS>))
   {
      $boot =~ /\n/ and chop $boot;
      $countreplicas1++;
      $count_species=0;
      open (VVOD, "<$boot") or croak (sprintf ("cannot open $boot for reading"));
      $read = 0;
      while (defined($string = <VVOD>))
      {
         $read++;
         if ($read == 1)
         {	 
            if ($string =~ /^\s*\d+\s+\d+\n/)
            {
	       while (defined($string = <VVOD>) and $string =~ /^\s*([^\s]+)\s+([^\s]+)\s*\n/ )
	       {
	          $count_species++;
	          $orig_names{$1} = $count_species;
	       }
            }
	    mark1: 
            if ($string =~ /^\>([^\s]+)\s*\n/)
            {
	       $seq = "";
	       $count_species++;
	       $orig_names{$1} = $count_species;
	       while (defined($string = <VVOD>) and $string !~ /^>/ and $string =~ /\n/)
	       {
	       }
	       if (defined $string){$string =~ /^\>.+/ and goto mark1;}
	    }
         }
      }
      close VVOD;   
   }      
   close INPUTOBS;
}

if ($list_mode == 0)
{
   $count_species=0; $read = 0;  $countreplicas1 = 0;
   open (INPUTOBS, "<$infile"); 	
   if ($list_format eq "fasta")
   {  
      while (defined($string = <INPUTOBS>))
      {
         $count_species=0;
	 mark2:
	 if ($string =~ /^\>([^\s]+)\s*\n/)
         {
            $count_species++;
	    $orig_names{$1} = $count_species;
            while (defined($string = <INPUTOBS>) and $string !~ /^>/ and $string =~ /\n/ and $string !~ /^\s*\n/)
            {
            }
            if (defined $string)
	    {
	       if ($string =~ /^\s*\n/)
	       {
	           $count_species = 0;
	       }
	       if ($string =~ /^>/)
	       {
	           goto mark2;
	       }
	    }   
	    else
	    {
	       $count_species = 0;
	    }
         }
      }
   }   
   if ($list_format eq "phylip")
   {
      while (defined($string = <INPUTOBS>))
      {
	 if ($string =~ /^\s*\d+\s+\d+\n/)
	 {
	    $count_species = 0;
	 }
	 if ($string =~ /^\s*([^\s]+)\s+([^\d\s]+)\s*\n/)
	 {
	    $count_species++;
	    $orig_names{$1} = $count_species;
	 }
      }
   }
}      
close INPUTOBS;
$laufnummer = 0;
for $key (sort {$a cmp $b} keys %orig_names)
{
   $laufnummer++;
   $Names_taxa{$key} = $laufnummer;
}



$countreplicas1 = 0;
$count_species = 0;
if ($list_mode == 1)
{
   open (INPUTOBS, "<$infile"); 	
   while (defined($boot = <INPUTOBS>))
   {
      $boot =~ /\n/ and chop $boot;
      $countreplicas1++;
      $count_species=0;
      open (VVOD, "<$boot") or croak (sprintf ("cannot open $boot for reading"));
      $read = 0;
      while (defined($string = <VVOD>))
      {
         $read++;
         if ($read == 1)
         {	 
            if ($string =~ /^\s*\d+\s+\d+\n/)
            {
	       while (defined($string = <VVOD>) and $string =~ /^\s*([^\s]+)\s+([^\s]+)\s*\n/ )
	       {
	          $count_species++;
	          $Xone{$countreplicas1}{$Names_taxa{$1}} = uc $2;
		  $seq = $2;
	       }
            }
	    mark10: 
            if ($string =~ /^\>([^\s]+)\s*\n/)
            {
	       $seq = "";
	       $count_species++;
	       while (defined($string = <VVOD>) and $string !~ /^>/ and $string =~ /\n/)
	       {
	          chop $string;
	          $seq = $seq.$string;
	       }
	       $Xone{$countreplicas1}{$Names_taxa{$1}} = uc $seq;	   

	       if (defined $string){$string =~ /^\>.+/ and goto mark10;}
	    }
         }
      }
      close VVOD;   
   }      
   close INPUTOBS;
}

if ($list_mode == 0)
{
   $count_species=0; $read = 0;  $countreplicas1 = 0;
   open (INPUTOBS, "<$infile"); 	
   if ($list_format eq "fasta")
   {  
      while (defined($string = <INPUTOBS>))
      {
         $count_species=0;
	 mark20:
	 if ($string =~ /^\>([^\s]+)\s*\n/)
         {
            $index_names{$1}++;
	    $seq = "";
            $count_species++;
            while (defined($string = <INPUTOBS>) and $string !~ /^>/ and $string =~ /\n/ and $string !~ /^\s*\n/)
            {
               chop $string;
               $string =~ s/ *$//;
               $seq = $seq.$string;
            }
            $Xone{$index_names{$1}}{$Names_taxa{$1}} = uc $seq;

            if (defined $string)
	    {
	       if ($string =~ /^\s*\n/)
	       {
	           $count_species = 0;
	       }
	       if ($string =~ /^>/)
	       {
	           goto mark20;
	       }
	    }   
	    else
	    {
	       $count_species = 0;
	    }
         }
      }
   }   
   if ($list_format eq "phylip")
   {
      while (defined($string = <INPUTOBS>))
      {
	 if ($string =~ /^\d+\s+\d+\n/)
	 {
	    $count_species = 0;
	    $countreplicas1++;
	 }
	 if ($string =~ /^\s*([^\s]+)\s+([^\d\s]+)\s*\n/)
	 {
	    $count_species++;
	    $Xone{$countreplicas1}{$Names_taxa{$1}} = uc $2;
	    $seq = $2;	       
	 }
      }
   }
}      
close INPUTOBS;

$length_inputaln_orig = length $seq;
$read = 0;
if ($list_mode == 0 and $list_format eq "fasta")
{  
   $count_species = keys %index_names;
   for $key (keys %index_names)
   {
      $countreplicas1 = $index_names{$key};
   }
}


if ($seqtype eq "dna")
{
   $alphabet =  "A T G C";
}
if ($seqtype eq "prot")
{
   $alphabet =  "A R N D C E Q G H I L K M F P S T W Y V";
}
if ($seqtype eq "ry")
{
   $alphabet =  "R Y";
}
@alphabet = split ' ', $alphabet;

#checking for presence of non-specific characters in the observed alignment
if ($countreplicas1 == 1)
{
   for $d (0.. $length_inputaln_orig - 1) 
   {
      $good = 0;
      for  $OTU (sort {$a <=> $b} keys %{$Xone{1}})
      {
         $char1 = substr $Xone{1}{$OTU}, $d, 1;
         for $char2 (@alphabet)
         {
	    if ($char1 eq $char2)
	    {
	       $good++;
	    }
         }
      }   
      if ($good == $count_species)
      {
         push @valid, $d;
      }
      else
      {
         $stellung = $d+1;
         print "Position $stellung in the observed alignment contains character(s) which are not specific for the selected alphabet\n";
	 print "Quitting...\n";
	 die;
      }
   }	  
}
else
{
   for $d (0.. $length_inputaln_orig - 1) 
   {
      push @valid, $d;
   }
}


open (TREELIST, "<$treelist") or croak (sprintf ("cannot open $treelist for reading")); 	
while (defined($treename = <TREELIST>))
{
   chop $treename;
   $counttrees++;
   open (VVOD, "<$treename") or croak (sprintf ("cannot open $tree for reading"));
   while (defined($tree = <VVOD>))
   {
      undef %hash; undef %taxa; undef %clades;
      @tree = split /([()])/, $tree;
      $start = 0;
      $count1 = 0;
      $count2 = 0;
      for $i (0.. $#tree)
      {
         if ($tree[$i] eq "(")
         {
            $start = $i;		 
            $count1 = 0; $count2 = 0; 
            for $ii ($start.. $#tree)
            {
               if ($tree[$ii] eq "(")
               {
                  $count1++;
               }
               if ($tree[$ii] eq ")")
               {
                  $count2++;
	          $end = $ii;
               }
               if ($count1 != 0)
               {
                  if ($count2 != 0)
                  {
                     if ($count1 == $count2)
	             {
	                 $hash{$start} = $end;
		         last;
	             }   
	          }
               }	 
            }
         }
      }
      for $i (sort {$a <=> $b} keys %hash)
      {
         for $ii ($i..  $hash{$i})
         {
            if ($tree[$ii] =~ /[a-zA-Z]/)
            {
               @parts = split /[:,]/, $tree[$ii];	
               for $j (0..  $#parts)
               {
	          if ($parts[$j] =~ /[a-zA-Z]/)
	          {
	             $taxa{$parts[$j]} = "";
	             $clades{$i}{$parts[$j]} = "";
	          }
               }	    
            }
         }
      }
      my $hash_count = keys %taxa;
      for $i (sort {$a <=> $b} keys %clades)
      {
         for $ii (keys %taxa)
         {
            unless (exists $clades{$i}{$ii})
            {
              $new = "-$i";
	      $clades{$new}{$ii} = "";
            }
         }
      }   
      for $i (sort {$a <=> $b} keys %clades)
      {
         $leafs_count = keys %{$clades{$i}};
         if ($hash_count - $leafs_count < 2)
         {
            delete $clades{$i};
         }
         if ($leafs_count == 1)
         {
            delete $clades{$i};
         }
      }   
      for $i (sort {$a <=> $b} keys %clades)
      {	  
         $leafs_count1 = keys %{$clades{$i}};
         for $ii (sort {$a <=> $b} keys %clades)
         {
            if ($i < $ii)
            {  
	       $leafs_count2 = keys %{$clades{$ii}};
  	       if ($leafs_count1 + $leafs_count2 == $hash_count)
  	       {
	          $total_leafs = 0;
	          for $j (keys %taxa)
	          {
	             if (exists $clades{$i}{$j} or exists $clades{$ii}{$j})
	             {
  		        $total_leafs++;
	             }
	          }
  	          if ($total_leafs == $hash_count)
  	          {
	             $count_splits++;
  	             for $k (keys %{$clades{$i}})
  	             {
		        $splits1{$counttrees}{$count_splits}{$Names_taxa{$k}} = "";
			
	             }		  
  	             for $k (keys %{$clades{$ii}})
  	             {
                        $splits2{$counttrees}{$count_splits}{$Names_taxa{$k}} = "";
			
	             }		  
  	          }
	       }
            }
         }
      }   
   }
   close VVOD;   
}      
close TREELIST;
for $tr (sort {$a <=> $b} keys %splits1)
{
   for $spl (keys %{$splits1{$tr}})
   {
      undef @gr;
      for $OTU (keys %{$splits1{$tr}{$spl}})
      {
         push @gr, $OTU;
      }	 
      @gr = sort { $a <=> $b } @gr;
      if ($gr[0] == 1)
      {
         @gr1 = @gr;
         $OTULIST1 = join( ' ' , @gr1) ;
      }
      else
      {
         @gr2 = @gr;
         $OTULIST2 = join( ' ' , @gr2) ;
      }	 
      undef @gr;
      for $OTU (keys %{$splits2{$tr}{$spl}})
      {
	 push @gr, $OTU;
      }
      @gr = sort { $a <=> $b } @gr;
      if ($gr[0] == 1)
      {
         @gr1 = @gr;
         $OTULIST1 = join( ' ' , @gr1) ;
      }
      else
      {
         @gr2 = @gr;
         $OTULIST2 = join( ' ' , @gr2) ;
      }	 
      if (exists $selected_splits{$OTULIST1}{$OTULIST2})
      {

      }
      else
      {
         $counted++;
         $SPLITS1{$counted} = [ @gr1 ];
         $SPLITS2{$counted} = [ @gr2 ];
         $selected_splits{$OTULIST1}{$OTULIST2} = $counted;
         $selected_splits{$OTULIST2}{$OTULIST1} = $counted;
         foreach $OTU1 (@gr1)
         {
	    foreach $OTU2 (@gr1)
	    {
               $name = "$OTU1 $OTU2";
	       if ($OTU1  <  $OTU2)
               {
     	          $members1{$counted}++;
               }         
	    } 
         }
         foreach $OTU1 (@gr2)
         {
	    foreach $OTU2 (@gr2)
	    {
               $name = "$OTU1 $OTU2";
	       if ($OTU1  <  $OTU2)
               {
     	          $members2{$counted}++;
               }         
	    } 
         }
         foreach $OTU1 (@gr1)
         {
            foreach $OTU2 (@gr2)
	    {
               $name = "$OTU1 $OTU2";
     	       $members3{$counted}++;
	    } 
         }
      }
   }
}

for $split (sort {$a <=> $b} keys  %SPLITS1)
{
  $amount1 = 0;
  $poradok1 = 0; 
  foreach $OTU (@{$SPLITS1{$split}})
  {
      $poradok1++;
      $amount1++;
  }
  $amount2 = 0;
  $poradok2 = 0; 
  foreach $OTU (@{$SPLITS2{$split}})
  {
      $poradok2++;
      $amount2++
  }
  $taxon_number1 = $#{$SPLITS1{$split}} + 1;
  $taxon_number2 = $#{$SPLITS2{$split}} + 1;  
  @cathegory =  ($taxon_number1 , $taxon_number2);
  @cathegory = sort {$a <=> $b} @cathegory;
  push @{$split_cat{$cathegory[0]}}, $split;
}



$correct_number_of_internal_branches = $count_species - 3;
for $tr (sort {$a <=> $b} keys %splits1)
{
   $internal_branches = 0;
   for $spl (sort {$a <=> $b} keys %{$splits1{$tr}})
   {
      $internal_branches++;
   }
   if ($internal_branches != $correct_number_of_internal_branches)
   {
      print "Tree parsing error. Tree parser estimated that tree number $tr in the input list has $internal_branches, should be $correct_number_of_internal_branches\n";
      print  "Please use trees in plain newick format, e.g. saved by the SEAVIEW editor\n";
      die;
   }
} 

for $rep (1.. $countreplicas1)
{
   print "handling observed MSA #$rep\n";
   undef %PAT; 
   undef %PAT1; 
   undef @first_patterns;
   undef %FP; 
   for  $OTU (1.. $count_species)
   {
      $Xcurrent{$OTU} = $Xone{$rep}{$OTU};
   }
   $comparisons = 0;
   $rep_length = length $Xcurrent{1};
   for  $OTU1 (1.. $count_species)
   {
      for  $OTU2 (1.. $count_species)
      {
         if ($OTU1 != $OTU2)
	 {
	    $comparisons++;
	 }
      }
   }
   $total_pairs = $comparisons * $rep_length;
   for $d (@valid)
   {
      $pattern = "";
      for $OTU (1.. $count_species)
      {
         $char = substr $Xcurrent{$OTU}, $d, 1;
	 $pattern = $pattern.$char;      
      }
      $PAT{$d} = $pattern;
   }
   for $pos (sort {$a cmp $b} keys %PAT)
   {      
      push @{$PAT1{$PAT{$pos}}}, $pos;
   }
   foreach $pat (  keys %PAT1)
   {         
      push @first_patterns, $PAT1{$pat}[0];
      for $i (0.. $#{$PAT1{$pat}})
      {
         push @{$FP{$PAT1{$pat}[0]}}, $PAT1{$pat}[$i]; 
      }	  
   }
   undef %xnonrandom_prob; undef %xnonrandom_probA; 
   for  $OTU1 (1.. $count_species)
   {
      $fff = $Xcurrent{$OTU1};
      for $pat (keys %PAT1)
      {	  
         $arr_memb = @{$PAT1{$pat}};
         $place = $PAT1{$pat}[0];
     	 $char1 = substr $fff, $place, 1;
	 for  $OTU2 (1.. $count_species)
     	 {
     	    if ($OTU1 < $OTU2)
     	    {
     	       $char2 = substr $Xcurrent{$OTU2}, $place, 1;
   	       for (1.. $arr_memb)
   	       {
     		  $xnonrandom_probA{$char1}{$char2}++;
     	       }   
   	    }
   	 }    
      }
   }
   for $key1 (keys %xnonrandom_probA)
   {  
      for $key2 (keys %{$xnonrandom_probA{$key1}})
      {  
   	 $xnonrandom_prob{$key1}{$key2} = $xnonrandom_probA{$key2}{$key1} + $xnonrandom_probA{$key1}{$key2}; 
      }
   }  
   for $p (0.. $#alphabet)
   {
      for $pp (0.. $#alphabet)
      {
   	 unless (exists $xnonrandom_prob{$alphabet[$p]}{$alphabet[$pp]})
   	 {
   	    $xnonrandom_prob{$alphabet[$p]}{$alphabet[$pp]} = 0;
   	 }
      }   
   }
   for $p (0.. $#alphabet)
   {
      for $pp (0.. $#alphabet)
      {
	 push @{$m_obs{$rep}}, $xnonrandom_prob{$alphabet[$p]}{$alphabet[$pp]};
      }
   }
   undef %prob2char_1;
   for $key1 (sort {$a cmp $b} keys %xnonrandom_prob)
   {  
      for $key2 (sort {$a cmp $b} keys %{$xnonrandom_prob{$key1}})
      {  
   	 $char_aln = "$key1$key2";
     	 $xnonrandom_prob{$char_aln} = $xnonrandom_prob{$key1}{$key2}/$total_pairs;
         delete $xnonrandom_prob{$key1}{$key2};

         if ($key1 ne $key2)
	 {
	    $prob2char_1{$char_aln} = $xnonrandom_prob{$char_aln}/$PVAL;	    
	 }
	 else
	 {
	    $prob2char_1{$char_aln} = $xnonrandom_prob{$char_aln};	      
	 }
      }
   }

   for $cat (sort {$a <=> $b} keys  %split_cat)
   {
      undef %hash_record1_1;
      undef %hash_record2_1;
      for $i (0.. $#{$split_cat{$cat}})
      {
         $split = $split_cat{$cat}[$i];$SCORE_A1 = 0;$SCORE_B1 = 0;
         foreach (@first_patterns)
         {               
	    undef %hash_index1; 
	    undef %hash_index2;
	    $split_signature1 = "";
	    $split_signature2 = "";
	    foreach $OTU (@{$SPLITS1{$split}})
	    {
	       $char = substr $Xcurrent{$OTU}, $_, 1;
	       $hash_index1{$char}++;
	    }
	    foreach $OTU (@{$SPLITS2{$split}})
	    {
	       $char = substr $Xcurrent{$OTU}, $_, 1;
	       $hash_index2{$char}++;
	    }
            for $char (@alphabet)
            {
               if (exists $hash_index1{$char})
	       {
	          $split_signature1.=" $hash_index1{$char}";
	       }
	       else
	       {
	          $split_signature1.=" 0";
	       }   
            }
            for $char (@alphabet)
            {
               if (exists $hash_index2{$char})
	       {
	          $split_signature2.=" $hash_index2{$char}";
	       }
	       else
	       {
	          $split_signature2.=" 0";
	       }   
            }
            singular();
            foreach ( @{$FP{$_}} )
	    {
	       $SCORE_A1+=$SCORE_A;
   	       $SCORE_B1+=$SCORE_B;
	    }
         }
	 $score_full1{"TTC"}{$rep}{"$split#1"} = $SCORE_A1/$rep_length;	  
	 $score_full1{"TTC"}{$rep}{"$split#2"} = $SCORE_B1/$rep_length;
         $SCORE_A1 = 0;  $SCORE_B1 = 0;
      }
   }
}
undef %SCORE_A1;  undef %SCORE_B1; undef %SCORE_A; undef %SCORE_B;  


#################################################################
#reading replicates
#################################################################

while (defined($boot = <INPUTREP>))
{
   if ($boot =~ /^\s*\d+\s+\d+/)
   {
      $list_mode = 0; $list_format = "phylip";
   }
   elsif ($boot =~ /^\>/)
   {
      $list_mode = 0; $list_format = "fasta";
   }
   else
   {
      $list_mode = 1;
   }
   last;
}
close  INPUTREP;

$countreplicas = 0;
$count_species = 0;
if ($list_mode == 1)
{
   open (INPUTREP, "<$bootlist"); 	
   while (defined($boot = <INPUTREP>))
   {
      $boot =~ /\n/ and chop $boot;
      $countreplicas++;
      $count_species=0;
      open (VVOD, "<$boot") or croak (sprintf ("cannot open $boot for reading"));
      $read = 0;
      while (defined($string = <VVOD>))
      {
         $read++;
         if ($read == 1)
         {	 
            if ($string =~ /^\s*\d+\s+\d+\n/)
            {
	       while (defined($string = <VVOD>) and $string =~ /^\s*([^\s]+)\s+([^\s]+)\s*\n/ )
	       {
	          $count_species++;
	          $Xtest{$countreplicas}{$Names_taxa{$1}} = $2;
		  $seq = $2;
	       }
            }
	    mark100: 
            if ($string =~ /^\>([^\s]+)\s*\n/)
            {
	       $seq = "";
	       $count_species++;
	       while (defined($string = <VVOD>) and $string !~ /^>/ and $string =~ /\n/)
	       {
	          chop $string;
	          $seq = $seq.$string;
	       }
	       $Xtest{$countreplicas}{$Names_taxa{$1}} = $seq;	   

	       if (defined $string){$string =~ /^\>.+/ and goto mark100;}
	    }
         }
      }
      close VVOD;   
   }      
   close INPUTREP;
}

if ($list_mode == 0)
{
   $count_species=0; $read = 0;  $countreplicas = 0;
   open (INPUTREP, "<$bootlist"); 	
   if ($list_format eq "fasta")
   {  
      while (defined($string = <INPUTREP>))
      {
         $count_species=0;
	 mark200:
	 if ($string =~ /^\>([^\s]+)\s*\n/)
         {
            $index_names{$1}++;
	    $seq = "";
            $count_species++;
            while (defined($string = <INPUTREP>) and $string !~ /^>/ and $string =~ /\n/ and $string !~ /^\s*\n/)
            {
               chop $string;
               $string =~ s/ *$//;
               $seq = $seq.$string;
            }
            $Xtest{$index_names{$1}}{$Names_taxa{$1}} = $seq;

            if (defined $string)
	    {
	       if ($string =~ /^\s*\n/)
	       {
	           $count_species = 0;
	       }
	       if ($string =~ /^>/)
	       {
	           goto mark200;
	       }
	    }   
	    else
	    {
	       $count_species = 0;
	    }
         }
      }
   }   
   if ($list_format eq "phylip")
   {
      while (defined($string = <INPUTREP>))
      {
	 if ($string =~ /^\s*\d+\s+\d+\n/)
	 {
	    $count_species = 0;
	    $countreplicas++;
	 }
	 if ($string =~ /^\s*([^\s]+)\s+([^\d\s]+)\s*\n/)
	 {
	    $count_species++;
	    $Xtest{$countreplicas}{$Names_taxa{$1}} = $2;
	    $seq = $2;
	 }
      }
   }
}      
close INPUTREP;


if ($list_mode == 0 and $list_format eq "fasta")
{  
   $count_species = keys %index_names;
   for $key (keys %index_names)
   {
      $countreplicas = $index_names{$key};
   }
}



for $rep (1.. $countreplicas)
{
   print "handling TM-based MSA #$rep\n";
   undef %PAT; 
   undef %PAT1; 
   undef @first_patterns;
   undef %FP; 
   for  $OTU (1.. $count_species)
   {
      $Xcurrent{$OTU} = uc $Xtest{$rep}{$OTU};
   }
   $rep_length = length $Xcurrent{1};
   $comparisons = 0;
   for  $OTU1 (1.. $count_species)
   {
      for  $OTU2 (1.. $count_species)
      {
         if ($OTU1 != $OTU2)
	 {
	    $comparisons++;
	 }
      }
   }
   $total_pairs = $comparisons * $rep_length;
   for $d (@valid)
   {
      $pattern = "";
      for $OTU (1.. $count_species)
      {
         $char = substr $Xcurrent{$OTU}, $d, 1;
	 $pattern = $pattern.$char;      
      }
      $PAT{$d} = $pattern;
   }
   for $pos (sort {$a cmp $b} keys %PAT)
   {      
      push @{$PAT1{$PAT{$pos}}}, $pos;
   }
   foreach $pat (  keys %PAT1)
   {         
      push @first_patterns, $PAT1{$pat}[0];
      for $i (0.. $#{$PAT1{$pat}})
      {
         push @{$FP{$PAT1{$pat}[0]}}, $PAT1{$pat}[$i]; 
      }	  
   }
   undef %xnonrandom_prob; undef %xnonrandom_probA; 
   for  $OTU1 (1.. $count_species)
   {
      $fff = $Xcurrent{$OTU1};
      for $pat (keys %PAT1)
      {	  
         $arr_memb = @{$PAT1{$pat}};
         $place = $PAT1{$pat}[0];
     	 $char1 = substr $fff, $place, 1;
	 for  $OTU2 (1.. $count_species)
     	 {
     	    if ($OTU1 < $OTU2)
     	    {
     	       $char2 = substr $Xcurrent{$OTU2}, $place, 1;
   	       for (1.. $arr_memb)
   	       {
     		  $xnonrandom_probA{$char1}{$char2}++;
     	       }   
   	    }
   	 }    
      }
   }
   for $key1 (keys %xnonrandom_probA)
   {  
      for $key2 (keys %{$xnonrandom_probA{$key1}})
      {  
   	 $xnonrandom_prob{$key1}{$key2} = $xnonrandom_probA{$key2}{$key1} + $xnonrandom_probA{$key1}{$key2}; 
      }
   }  
   for $p (0.. $#alphabet)
   {
      for $pp (0.. $#alphabet)
      {
   	 unless (exists $xnonrandom_prob{$alphabet[$p]}{$alphabet[$pp]})
   	 {
   	    $xnonrandom_prob{$alphabet[$p]}{$alphabet[$pp]} = 0;
   	 }
      }   
   }
   for $p (0.. $#alphabet)
   {
      for $pp (0.. $#alphabet)
      {
	 push @{$m_rep{$rep}}, $xnonrandom_prob{$alphabet[$p]}{$alphabet[$pp]};
      }
   }
   undef %prob2char_1;
   for $key1 (sort {$a cmp $b} keys %xnonrandom_prob)
   {  
      for $key2 (sort {$a cmp $b} keys %{$xnonrandom_prob{$key1}})
      {  
   	 $char_aln = "$key1$key2";
     	 $xnonrandom_prob{$char_aln} = $xnonrandom_prob{$key1}{$key2}/$total_pairs;
         delete $xnonrandom_prob{$key1}{$key2};

         if ($key1 ne $key2)
	 {
	    $prob2char_1{$char_aln} = $xnonrandom_prob{$char_aln}/$PVAL;	    
	 }
	 else
	 {
	    $prob2char_1{$char_aln} = $xnonrandom_prob{$char_aln};	      
	 }
      }
   }

   for $cat (sort {$a <=> $b} keys  %split_cat)
   {
      undef %hash_record1_1;
      undef %hash_record2_1;
      for $i (0.. $#{$split_cat{$cat}})
      {
         $split = $split_cat{$cat}[$i];$SCORE_A1 = 0;$SCORE_B1 = 0;
         foreach (@first_patterns)
         {               
	    undef %hash_index1; 
	    undef %hash_index2;
	    $split_signature1 = "";
	    $split_signature2 = "";
	    foreach $OTU (@{$SPLITS1{$split}})
	    {
	       $char = substr $Xcurrent{$OTU}, $_, 1;
	       $hash_index1{$char}++;
	    }
	    foreach $OTU (@{$SPLITS2{$split}})
	    {
	       $char = substr $Xcurrent{$OTU}, $_, 1;
	       $hash_index2{$char}++;
	    }
            for $char (@alphabet)
            {
               if (exists $hash_index1{$char})
	       {
	          $split_signature1.=" $hash_index1{$char}";
	       }
	       else
	       {
	          $split_signature1.=" 0";
	       }   
            }
            for $char (@alphabet)
            {
               if (exists $hash_index2{$char})
	       {
	          $split_signature2.=" $hash_index2{$char}";
	       }
	       else
	       {
	          $split_signature2.=" 0";
	       }   
            }
            singular();
            foreach ( @{$FP{$_}} )
	    {
	       $SCORE_A1+=$SCORE_A;
   	       $SCORE_B1+=$SCORE_B;
	    }
         }
	 $RESULT_REP3_LINK{"TTC $split#1 $rep"} = $SCORE_A1/$rep_length;  
	 $RESULT_REP3_LINK{"TTC $split#2 $rep"} = $SCORE_B1/$rep_length;  
	 $RESULT_REP4_LINK{"TTC $split#1"}++; 
	 $RESULT_REP4_LINK{"TTC $split#2"}++; 
         $SCORE_A1 = 0;  $SCORE_B1 = 0;
      }
   }
}

###############################
for $rep (keys %m_obs)
{
   for $i (0.. $#{$m_obs{$rep}})
   {
      if ($m_obs{$rep}[$i] == 0)
      {
         print  "msa #$rep representing empirical model violates GGg/GG requirement (positive values of all counts of alignments of character states)\n";
         die;
      }
   }
}   

if ($fullGG == 1)
{
   for $rep (keys %m_rep)
   {
      for $i (0.. $#{$m_rep{$rep}})
      {
         if ($m_rep{$rep}[$i] == 0)
         {
            print  "msa #$rep representing simulation model violates GG requirement (positive values of all counts of alignments of character states)\n";
            die;
         }
      }
   }   
}   
else
{
   for $i (0.. $#{$m_rep{1}})
   {
      $passt = 0;
      for $rep (keys %m_rep)
      {
         if ($m_rep{$rep}[$i] > 0)
         {
            $passt = 1;
         }
      }
      if ($passt == 0)
      {
         print "distribution of replicates representing simulation model violates GGg requirement\n(positive value of each count of alignments of character states should be encountered at least once in any replicate)\n";
         die;
      }
   }   
}   
   
   
   
   
#   for $rep1 (keys %m_rep)
#   {
#      if ($#{$m_obs{$rep}} != $#{$m_rep{$rep1}})
#      {
#         print  "requirement of co-occurrence of all alignments among character in the datasets representing SM and TM is violated\n";
#         die;
#      }
#   }


for $rep (keys %m_obs)
{
   if ($rep == 1)
   {
      for $i (0.. $#{$m_obs{$rep}}) #here
      {
         push @bins, $i;
      }  
   }
   for $i (0.. $#{$m_obs{$rep}})
   {
      $val = $m_obs{$rep}[$i];
      $number_sites_obs{$rep}+=$val;
      $BIN{$rep}{$i} = $val;
   }
}
$np = 0;
$nrep = $countreplicas;
for $rep (keys %m_rep)
{
   $np++;

   for $i (0.. $#{$m_rep{$rep}})
   {
      $val = $m_rep{$rep}[$i];
      $BIN_RT{$i}+=$val;
      $number_sites_rep{$rep}+=$val;
      $BIN_RS{$rep}{$i} = $val;
   }
}


for $rep1 (keys %number_sites_obs)
{
   for $rep2 (keys %number_sites_rep)
   {
      if ($number_sites_obs{$rep1} != $number_sites_rep{$rep2})
      {
         print  "Alignments representing SM and TM have different lengths\n";
         die;
      }
      else
      {
        $ns = $number_sites_obs{$rep1};
      }
   }
}





for $bin (keys %BIN_RT)
{
   $BIN_RT{$bin} = $BIN_RT{$bin}/$countreplicas;
}
for $rep0 (sort {$a <=> $b} keys %BIN)
{
   undef %BIN_REF;
   for $bin (keys %{$BIN{$rep0}})
   {
      $BIN_REF{$bin} = ($BIN_RT{$bin} + $BIN{$rep0}{$bin})/2;
   }


   $t_m = 0;
   for $i (0.. $#bins)
   {
      $part1_0 = $BIN_RT{$bins[$i]} * log ($BIN_RT{$bins[$i]});
      $t_m+=$part1_0;
   }
   $t_m = -log($ns) + (1/$ns)*$t_m;

   $t_y = 0;
   for $i (0.. $#bins)
   {
      $part1_0 = $BIN{$rep0}{$bins[$i]} * log ($BIN{$rep0}{$bins[$i]});
      $t_y+=$part1_0;
   }
   $t_y = -log($ns) + (1/$ns)*$t_y;

   $t_ref = 0;
   for $i (0.. $#bins)
   {
      $part1_0 = $BIN_REF{$bins[$i]} * log ($BIN_REF{$bins[$i]});
      $t_ref+=$part1_0;
   }
   $t_ref = -log($ns) + (1/$ns)*$t_ref;


   $GGg = 2*$ns*2*(  ($t_m + $t_y)/2 - $t_ref);
   push @{$REPORT{GGg}}, $GGg;


   if ($fullGG == 1)
   {
      $part2 = 0;
      for $replicas (1.. $countreplicas)
      {
         $part1 = 0;
         for $i (0.. $#bins)
         {
            $part1_0 = $BIN_RS{$replicas}{$bins[$i]} * log ($BIN_RS{$replicas}{$bins[$i]});
            $part1+=$part1_0;
      
         }
         $part1 = -log($ns) + (1/$ns)*$part1;
         $part2+=$part1;
      }
      $GGp = 2*$ns*(  ((1/$np)*$part2)  - $t_m );
      $GG = $GGp + $GGg;
      push @{$REPORT{GGf}}, $GG;
   }
}

for $FNC (keys %score_full1)
{
   for $rep0 (sort {$a <=> $b} keys %{$score_full1{$FNC}}  )
   {
      $mean_zsp0 = 0;
      for $rep (1.. $nrep) 
      {
	 undef @observ_for_regress;
         undef @replic_for_regress;
	 for $link (keys %{$score_full1{$FNC}{$rep0}}) 
         {
	    $obs_value = $score_full1{$FNC}{$rep0}{$link};
            $test_value = $RESULT_REP3_LINK{"$FNC $link $rep"};
	    push @observ_for_regress,  $obs_value;
            push @replic_for_regress,  $test_value;
         }
         rank();
	 if ($spearman > 0.999999999999999)
	 {
	    $z_rnk = 0.5*(log((1+0.999999999999999)/(1-0.999999999999999)));
	 }
	 else
	 {
	    $z_rnk = 0.5*(log((1+$spearman)/(1-$spearman)));
	 }
         $mean_zsp0+=$z_rnk;
      }
      $mean_zsp0/=$nrep;
      $mean_zsp0 = (exp(2*$mean_zsp0) - 1)/(exp(2*$mean_zsp0) + 1);
      push @{$REPORT{$FNC}}, 1-$mean_zsp0;
   }       
} 


for $FNC (sort {$a cmp $b} keys %REPORT)
{
   print VIVOD "Tvalue_$FNC:";
   for $i (0.. $#{$REPORT{$FNC}})
   {
      print VIVOD " $REPORT{$FNC}[$i]";
   }
   print VIVOD "\n";
}

################################################################################à
#Spearman
#################################################################################

sub rank
{

   undef  %hash_ranks1; undef %hash_dist1; 
   undef  %hash_ranks2; undef %hash_dist2; 
   undef @tied;undef @tied1;undef @oboima;$krugi = 0;$sum_of_ranks = 0;
   push @oboima, 0;
   for $i (0.. $#replic_for_regress)
   {
      $hash_dist1{$i} = $replic_for_regress[$i];
   }
   for $i (sort {$hash_dist1{$b} <=> $hash_dist1{$a}} keys %hash_dist1)
   {
      $krugi++;
      if ($krugi > 1 and $hash_dist1{$i} eq $hash_dist1{$previous})
      {
           if ($oboima[$#oboima] == 0)
           {
              push @tied, $previous; push @tied1, $krugi-1;
              push @tied, $i;	     push @tied1, $krugi;
           }
           if ($oboima[$#oboima] == 1 and $tied[$#tied] == $previous)
           {
              push @tied, $i; push @tied1, $krugi;
           }
           push @oboima, 1;
      }
      else
      {
         if ($oboima[$#oboima] == 1)
         {
            $sum_of_ranks = 0;
            for $ii (0.. $#tied1)
            {
               $sum_of_ranks =  $sum_of_ranks + $tied1[$ii];
            }
            for $ii (0.. $#tied)
            {
               $hash_ranks1{$tied[$ii]} = $sum_of_ranks/($#tied+1);
            }
            undef @tied;undef @tied1;
         }
         push @oboima, 0;
         $hash_ranks1{$i} = $krugi;
      }
      $previous = $i;
   }
   if ($oboima[$#oboima] == 1)
   {
      $sum_of_ranks = 0;
      for $ii (0.. $#tied1)
      {
         $sum_of_ranks =  $sum_of_ranks + $tied1[$ii];
      }
      for $ii (0.. $#tied)
      {
         $hash_ranks1{$tied[$ii]} = $sum_of_ranks/($#tied+1);
      }
      undef @tied;undef @tied1;
   }
   $krugi = 0;
   
   undef @tied;undef @tied1;undef @oboima;$krugi = 0;$sum_of_ranks = 0;
   push @oboima, 0;

   for $i (0.. $#observ_for_regress)
   {
      $hash_dist2{$i} = $observ_for_regress[$i];
   }
   for $i (sort {$hash_dist2{$b} <=> $hash_dist2{$a}} keys %hash_dist2)
   {
      $krugi++;
      if ($krugi > 1 and $hash_dist2{$i} eq $hash_dist2{$previous})
      {
           if ($oboima[$#oboima] == 0)
           {
              push @tied, $previous; push @tied1, $krugi-1;
              push @tied, $i;	     push @tied1, $krugi;
           }
           if ($oboima[$#oboima] == 1 and $tied[$#tied] == $previous)
           {
              push @tied, $i; push @tied1, $krugi;
           }
           push @oboima, 1;
      }
      else
      {
          if ($oboima[$#oboima] == 1)
          {
             $sum_of_ranks = 0;
             for $ii (0.. $#tied1)
             {
        	$sum_of_ranks =  $sum_of_ranks + $tied1[$ii];
             }
             for $ii (0.. $#tied)
             {
        	$hash_ranks2{$tied[$ii]} = $sum_of_ranks/($#tied+1);
             }
             undef @tied;undef @tied1;
          }
          push @oboima, 0;
          $hash_ranks2{$i} = $krugi;
      }
      $previous = $i;
   }
   if ($oboima[$#oboima] == 1)
   {
      $sum_of_ranks = 0;
      for $ii (0.. $#tied1)
      {
         $sum_of_ranks =  $sum_of_ranks + $tied1[$ii];
      }
      for $ii (0.. $#tied)
      {
         $hash_ranks2{$tied[$ii]} = $sum_of_ranks/($#tied+1);
      }
      undef @tied;undef @tied1;
   }
   $krugi = 0;  	   
   $sumX = 0;$sumY = 0;$sumXY = 0;$sumXX = 0;$sumYY = 0;
   $nr1 = $#replic_for_regress+1;
   for $i (0.. $#replic_for_regress)
   {
      $sumX = $sumX + $hash_ranks1{$i};
      $sumY = $sumY + $hash_ranks2{$i};
      $XY = $hash_ranks1{$i} * $hash_ranks2{$i};
      $sumXY = $sumXY + $XY;
      
      $XX = $hash_ranks1{$i} * $hash_ranks1{$i};
      $sumXX = $sumXX + $XX;
      
      $YY = $hash_ranks2{$i} * $hash_ranks2{$i};
      $sumYY = $sumYY + $YY;
   }
   if (($nr1*$sumXX - $sumX*$sumX) >= 0)
   {
      $toroot1 = $nr1*$sumXX - $sumX*$sumX;
   }
   else
   {
      $toroot1 = 0 - ($nr1*$sumXX - $sumX*$sumX);
   }
   if (($nr1*$sumYY - $sumY*$sumY) >= 0)
   {
      $toroot2 = $nr1*$sumYY - $sumY*$sumY;
   }
   else
   {
      $toroot2 = 0 - ($nr1*$sumYY - $sumY*$sumY);
   }
   $spearman = ($nr1*$sumXY - $sumX*$sumY)/(sqrt($toroot1)*sqrt($toroot2));
}


sub singular
{
   if (exists $hash_record1_1{$split_signature1})							     	   
   {													     	   
      $result_input1_1 = $hash_record1_1{$split_signature1};						     	   
   }													     	   
   else 												     	   
   {													     	   
      undef $result_input1_1; 										     	   
      foreach  $OTU1 (@{$SPLITS1{$split}})								     	   
      { 												     	   
   	 foreach  $OTU2 (@{$SPLITS1{$split}})								     	   
   	 {												     	   
   	    if ($OTU1  <  $OTU2)				     	   
   	    {												     	   
   	       $char1 = substr $Xcurrent{$OTU1}, $_, 1; 						     	  
   	       $char2 = substr $Xcurrent{$OTU2}, $_, 1; 	   					     	  
   	       $char_aln = "$char1$char2";
   	       $result_input1_1+=$prob2char_1{$char_aln};						     	  
   	    }												     	  
   	 }												     	   
      } 												     	   
      $hash_record1_1{$split_signature1} = $result_input1_1;						     	   
   }													     	   
   if (exists $hash_record1_1{$split_signature2})								   
   {													     	   
      $result_input2_1 = $hash_record1_1{$split_signature2};						     	   
   }													     	   
   else 												     	   
   {													     	   
      undef $result_input2_1; 										     	   
      foreach  $OTU1 (@{$SPLITS2{$split}})								     	   
      { 												     	   
   	 foreach  $OTU2 (@{$SPLITS2{$split}})								     	   
   	 {												     	   
   	    if ($OTU1  <  $OTU2) 						     	   
   	    {												     	   
   	       $char1 = substr $Xcurrent{$OTU1}, $_, 1; 						     	  
   	       $char2 = substr $Xcurrent{$OTU2}, $_, 1; 						     	  
   	       $char_aln = "$char1$char2";
   	       $result_input2_1+=$prob2char_1{$char_aln};						     	  
   	    }												     	  
   	 }												     	   
      } 												     	   
      $hash_record1_1{$split_signature2} = $result_input2_1;						     	   
   }													     	   
   $cypher0 = "$split_signature1#$split_signature2";							     	   
   if (exists $hash_record2_1{$cypher0})
   {													     	   
      $result_input3_1 = $hash_record2_1{$cypher0};							     	   
   }													     	   
   else 												     	   
   {													     	   
      undef $result_input3_1; 										     	   
      foreach  $OTU1 (@{$SPLITS1{$split}})								     	   
      { 												     	   
   	 foreach  $OTU2 (@{$SPLITS2{$split}})								     	   
   	 {												     	   
   	    $char1 = substr $Xcurrent{$OTU1}, $_, 1;							     	  
   	    $char2 = substr $Xcurrent{$OTU2}, $_, 1;							     	  
   	    $char_aln = "$char1$char2";
   	    $result_input3_1+=$prob2char_1{$char_aln};							     	  
   	 }												     	   
      } 												     	   
      $cypher1 = "$split_signature1#$split_signature2"; 						     	   
      $cypher2 = "$split_signature2#$split_signature1"; 						     	   
      $hash_record2_1{$cypher1} = $result_input3_1;								   
      $hash_record2_1{$cypher2} = $result_input3_1;								   
   }													     	   
   $memb1 = $members1{$split};										     	   
   $memb2 = $members2{$split};										     	   
   $memb3 = $members3{$split};										     	   
   $result_input1_1/=$memb1;											   
   $result_input2_1/=$memb2;											   
   $result_input3_1/=$memb3;											   


   $value1_1 = $result_input1_1/$result_input3_1;								 
   $value2_1 = $result_input2_1/$result_input3_1;								 


   $SCORE_A  = $value1_1;										   
   $SCORE_B  = $value2_1;									     	   
}
