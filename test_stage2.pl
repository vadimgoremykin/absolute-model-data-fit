#!/usr/bin/perl -w    
my $MAN = <<TEXT;
test_stage2.pl (v. 1.0) script calculates the test values based on the data produced by test_stage1.pl
script and provides a summary of the test results.

synopsis: test_stage2.pl argument1 argument2 > result

argument1: a file produced using the following command line:
"grep Tvalue *.extension > file", wherein *.extension is a common extension of the
output files produced by the test_stage1.pl script

argument2 (optional): -g 
Specifies the full Gelfand-Ghosh (GG = GGg + GGp) statistic instead of the default GGg statistic
to be used in calculation of substitution model fit. If GG statistic was
not calculated by test_stage1.pl script, the option should not be used.

Note 1: The script requires test_stage1.pl output file names to be in a certain format.
If the empirical model is represented by a biological dataset then the file names should be in
the following format: *simulation_model[MODELNAME]*.extension
If the empirical model(s) is/are represented by a distribution of replicates then the file names should
be in the following format: *simulation_model[MODELNAME1].vs.empirical_model[MODELNAME2]*.extension
wherein "*" is a wildcard character. MODELNAME and MODELNAME1 designate the evolutionary model 
which fit to the observed data should be estimated (SM). MODELNAME2 designate the evolutionary model 
used to generate the observed data (EM).

Each MODELNAME, MODELNAME1 and MODELNAME2 should be 8 characters in length, of which the first 
two are digits (tree code), designating model tree topology and the rest 6 characters designate 
the substitution model used to generate the corresponding replicates, for example: 10GTR+GI

Note 2: In the comparisons involving representation of the observed data by replicates:
- The set of model tree topologies (T set) should be the same for all sets of SMs assuming different
  substitution model specifications. 
- The set of SMs which shares substitution model specification with EM(s) should be included into analysis.
- All EMs should share a distinct substitution model component (e.g. GTR+GI). 
- Each of tree codes should designate the same, distinct tree topology. 
- EM model tree topologies should be present in the T set.
TEXT

$#ARGV>=0 || die ($MAN);
$infile = $ARGV[0];

if ($#ARGV == 1)
{
   if ($ARGV[1] eq "-g")
   {
      $EM_FIT = "Aal";
   }
}
else
{
   $EM_FIT = "Ade";
}


-e $infile || die (sprintf ("cannot find %s file", $infile));
open (VVOD, "<$infile") or die (sprintf ("cannot open %s file for reading", $infile));
while (defined($stroka = <VVOD>))
{
   if ($stroka =~ /simulation_model(..)(......).vs.empirical_model(..)(......).*Tvalue_(...): (.+)\n/)	
   {
      $mode = 2;
      @values = split / /, $6;
      $EM = "$3$4";
      $subst_SM = $2;
      $subst_SM_tree = $1;
      if ($5 eq $EM_FIT)
      {
         $subst_model_fit_all{$EM}{$subst_SM}{$subst_SM_tree} = [ @values ];
      }
   }
   elsif ($stroka =~ /simulation_model(..)(......).*Tvalue_(...): (.+)\n/)
   {
      $mode = 1;
      $subst_SM = $2;
      $subst_SM_tree = $1;
      if ($3 eq $EM_FIT)
      {
         $subst_model_fit_all{$subst_SM}{$subst_SM_tree} = $4;
      }
   }
}      
close VVOD;
if ($mode == 2)
{
   $smfit_defined = 0;
   for $EM (keys %subst_model_fit_all)
   {  
      for $subst_SM (keys %{$subst_model_fit_all{$EM}})
      {  
         $smfit_defined = 1;
         for $EM_replicate (0.. $#values)
	 {
	    $subst_model_fit = 0;$trees = 0;
	    for $tree_SM (keys %{$subst_model_fit_all{$EM}{$subst_SM}})
            {  
               $trees++;
	       $subst_model_fit+=$subst_model_fit_all{$EM}{$subst_SM}{$tree_SM}[$EM_replicate];
	    
            }
	    $subst_model_fit= abs $subst_model_fit/$trees;
	    push @{$subst_model_fit{$EM}{$subst_SM}}, $subst_model_fit;
	 }
      }
   }
}
if ($mode == 1)
{
   for $subst_SM (keys %subst_model_fit_all)
   {  
      $subst_model_fit = 0;$trees = 0;
      for $tree_SM (keys %{$subst_model_fit_all{$subst_SM}})
      {  
         $trees++;
         $subst_model_fit+=$subst_model_fit_all{$subst_SM}{$tree_SM};
      }
      $subst_model_fit= abs $subst_model_fit/$trees;
      $subst_model_fit{$subst_SM} = $subst_model_fit;
   }
}

open (VVOD, "<$infile") or die (sprintf ("cannot open %s file for reading", $infile));
$out_full = "$infile.final_test_values";
if ($mode == 2)
{
   if ($smfit_defined == 1)
   {
      open (VIVODT, ">$out_full") or die (sprintf ("cannot open %s file for writing", $out_full));
   }
}
while (defined($stroka = <VVOD>))
{
   if ($stroka =~ /simulation_model(..)(......).vs.empirical_model(..)(......).*Tvalue_(...): (.+)\n/)	
   {
      if ($5 ne "Aal" and $5 ne "Ade")
      {
      $mode = 2;
      @values = split / /, $6;
      $EM = "$3$4";
      $subst_SM = $2;
      if ($smfit_defined == 1)
      {
         print VIVODT "simulation_model$1$2.vs.empirical_model$3$4:Tvalue_FIN:";
      }
      for $i (0.. $#values)
      {
         if ($smfit_defined == 1)
	 {
	    $values[$i] = $subst_model_fit{$EM}{$subst_SM}[$i]*$values[$i];#BAD
	    print VIVODT " $values[$i]";
	    
	 }
      }
      if ($smfit_defined == 1)
      {
         print VIVODT "\n";
      }
      
      
      
      $hash_R{$5}{$3}{$2}{$1} = [ @values ];
      $mean_value = 0;
      
      for $i (0.. $#values)
      {
	 $mean_value+=$values[$i];
	 $hash_single1{$5}{"$3$4"}{$i}{"$1$2"} = $values[$i];
	 $hash_pairwise_boot{$5}{"$3$4"}{"$1$2"} = 0;
      }
      $mean_value/=($#values + 1);

      $skolko_rep = ($#values + 1);
      $hash{$5}{"$3$4"}{$mean_value} = "$1$2";
      $hashX{$5}{"$3$4"}{$2}{$mean_value} = "$1$2";
      }
   }
   elsif ($stroka =~ /simulation_model(..)(......).*Tvalue_(...): (.+)\n/)
   {
      $mode = 1;
      $hash{$3}{$2}{"$1$2"} = $4;
   }
}
if ($mode == 1)
{
   for $key1 (sort {$a cmp $b} keys %hash)
   {  
      $i=0;
      if ($key1 ne "Aal" and $key1 ne "Ade")
      {
      
         for $key2 (keys %{$hash{$key1}})
         {  
	    for $key3 (keys %{$hash{$key1}{$key2}})
            {  
               $hash_to_print{$key1}{$key3}=   $hash{$key1}{$key2}{$key3}*$subst_model_fit{$key2};
            }
         }
         print "Models sorted in ascending order of test values:\n\n";
         for $curr (sort {$hash_to_print{$key1}{$a} <=> $hash_to_print{$key1}{$b}} keys %{$hash_to_print{$key1}})
         {  
            $i++;
            print "model: $curr test value: $hash_to_print{$key1}{$curr}\n";
            $i == 1 and $output{$key1}{$curr} = $hash_to_print{$key1}{$curr};	 
         }
         print "\n";
      }
   }
   for $key1 (sort {$a cmp $b} keys %output)
   {  
      for $key2 (keys %{$output{$key1}})
      {  
         print "Best model is $key2, test value: $output{$key1}{$key2}\n";
      }
   }
}   

if ($mode == 2)
{



#######################################################################
#these loops compute R values for model sets assuming distinct SUBST model
#######################################################################
for $FNC (keys %hash_R)
{  
   for $CT (keys %{$hash_R{$FNC}})
   {  
      for $SM (keys %{$hash_R{$FNC}{$CT}})
      {
	 for $TT (keys %{$hash_R{$FNC}{$CT}{$SM}})
         {
	     $count_ultimate{$FNC}{"$CT$SM"}{"$TT$SM"} = 0;
	 }
      }
   }    
} 
for $FNC (sort {$a cmp $b} keys %hash_R)
{  
   for $CT (sort {$a cmp $b} keys %{$hash_R{$FNC}})
   {  
      for $SM (sort {$a cmp $b} keys %{$hash_R{$FNC}{$CT}})
      {
	for $TT (sort {$a cmp $b} keys %{$hash_R{$FNC}{$CT}{$SM}})
	 {
	   for $i (0.. $#{$hash_R{$FNC}{$CT}{$SM}{$TT}})
	   {
	      if ($hash_R{$FNC}{$CT}{$SM}{$CT}[$i] < $hash_R{$FNC}{$CT}{$SM}{$TT}[$i])
	      {
		 $count_ultimate{$FNC}{"$CT$SM"}{"$TT$SM"}++;
	      }
	   }
	}
      }
   }	
} 
for $FNC (sort {$a cmp $b} keys %count_ultimate)
{  
   for $CTSM (sort {$a cmp $b} keys %{$count_ultimate{$FNC}})
   {  
      for $key3 (sort {$count_ultimate{$FNC}{$CTSM}{$a} <=> $count_ultimate{$FNC}{$CTSM}{$b}} keys %{$count_ultimate{$FNC}{$CTSM}}  ) 
      {
         $R = ($count_ultimate{$FNC}{$CTSM}{$key3}/$skolko_rep)*100;
	 $R = sprintf ("%.1f", $R);   
	 $count_ultimate{$FNC}{$CTSM}{$key3} = $R;  
      }
   }
}

#######################################################################
#these loops compute R values for separation of each EM from all other models
#######################################################################
for $key1 (sort {$a cmp $b} keys %hash_single1)
{  
   for $key2 (sort {$a cmp $b} keys %{$hash_single1{$key1}})
   {  
      for $key3 (sort {$a <=> $b} keys %{$hash_single1{$key1}{$key2}})
      {
	 for $key4 (sort {$a cmp $b} keys %{$hash_single1{$key1}{$key2}{$key3}})
	 {
	    $test_value = $hash_single1{$key1}{$key2}{$key3}{$key4};
	    $corr_value = $hash_single1{$key1}{$key2}{$key3}{$key2};
	    if ($corr_value <  $test_value)
	    {
	       $hash_pairwise_boot{$key1}{$key2}{$key4}++;
	    }
	 }
      }
   }	
} 
for $key1 (sort {$a cmp $b} keys %hash_pairwise_boot)
{  
   for $key2 (sort {$a cmp $b} keys %{$hash_pairwise_boot{$key1}})
   {  
      for $key3 (sort {$a cmp $b} keys %{$hash_pairwise_boot{$key1}{$key2}}) 
      {

         $bootstrap = ($hash_pairwise_boot{$key1}{$key2}{$key3}/$skolko_rep)*100;
	 $bootstrap = sprintf ("%.1f", $bootstrap);
	 $hash_pairwise_boot{$key1}{$key2}{$key3} = $bootstrap;
      }
   }
}      	 

print "PART 1. Comparison of all simulation models (SMs) to each empirical model (EM)\n\n";

print "MS1 values indicate percentage of times when the SM assuming EM tree topology + EM substitution model\nshowed better fit to EM-based replicates in comparison to any other SM.\n\n";
print "Shown are the mean test values estimated in comparisons of each SM to each EM-based replicate.\n\n";

for $key1 (sort {$a cmp $b} keys %hash)
{  
   $key1 =~ /(...)/;
   for $key2 (sort {$a cmp $b} keys %{$hash{$key1}})
   {  
      for $key3 (sort {$a <=> $b} keys %{$hash{$key1}{$key2}})
      {
         
	 if ($key2 eq $hash{$key1}{$key2}{$key3})
	 {
	     $znach0 = sprintf ("%.10f", $key3);
	     $length = length scalar $znach0;
	     $rest = 30 - $length;
	     $add = "";
	     for $k (1.. $rest){$add.=" ";}
	     
	     print "EM:$key2 SM:$hash{$key1}{$key2}{$key3} test value:$znach0$add MS1: N/A\n";
	 }
	 else
	 {
	     $znach0 = sprintf ("%.10f", $key3); 
	     $length = length scalar $znach0;
	     $rest = 30 - $length;
	     $add = "";
	     for $k (1.. $rest){$add.=" ";}
	     print "EM:$key2 SM:$hash{$key1}{$key2}{$key3} test value:$znach0$add MS1: $hash_pairwise_boot{$key1}{$key2}{$hash{$key1}{$key2}{$key3}}\n";
	 }
      }
      print "\n";
      
   }
}              


######################################################################à
#sorting by subst.model: correct constraint vs. other constraints
#####################################################################à

print "PART 2. Comparison of fit of simulation models sharing common substitution model scheme\nbut assuming different model tree topologies\n\n";

print "MS2 values indicate percentage of times when the SM assuming EM tree topology and a certain\nsubstitution model scheme, (termed 'preferred model') showed better fit to EM-based replicates\nin comparison to any other SM assuming the same substitution model scheme as the preferred model.\n\n";
print "Shown are the mean test values estimated in comparisons of each SM to each EM-based replicate.\n\n";

for $key1 (sort {$a cmp $b} keys %hashX)
{  
   $key1 =~ /(...)/;
   for $key2 (sort {$a cmp $b} keys %{$hashX{$key1}})
   {  
      $key2 =~ m/^(\d\d)/;
      $true_tree = $1;
      
      
      for $key3 (sort {$a cmp $b} keys %{$hashX{$key1}{$key2}})
      {
	 $mean = 0;
	 
	 
	 undef @Rvalues; 
	 for $key4 (sort {$a <=> $b} keys %{$hashX{$key1}{$key2}{$key3}})
         {
            $mean+=$key4;
	    
	    $simulation_model = $hashX{$key1}{$key2}{$key3}{$key4};
	    
	    $znach = sprintf ("%.10f", $key4); 
	    $SM_RT = "$true_tree$key3";
	    
	    
	    $hashX{$key1}{$key2}{$key3}{$key4} =~ m/^(\d\d)/;
            $test_tree = $1;

	    
	    if ($true_tree eq $test_tree)
	    {
	       $length = length scalar $znach;
	       $rest = 30 - $length;
	       $add = "";
	       for $k (1.. $rest){$add.=" ";}
	       print "EM:$key2 SM:$hashX{$key1}{$key2}{$key3}{$key4} test value:$znach$add MS2: N/A\n";
	    }
	    else
	    {
	       $length = length scalar $znach;
	       $rest = 30 - $length;
	       $add = "";
	       for $k (1.. $rest){$add.=" ";}
	       print "EM:$key2 SM:$hashX{$key1}{$key2}{$key3}{$key4} test value:$znach$add MS2: $count_ultimate{$key1}{$SM_RT}{$simulation_model}\n";
	    }
	    if ($hashX{$key1}{$key2}{$key3}{$key4} eq $SM_RT)
	    {
	       $topo_exact_Tvalue_for_SM_RT{$key1}{$key2}{$key3}  =  $znach;
	    }
	    push @Rvalues, $count_ultimate{$key1}{$SM_RT}{$simulation_model};
	 }
         @Rvalues = sort {$a <=> $b} @Rvalues;
         $topo_exact_error{$key1}{$key2}{$key3}  =  $Rvalues[1];
	 print "\n";   
      }
      print "\n";
      
   }
}              


for $key1 (sort {$a cmp $b} keys %hashX)
{  
   for $key2 (sort {$a cmp $b} keys %{$hashX{$key1}})
   {  
      
      $key2 =~ m/^(\d\d)/;
      $true_tree = $1;
      for $key3 (sort {$a cmp $b} keys %{$hashX{$key1}{$key2}})
      {
         $topo_exact{$key1}{$key2}{$key3} = 0;
	 for $key4 (sort {$a <=> $b} keys %{$hashX{$key1}{$key2}{$key3}})
         {
            $hashX{$key1}{$key2}{$key3}{$key4} =~ m/^(\d\d)/;
	    $test_tree = $1;
	    if ($true_tree eq $test_tree)
	    {
	       
	       
	       $topo_exact{$key1}{$key2}{$key3}++;
	    
	    }
	    last;
	 }   
      }
   }
}              
print "PART 3. Summary of the results\n\n";

print "-Subtable 1 shows success (1) or failure (0) of identification of a preferred SM\n- sharing model tree topology with the EM - by the mean test values in comparison to all SMs\n";
print "which share a common substitution model scheme (shown above the subtable) with the preferred SM.\n\n";

print "-Subtable 2 shows the mean test values for the preferred SMs assuming substitution\nmodel components shown above the subtable.\n\n";

print "-Subtable 3 shows the worst model separation (MS2) values between the preferred SM\nand any other SM which assumes the same substitution model scheme (shown above the subtable)\nas the preferred SM, but different model tree topology.\n";
print "See PART2 of the report (above) for the full list of MS2 values.\n\n";

for $key1 (sort {$a cmp $b} keys %topo_exact)#function
{  
   $key1 =~ /(...)/;
   $mean_exact = 0;$mean_exact_times = 0;
   undef %mean_test_CMS;
   for $key2 (sort {$a cmp $b} keys %{$topo_exact{$key1}})
   {  
      for $key3 (sort {$a cmp $b} keys %{$topo_exact{$key1}{$key2}})
      {
          #$rounded = sprintf ("%.3f", $topo_exact_Tvalue_for_SM_RT{$key1}{$key2}{$key3});
	  $mean_test_CMS{$key3}+=$topo_exact_Tvalue_for_SM_RT{$key1}{$key2}{$key3};
      }
   }


   print "Subtable 1\n";
   print "\t\t";
   for $SM (sort { $mean_test_CMS{$a} <=> $mean_test_CMS{$b} } keys %mean_test_CMS)
   {
      print "$SM\t"; 
   }
   print "\n";
   for $key2 (sort {$a cmp $b} keys %{$topo_exact{$key1}})
   {  
      print "EM:$key2\t";
      for $key3 (sort { $mean_test_CMS{$a} <=> $mean_test_CMS{$b} } keys %mean_test_CMS)
      {
          print "$topo_exact{$key1}{$key2}{$key3}\t";
      }
      print "\n";
   }
   print "\n";
   
   print "Subtable 2\n";
   print "\t\t";
   for $SM (sort { $mean_test_CMS{$a} <=> $mean_test_CMS{$b} } keys %mean_test_CMS) 
   {
      print "$SM\t"; 
   }
   print "\n";
   for $key2 (sort {$a cmp $b} keys %{$topo_exact{$key1}})
   {  
      print "EM:$key2\t";
      for $key3 (sort { $mean_test_CMS{$a} <=> $mean_test_CMS{$b} } keys %mean_test_CMS)
      {
          $rounded = sprintf ("%.f", $topo_exact_Tvalue_for_SM_RT{$key1}{$key2}{$key3});
	  
	  print "$rounded\t";
      }
      print "\n";
   }
   print "\n";   

   print "Subtable 3\n";
   print "\t\t";
   for $SM (sort { $mean_test_CMS{$a} <=> $mean_test_CMS{$b} } keys %mean_test_CMS) 
   {
      print "$SM\t"; 
   }
   print "\n";
   for $key2 (sort {$a cmp $b} keys %{$topo_exact{$key1}})
   {  
      print "EM $key2\t";
      for $key3 (sort { $mean_test_CMS{$a} <=> $mean_test_CMS{$b} } keys %mean_test_CMS)
      {
          $mean_exact_times++;
	  $mean_exact+=$topo_exact_error{$key1}{$key2}{$key3};
	  $rounded = sprintf ("%.f", $topo_exact_error{$key1}{$key2}{$key3});
	  
	  print "$rounded\t"; 
      }
      print "\n";
   }
   $mean_exact/=$mean_exact_times;
   $mean_exact= sprintf ("%.2f", $mean_exact);
   print "\nThe mean over values in above subtable is: $mean_exact\n\n";
   $mean_exact{$key1} = $mean_exact;   
}


}
