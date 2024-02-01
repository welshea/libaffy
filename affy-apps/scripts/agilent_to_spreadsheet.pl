#!/usr/bin/perl -w

#
# Convert Agilent .txt files into a spreadsheet of intensities.
# Sets bad spots (bg-subtracted signal <= -100) to 0.
# If a spot is bad in one channel, flag it as bad in the other channel too.
# The minimum of the remaining "good" intensities is then used to cause all
#  points to be shifted so that the minimum "good" signal -> 1.
#
# When mean signal is not available, presumably due to bugs in certain versions
#  of the Agilent scanning software (?), median signal is used instead.
#
# iron_generic, part of the libaffy software package, implements an RMA
#  background subtraction that treats 0's as missing values and handles them
#  appropriately.  I would recommend using RMA background subtraction on
#  the intensities that result from this perl script that you are currently
#  reading.
#
#  It looks like on 2011-12-26, ArrayExpress automatically converted
#   deposited Agilent raw files into some mis-formatted header nonsense.
#    Example accessions: E-TABM-337, E-TABM-386
#
#  It looks like a different corruption took place on 2012-02-04 on at least
#   one accession I have looked at: E-MTAB-271
#
#  I've added some code to try to deal with this....


use Scalar::Util qw(looks_like_number);


sub is_number
{
    # use what Perl thinks is a number first
    # this is purely for speed, since the more complicated REGEX below should
    #  correctly handle all numeric cases
    if (looks_like_number($_[0]))
    {
        # Perl treats infinities as numbers, Excel does not.
        #
        # Perl treats NaN or NaNs, and various mixed caps, as numbers.
        # Weird that not-a-number is a number... but it is so that
        # it can do things like nan + 1 = nan, so I guess it makes sense
        #
        if ($_[0] =~ /^[-+]*(Inf|NaN)/i)
        {
            return 0;
        }
        
        return 1;
    }

    # optional + or - sign at beginning
    # then require either:
    #  a number followed by optional comma stuff, then optional decimal stuff
    #  mandatory decimal, followed by optional digits
    # then optional exponent stuff
    #
    # Perl cannot handle American comma separators within long numbers.
    # Excel does, so we have to check for it.
    # Excel doesn't handle European dot separators, at least not when it is
    #  set to the US locale (my test environment).  I am going to leave this
    #  unsupported for now.
    #
    if ($_[0] =~ /^([-+]?)([0-9]+(,[0-9]{3,})*\.?[0-9]*|\.[0-9]*)([Ee]([-+]?[0-9]+))?$/)
    {
        # current REGEX can treat '.' as a number, check for that
        if ($_[0] eq '.')
        {
            return 0;
        }
        
        return 1;
    }
    
    return 0;
}



$subtract_local_background_flag = 1;

$min_signal_cutoff = 1E-5;
$read_a_file_flag  = 0;

# read in command line arguments
for ($f = 0; $f < @ARGV; $f++)
{
    $filename = $ARGV[$f];
    
    # strip .txt file extension
    $filename_root = $filename;
    while ($filename_root =~ /\.txt$/i)
    {
        $filename_root =~ s/\.txt$//i;
    }
    
    # strip leading path
    while ($filename_root =~ /.*\//)
    {
        $filename_root =~ s/.*\///;
    }
    while ($filename_root =~ /.*\\/)
    {
        $filename_root =~ s/.*\\//;
    }
    
    $cy3_name = sprintf "%s_cy3", $filename_root;
    $cy5_name = sprintf "%s_cy5", $filename_root;

    open INFILE,  "$filename"  or die "can't open $filename\n";

    # skip down to header line
    while(defined($line=<INFILE>))
    {
        if ($line =~ /^FEATURES/)
        {
            last;
        }
        
        # weird file format from E-MTAB-271
        # it looks like ArrayExpress now imports the raw files into a database,
        #  then re-exports them as corrupted format files?
        if ($line =~ /CompositeSequence Identifier/)
        {
            last;
        }
        
        # other weird ArrayExpress-corrupted files
        #  ex: E-TABM-337, E-TABM-386
        #
        # I know ArrayExpress corrupted these, because I uploaded them myself
        #  before their current date stamps.
        if ($line =~ /^\s*metaColumn/)
        {
            last;
        }
        
        # Just in case ArrayExpress has yet other corrupted file formats,
        #  I'll check for this string too, which seems to be in common between
        #  the above corrupted cases
        if ($line =~ /\bFeature Extraction Software:/)
        {
            last;
        }
    }
    
    # header line
    %agilent_header_hash = ();
    $line =~ s/[\r\n]+//;
    $line =~ s/\"//g;
    @array = split /\t/, $line;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\s+/ /g;
        
        if ($array[$i] =~ /[A-Za-z0-9]/)
        {
            # strip the junk from the beginning of E-MTAB-271 files
            $array[$i] =~ s/^Feature Extraction Software://;
        
            $agilent_header_hash{$array[$i]} = $i;
        }
    }

    $feature_num_col  = $agilent_header_hash{'FeatureNum'};
    $probe_name_col   = $agilent_header_hash{'ProbeName'};
    $control_type_col = $agilent_header_hash{'ControlType'};
    $g_mean_col       = $agilent_header_hash{'gMeanSignal'};
    $g_bg_mean_col    = $agilent_header_hash{'gBGMeanSignal'};
    $r_mean_col       = $agilent_header_hash{'rMeanSignal'};
    $r_bg_mean_col    = $agilent_header_hash{'rBGMeanSignal'};
    
    if (!defined($feature_num_col))
    {
        die "ERROR -- FeatureNum column not found in $filename\n";
    }
    if (!defined($probe_name_col))
    {
        die "ERROR -- ProbeName column not found in $filename\n";
    }
    
    # Check for bugged Agilent software circa 2009 or so
    # There was a period in which it did not output MeanSignal
    # Use MedianSignal instead...
    if (!defined($g_mean_col) &&
        defined($agilent_header_hash{'gMedianSignal'}))
    {
        printf STDERR "WARNING -- %s missing gMeanSignal, using gMedianSignal\n",
            $filename;
        $g_mean_col = $agilent_header_hash{'gMedianSignal'};
    }
    if (!defined($r_mean_col) &&
        defined($agilent_header_hash{'rMedianSignal'}))
    {
        printf STDERR "WARNING -- %s missing rMeanSignal, using rMedianSignal\n",
            $filename;
        $r_mean_col = $agilent_header_hash{'rMedianSignal'};
    }

    $has_cy5_flag = 0;
    $has_cy3_flag = 0;

    if (defined($g_mean_col) && defined($g_bg_mean_col))
    {
        $has_cy3_flag = 1;
    }
    if (defined($r_mean_col) && defined($r_bg_mean_col))
    {
        $has_cy5_flag = 1;
    }
    
    if ($has_cy3_flag == 0 && $has_cy5_flag == 0)
    {
        die "ERROR -- missing both Cy3 and Cy5 channel data\n";
    }
    
    # Read in data
    while(defined($line=<INFILE>))
    {
        # skip blank lines
        if (!($line =~ /[A-Za-z0-9]/))
        {
            next;
        }

        $line =~ s/[\r\n]+//g;
        $line =~ s/\"//g;

        @array = split /\t/, $line;
        
        for ($i = 0; $i < @array; $i++)
        {
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
            $array[$i] =~ s/\s+/ /g;
            
            # deal with missing data and ArrayExpress-inserted nulls
            if (!($array[$i] =~ /\S/) ||
                $array[$i] =~ /^null$/i)
            {
                $array[$i] = '';
            }
        }

        $feature_num  = $array[$feature_num_col];
        $probe_name   = $array[$probe_name_col];
        $control_type = $array[$control_type_col];

        # skip Agilent control probes
        if ($control_type =~ /[A-Za-z0-9]/)
        {
            # GSE77380 uses 0/1 for data/control
            # don't skip if it is a number and that number is zero
            if (!(is_number($control_type) && $control_type == 0))
            {
                next;
            }
        }

        $read_a_file_flag = 1;

        if ($has_cy3_flag)
        {
            $g_mean      = $array[$g_mean_col];
            $g_bg_mean   = $array[$g_bg_mean_col];
        }

        if ($has_cy5_flag)
        {
            $r_mean      = $array[$r_mean_col];
            $r_bg_mean   = $array[$r_bg_mean_col];
        }
        
        $probeid = sprintf "%05d:%s", $feature_num, $probe_name;
        $probeid_hash{$probeid} = 1;
        
        $bad_spot_flag = 0;
        
        if ($has_cy3_flag)
        {
            $signal = $g_mean;
            if ($subtract_local_background_flag)
            {
                $signal = $g_mean - $g_bg_mean;
            }

#            if ($signal <= $min_signal_cutoff)
#            {
#                $signal = 0;
#            }
            
            $signal_hash{$cy3_name}{$probeid} = $signal;
            $temp_signal_hash{$cy3_name} = $signal;

            # there is something seriously wrong with this spot, so skip it
            if ($signal <= -100)
            {
                $bad_spot_flag = 1;
            }
        }

        if ($has_cy5_flag)
        {
            $signal = $r_mean;
            if ($subtract_local_background_flag)
            {
                $signal = $r_mean - $r_bg_mean;
            }

#            if ($signal <= $min_signal_cutoff)
#            {
#                $signal = 0;
#            }
            
            $signal_hash{$cy5_name}{$probeid} = $signal;
            $temp_signal_hash{$cy5_name} = $signal;

            # there is something seriously wrong with this spot, so skip it
            if ($signal <= -100)
            {
                $bad_spot_flag = 1;
            }
        }
        
        if ($bad_spot_flag)
        {
            # effectively zero it out, for later
            $signal = -9E99;

            if ($has_cy3_flag && defined($cy3_name))
            {
                $signal_hash{$cy3_name}{$probeid} = $signal;

                printf STDERR "WARNING -- bad spot\t%s\t%s\t%g\n",
                       $cy3_name, $probeid, $temp_signal_hash{$cy3_name};
            }
            if ($has_cy5_flag && defined($cy5_name))
            {
                $signal_hash{$cy5_name}{$probeid} = $signal;

                printf STDERR "WARNING -- bad spot\t%s\t%s\t%g\n",
                       $cy5_name, $probeid, $temp_signal_hash{$cy5_name};
            }
        }
        # update "good" spot minimums
        else
        {
            if ($has_cy3_flag)
            {
                if (!defined($min_signal_hash{$cy3_name}))
                {
                    $min_signal_hash{$cy3_name} =
                        $signal_hash{$cy3_name}{$probeid};
                }
                if ($signal_hash{$cy3_name}{$probeid} <
                    $min_signal_hash{$cy3_name})
                {
                    $min_signal_hash{$cy3_name} =
                        $signal_hash{$cy3_name}{$probeid};
                }
            }

            if ($has_cy5_flag)
            {
                if (!defined($min_signal_hash{$cy5_name}))
                {
                    $min_signal_hash{$cy5_name} =
                        $signal_hash{$cy5_name}{$probeid};
                }
                if ($signal_hash{$cy5_name}{$probeid} <
                    $min_signal_hash{$cy5_name})
                {
                    $min_signal_hash{$cy5_name} =
                        $signal_hash{$cy5_name}{$probeid};
                }
            }
        }
    }
   
    close INFILE;
}

@sample_array  = sort keys %signal_hash;
@probeid_array = sort keys %probeid_hash;


# DEBUG code
if (0)
{
  foreach $sample (@sample_array)
  {
    foreach $probeid (@probeid_array)
    {
        $signal = $signal_hash{$sample}{$probeid};
        
        if ($signal < 0)
        {
            printf "%g\n", $signal;
        }
    }
  }
}


# exit with example usage
if ($read_a_file_flag == 0)
{
    printf "No valid input files found\n";
    printf "Usage: agilent_to_spreadsheet.pl agilent_scans*.txt > unnorm_spreadsheet.txt\n";
    die;
}


# print header

printf "ProbeID";
foreach $sample (@sample_array)
{
    printf "\t%s", $sample;
}
printf "\n";



foreach $probeid (@probeid_array)
{
    # strip internal feature number identifier
    $probe_name = $probeid;
    $probe_name =~ s/^\d+://;

    printf "%s", $probe_name;
    
    foreach $sample (@sample_array)
    {
        $min_signal = $min_signal_hash{$sample};
        $signal     = ($signal_hash{$sample}{$probeid} - $min_signal) + 1;

        if ($signal <= $min_signal_cutoff)
        {
            $signal = 0;
        }
    
        printf "\t%g", $signal;
    }
    printf "\n";
}
