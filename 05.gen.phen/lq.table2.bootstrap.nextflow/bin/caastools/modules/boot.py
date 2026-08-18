#                      _              _     
#                     | |            | |    
#   ___ __ _  __ _ ___| |_ ___   ___ | |___ 
#  / __/ _` |/ _` / __| __/ _ \ / _ \| / __|
# | (_| (_| | (_| \__ \ || (_) | (_) | \__ \
#  \___\__,_|\__,_|___/\__\___/ \___/|_|___/


'''
A Convergent Amino Acid Substitution identification 
and analysis toolbox

Author:         Fabio Barteri (fabio.barteri@upf.edu)

Contributors:   Alejandro Valenzuela (alejandro.valenzuela@upf.edu)
                Xavier Farré (xfarrer@igtp.cat),
                David de Juan (david.juan@upf.edu).


MODULE NAME: boot.py
DESCRIPTION: bootstrap function
DEPENDENCIES: alimport.py, caas_id.py, pindex.py
CALLED BY: ct

'''


from modules.init_bootstrap import *
from modules.disco import process_position
from modules.caas_id import iscaas
from modules.alimport import *
from modules.hyper import calcpval_random

from os.path import exists
import functools


# FUNCTION fetch_caas():
# fetches caas per each thing


# FUNCTION infer_resampled_group_sizes()
# Infers and validates the nominal FG/BG sizes from resampled traits.

def infer_resampled_group_sizes(resampled_traits):

    fg_traits = set(resampled_traits.trait2fg.keys())
    bg_traits = set(resampled_traits.trait2bg.keys())
    trait_names = fg_traits.union(bg_traits)

    if len(trait_names) == 0:
        raise ValueError("the resampled traits file contains no valid cycles")

    if fg_traits != bg_traits:
        raise ValueError("each resampling cycle must contain both foreground and background species")

    if len(trait_names) != resampled_traits.cycles:
        raise ValueError("resampling cycle names must be unique and every line must define a valid cycle")

    group_sizes = set()

    for trait in trait_names:
        fg_species = resampled_traits.trait2fg[trait]
        bg_species = resampled_traits.trait2bg[trait]

        if len(fg_species) == 0 or len(bg_species) == 0:
            raise ValueError("foreground and background groups cannot be empty")

        if len(fg_species) != len(set(fg_species)) or len(bg_species) != len(set(bg_species)):
            raise ValueError("a species cannot occur more than once in the same resampling group")

        if len(set(fg_species).intersection(set(bg_species))) > 0:
            raise ValueError("foreground and background groups must not overlap within a resampling cycle")

        group_sizes.add((len(fg_species), len(bg_species)))

    if len(group_sizes) != 1:
        formatted_sizes = ", ".join([str(x[0]) + "/" + str(x[1]) for x in sorted(group_sizes)])
        raise ValueError("all resampling cycles must have the same FG/BG sizes; found " + formatted_sizes)

    return list(group_sizes)[0]


# FUNCTION filter_positions_by_hypergeometric_pvalue()
# Excludes alignment positions whose random CAAS p-value exceeds the threshold.

def filter_positions_by_hypergeometric_pvalue(sliced_object, fg_size, bg_size, threshold):

    positions_before_filtering = len(sliced_object.d)
    retained_positions = []

    for position in sliced_object.d:
        pvalue = calcpval_random(position, sliced_object.genename, fg_size, bg_size)

        if pvalue <= threshold:
            retained_positions.append(position)

    sliced_object.d = retained_positions

    return positions_before_filtering, len(retained_positions)

# FUNCTION filter_for_gaps()
# filters a trait for its gaps

def filter_for_gaps(max_bg, max_fg, max_all, gfg, gbg):

    out = True

    all_g = gfg + gbg
    
    if max_all != "NO" and all_g > int(max_all):
        out = False

    elif max_fg != "NO" and gfg > int(max_fg):
        out = False

    elif max_bg != "NO" and gbg > int(max_bg):
        out = False

    return out


# FUNCTION filter_for_missing()
# filters a trait for its gaps

def filter_for_missings(max_m_bg, max_m_fg, max_m_all, mfg, mbg):

    out = True

    all_m = mfg + mbg
    
    if max_m_all != "NO" and all_m > int(max_m_all):
        out = False

    elif max_m_fg != "NO" and mfg > int(max_m_fg):
        out = False

    elif max_m_bg != "NO" and mbg > int(max_m_bg):
        out = False

    return out



def caasboot(processed_position, genename, list_of_traits, maxgaps_fg, maxgaps_bg, maxgaps_all, maxmiss_fg, maxmiss_bg, maxmiss_all, minobserved_fg, minobserved_bg, cycles, admitted_patterns):

    a = set(list_of_traits)
    b = set(processed_position.trait2aas_fg.keys())
    c = set(processed_position.trait2aas_bg.keys())

    valid_traits = list(a.intersection(b).intersection(c))

    def filter_trait(trait, the_processed_position, max_bg_gaps, max_fg_gaps, max_all_gaps, max_bg_miss, max_fg_miss, max_all_miss, min_observed_fg, min_observed_bg):

        the_gfg = the_processed_position.trait2gaps_fg[trait]
        the_gbg = processed_position.trait2gaps_bg[trait]

        the_mfg = the_processed_position.trait2miss_fg[trait]
        the_mbg = processed_position.trait2miss_bg[trait]

        ### Minimum observed species filtering.

        # Observed species are members of the resampled FG/BG group that are
        # present in the alignment and do not carry a gap at this position.
        if len(the_processed_position.trait2ungapped_fg[trait]) < int(min_observed_fg):
            return False
        if len(the_processed_position.trait2ungapped_bg[trait]) < int(min_observed_bg):
            return False
    

        ### GAPS filtering.

        if maxgaps_fg != "NO" and processed_position.trait2gaps_fg[trait] > int(maxgaps_fg):
            return False
        if maxgaps_bg != "NO" and processed_position.trait2gaps_bg[trait] > int(maxgaps_fg):
            return False
        if maxgaps_all != "NO" and processed_position.trait2gaps_fg[trait] + processed_position.trait2gaps_bg[trait] > int(maxgaps_all):
           return False

        ### Missings filtering.
        
        if maxmiss_fg != "NO" and processed_position.trait2miss_fg[trait] > int(maxmiss_fg):
            return False
        if maxmiss_bg != "NO" and processed_position.trait2miss_bg[trait] > int(maxmiss_fg):
            return False
        if maxmiss_all != "NO" and processed_position.trait2miss_fg[trait] + processed_position.trait2miss_bg[trait] > int(maxmiss_all):
            return False

        else:
            aa_tag_fg_list = the_processed_position.trait2aas_fg[trait]
            aa_tag_fg_list.sort()
            aa_tag_fg = "".join(aa_tag_fg_list)

            aa_tag_bg_list = the_processed_position.trait2aas_bg[trait]
            aa_tag_bg_list.sort()
            aa_tag_bg = "".join(aa_tag_bg_list)

            tag = "/".join([aa_tag_fg, aa_tag_bg])

            check = iscaas(tag)

            if check.caas == True and check.pattern in admitted_patterns:
                return True
    
    output_traits = filter(functools.partial(
        filter_trait,
        the_processed_position = processed_position,
        max_bg_gaps = maxgaps_bg,
        max_fg_gaps = maxgaps_fg,
        max_all_gaps = maxgaps_all,

        max_bg_miss = maxmiss_bg,
        max_fg_miss = maxmiss_fg,
        max_all_miss = maxmiss_all,

        min_observed_fg = minobserved_fg,
        min_observed_bg = minobserved_bg,
    ),
                    valid_traits)



    output_traits = list(output_traits)

    # Return the line
    
    position_name = genename + "@" + str(processed_position.position)
    count = str(len(output_traits))

    traitline = ",".join(output_traits)
    empval = str(int(count)/cycles)

    # The hypergeometric p-value describes the alignment position for the
    # nominal FG/BG group sizes. It is shared by all positive bootstrap cycles
    # reported in traitline and is printed immediately after that field.
    pvalue_string = str(processed_position.hypergeometric_pvalue)
    outline = "\t".join([position_name, count, str(cycles), empval, traitline, pvalue_string])

    return outline
        

# FUNCTION disco_bootstrap()
# Launches the bootstrap in several lines. Returns a dictionary gene@position --> pvalue

def boot_on_single_alignment(trait_config_file, resampled_traits, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, min_fg_observed, min_bg_observed, the_admitted_patterns, output_file):


    # Step 3: processes the positions from imported alignment (process_position() from caas_id.py)
    processed_positions = list(map(functools.partial(process_position, multiconfig = resampled_traits, species_in_alignment = sliced_object.species), sliced_object.d))
    the_genename = sliced_object.genename
    print("caastools found", resampled_traits.cycles, "resamplings")

    # Calculate the same positional hypergeometric p-value used by
    # --filter_significant, including when that prefilter is disabled. The
    # value is attached to the processed position so caasboot can serialize it.
    fg_size, bg_size = infer_resampled_group_sizes(resampled_traits)
    for raw_position, processed_position in zip(sliced_object.d, processed_positions):
        processed_position.hypergeometric_pvalue = calcpval_random(
            raw_position,
            the_genename,
            fg_size,
            bg_size,
        )

    # Step 4: extract the raw caas

    output_lines = map(
        functools.partial(
            caasboot,
            list_of_traits = resampled_traits.alltraits,
            genename = the_genename,
            maxgaps_fg = max_fg_gaps,
            maxgaps_bg = max_bg_gaps,
            maxgaps_all = max_overall_gaps,

            maxmiss_fg = max_fg_miss,
            maxmiss_bg = max_bg_miss,
            maxmiss_all = max_overall_miss,

            minobserved_fg = min_fg_observed,
            minobserved_bg = min_bg_observed,

            admitted_patterns = the_admitted_patterns,
            cycles = resampled_traits.cycles) ,processed_positions
    )

    output_lines = list(output_lines)

    ooout = open(output_file, "w")

    for line in output_lines:
        print(line + "\t" + trait_config_file, file=ooout)
    
    ooout.close()

# FUNCTION pval()
# Returns a dictionary with the pvalue

def pval(bootstrap_result):
    with open(bootstrap_result) as h:
        thelist = h.read().splitlines()
    
    d = {}

    for line in thelist:
        try:
            c = line.split("\t")
            d[c[0]] = c[2]
        except:
            pass
    
    return d
