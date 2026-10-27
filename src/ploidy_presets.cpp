/// \file ploidy_presets.cpp
///
/// The text of each built-in ploidy table. Each table's header cites where its contig lengths and
/// PAR coordinates come from.

#include "ploidy_presets.hpp"

namespace vg {

const map<string, string>& ploidy_presets() {
    static const map<string, string> presets {
        {"chm13v2.XX", R"ploidy(# Ploidy of a female (XX) sample on T2T-CHM13v2.0: chrX diploid, so it has no line,
# chrY absent (ploidy 0, no calls), and chrM haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# Lengths: the CHM13v2.0 assembly (GCA_009914755.4); chrX also matches the CHM13#0#chrX path of
# the HPRC graph.
chrY	0	62460029	0
chrM	0	16569	1
)ploidy"},
        {"chm13v2.XY", R"ploidy(# Ploidy of a male (XY) sample on T2T-CHM13v2.0: chrX haploid outside its two
# pseudoautosomal regions (PARs), which are diploid, and chrY and chrM haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# chrY is haploid along its whole length, PARs included: a male sample's second PAR copy is
# counted on chrX, where the reads usually go because analysis sets mask chrY's PARs.
#
# Lengths: the CHM13v2.0 assembly (GCA_009914755.4); chrX also matches the CHM13#0#chrX path of
# the HPRC graph. The chrX PARs are as the T2T-Q100 HG002 truth set draws them: it is diploid in
# exactly these two blocks. T2T's own chm13v2.0_PAR.bed puts them at 0-2394410 and
# 153925834-154259566, which differs by 40 bp at PAR1's end and 169 bp at PAR2's start.
chrX	2394370	153926003	1
chrY	0	62460029	1
chrM	0	16569	1
)ploidy"},
        {"grch37.b37.XX", R"ploidy(# Ploidy of a female (XX) sample on GRCh37 with b37/hs37d5 names: X diploid, so it has no line,
# Y absent (ploidy 0, no calls), and MT haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# Lengths: the GRCh37 assembly report (GCA_000001405.1); MT is the 16,569 bp rCRS.
Y	0	59373566	0
MT	0	16569	1
)ploidy"},
        {"grch37.b37.XY", R"ploidy(# Ploidy of a male (XY) sample on GRCh37 with b37/hs37d5 names: X haploid outside its two
# pseudoautosomal regions (PARs), which are diploid, and Y and MT haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# Y is haploid along its whole length, PARs included: a male sample's second PAR copy is
# counted on X, where the reads usually go because analysis sets mask Y's PARs.
#
# Lengths: the GRCh37 assembly report (GCA_000001405.1); MT is the 16,569 bp rCRS. PARs: GRC pseudoautosomal regions, https://www.ncbi.nlm.nih.gov/grc/human,
# X:60,001-2,699,520 and X:154,931,044-155,260,560 (1-based, inclusive).
X	0	60000	1
X	2699520	154931043	1
X	155260560	155270560	1
Y	0	59373566	1
MT	0	16569	1
)ploidy"},
        {"grch37.ucsc.XX", R"ploidy(# Ploidy of a female (XX) sample on GRCh37 with UCSC hg19 names: chrX diploid, so it has no line,
# chrY absent (ploidy 0, no calls), and chrM haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# Lengths: the GRCh37 assembly report (GCA_000001405.1), and UCSC hg19.chrom.sizes for chrM,
# which in hg19 is the 16,571 bp NC_001807 rather than the 16,569 bp rCRS.
chrY	0	59373566	0
chrM	0	16571	1
)ploidy"},
        {"grch37.ucsc.XY", R"ploidy(# Ploidy of a male (XY) sample on GRCh37 with UCSC hg19 names: chrX haploid outside its two
# pseudoautosomal regions (PARs), which are diploid, and chrY and chrM haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# chrY is haploid along its whole length, PARs included: a male sample's second PAR copy is
# counted on chrX, where the reads usually go because analysis sets mask chrY's PARs.
#
# Lengths: the GRCh37 assembly report (GCA_000001405.1), and UCSC hg19.chrom.sizes for chrM,
# which in hg19 is the 16,571 bp NC_001807 rather than the 16,569 bp rCRS. PARs: GRC pseudoautosomal regions, https://www.ncbi.nlm.nih.gov/grc/human,
# chrX:60,001-2,699,520 and chrX:154,931,044-155,260,560 (1-based, inclusive).
chrX	0	60000	1
chrX	2699520	154931043	1
chrX	155260560	155270560	1
chrY	0	59373566	1
chrM	0	16571	1
)ploidy"},
        {"grch38.XX", R"ploidy(# Ploidy of a female (XX) sample on GRCh38: chrX diploid, so it has no line,
# chrY absent (ploidy 0, no calls), and chrM haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# Lengths: the GRCh38 assembly report (GCA_000001405.15).
chrY	0	57227415	0
chrM	0	16569	1
)ploidy"},
        {"grch38.XY", R"ploidy(# Ploidy of a male (XY) sample on GRCh38: chrX haploid outside its two
# pseudoautosomal regions (PARs), which are diploid, and chrY and chrM haploid.
#
# CHROM START END PLOIDY, 0-based half-open, with CHROM spelled as in the output VCF. A position no
# line covers is diploid. Ploidy 0 means the sample has no copy of the region, so vg makes no calls
# there.
#
# chrY is haploid along its whole length, PARs included: a male sample's second PAR copy is
# counted on chrX, where the reads usually go because analysis sets mask chrY's PARs.
#
# Lengths: the GRCh38 assembly report (GCA_000001405.15). PARs: GRC pseudoautosomal regions, https://www.ncbi.nlm.nih.gov/grc/human,
# chrX:10,001-2,781,479 and chrX:155,701,383-156,030,895 (1-based, inclusive).
chrX	0	10000	1
chrX	2781479	155701382	1
chrX	156030895	156040895	1
chrY	0	57227415	1
chrM	0	16569	1
)ploidy"}
    };
    return presets;
}

}
