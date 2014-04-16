multiwiggle
===========

Tool to create and query files containing multiple wiggle plots.

It accepts WIG and BEDGRAPH file formats.

Example:

mwiggle create mytest.mw *.bedgraph -t 9060 -a NCBI37 -d "My test tracks"

run mwiggle in the command line for help


to install
==========
make

make test

add to the path

export PATH=`pwd`/bin:$PATH


add libs to the path

export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH

Notes on integration with Ensembl Genomes
=================================
MW support has been added  from release/eg/22

eg-web-common/modules/Bio/EnsEMBL/GlyphSet/mw.pm
eg-web-common/modules/EnsEMBL/Web/ImageConfig.pm


add path to the mwiggle installation in SiteDefs.pm, e.g
   $SiteDefs::MWIGGLE_DIR= '/usr/local/mwiggle';

configure MW tracks in species ini files, e.g

[ENSEMBL_INTERNAL_MW_SOURCES]

PhyloCSF = transcript


[PhyloCSF]

source_name = PhyloCSF

source_type = MW

display = tiling

description = PhyloCSF: a comparative genomics method to distinguish protein coding and non-coding regions ( Study <a 
href='http://www.ncbi.nlm.nih.gov/pubmed/21685081'>21685081</a> )

source_url = /nfs/public/rw/ensembl/data/vectorbase/PhyloCSF.mw

score_colour = coral

strand = b

track_order = [1 3 5 0 2 4 6]

track_colour = [burlywood1 burlywood1 burlywood1 burlywood1 burlywood1 burlywood1 chocolate3]




