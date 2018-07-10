#!/bin/bash
script_dir=path_to_script
site_file=dm3_DpnII.txt

mnd=$1
o=$(basename $mnd -mnd.txt)

cat $mnd | parallel --pipe -k -j 8 perl ${script_dir}/contact_fragment_get.pl -s $site_file >$o.cf
cat $o.cf | parallel --pipe -k -j 8 perl ${script_dir}/calcContactFragGC.pl >$o.gc
