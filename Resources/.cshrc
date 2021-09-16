########################################################################
#
# Filname : .cshrc
#
# Private .cshrc 
#
# This file should be copied to the users HOME directory 
#
########################################################################
#
set ext=cshrc
set LoginPath=/adm/login
#
# Set personal configuration here.
# Read the file /adm/login/VarHelp.cshrc for information.
#========================================================
# 
#
# DO NOT CHANGE THIS LINE
source $LoginPath/DEFAULT.$ext
#
# Personal variables after this line :
#========================================================
#
#setenv PRINTER grb3u
# 
#if ( $?prompt == 0 || $?VUE != 0 ) exit
#
# Personal aliases after this line :
#========================================================
#

#source /prog/res/komodo/stable-py3/enable.csh   # komodo as of june 2018 (Joakim Hove - Yammer 15/6-18). replaced by WRA36b below
#source /project/res/komodo/testing/enable.csh   # hot from the press..
source /project/RCP/active/venv/agy/WRA36b/bin/activate.csh    # based on /prog/res/komodo/2021.08-py36/enable.csh

# Then add *my* stuff. Keep fingers crossed..
source $HOME/.mycshrc



# OLD STUFF
#source /project/res/SDP_cshrc                  # "geriljavirksomhet under radaren..." Per Arne Slotte. should still be active according to Hove 15/6-18
# It is *very* confusing with these 'offical' files that needs to be source'd.
# The order of them is *not* irrelevant, it turns out.
# What I found, to make both jupyter and pylab work, was to use the komodo version of python,
# but make sure pythonpath does not contain 'my' jupyter installation.
# see pip_ in .mycshrc
#setenv PYTHON_VERSION "2.7.14"                 # needs this to make pip work...
#source /prog/sdpsoft/env.csh                   # needs this to make pip work... and have also been recommended to do this anyway...

# NO CHANGES AFTER THIS LINE
source $LoginPath/END.$ext
