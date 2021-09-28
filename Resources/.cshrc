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

# notes:
# n1: komodo as of june 2018 (Joakim Hove - Yammer 15/6-18).
#     i read this could be replaced by WRA36b below, but it turns out
#     that when using ssh i miss some important elements in PATH
#     (like the ones needed for ert etc.)
# n2: as of sep 2021, i started using virtual environmens (venv).
#     i chose to put them in a public place so other may use them.
#     i also don't use the --system-site-packages option for the
#     WRA36b. it should also include on /prog/res/komodo/2021.08-py36/enable.csh,
#     but see n1

source /prog/res/komodo/stable-py3/enable.csh                  # n1
#source /project/res/komodo/testing/enable.csh                  # hot from the press..

source /project/RCP/active/venv/agy/WRA36b/bin/activate.csh    # n2

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
