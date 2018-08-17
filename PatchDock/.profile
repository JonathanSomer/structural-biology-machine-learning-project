#
# This is the default standard .profile provided to sh users.
# They are expected to edit it to meet their own needs.
#
# The commands in this file are executed when an sh user first
# logs in.
#
# $Revision: 1.10 $
#

# Set the interrupt character to Ctrl-c and do clean backspacing.
if [ -t 0 ]
then
	stty intr '^C' echoe 
fi

# Set the TERM environment variable
eval `tset -s -Q`

# Set the default X server.
if [ ${DISPLAY:-setdisplay} = setdisplay ]
then
    if [ ${REMOTEHOST:-islocal} != islocal ]
    then
        DISPLAY=${REMOTEHOST}:0
    else
        DISPLAY=:0
    fi
    export DISPLAY
fi


# List files in columns if standard out is a terminal.
ls()	{ if [ -t ]; then /bin/ls -C $*; else /bin/ls $*; fi }
