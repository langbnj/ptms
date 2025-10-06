#!/bin/sh

if [[ $@ == "" ]]
then
	# echo
	# echo "--------------------------------------------------------------------------------";
	# echo "> No command given, showing queue instead:";
	# echo "--------------------------------------------------------------------------------";
	echo
	/home/blang1/scripts/showqueue.pl
	echo
	exit
fi

if [ ! -z $1 ]
then
	# # Absolute path to this script. /home/user/bin/foo.sh
	# DIR=$PWD
	# echo "DIR $DIR"
	# # Absolute path this script is in. /home/user/bin
	# DIR=`dirname $DIR`
	# echo "DIR $DIR"
	# # Strip /home/blang1/
	# DIR=`echo $DIR | sed 's/^\/home\/lang\/home\/*//'`
	# echo "DIR $DIR"
	# # Get Name
	# NAME=`echo $@`
	# echo "NAME $NAME"
	# # Combine Path and Name
	# NAME="${DIR}/${NAME}"
	# echo "NAME $NAME"
	# # Strip non-alphanumerics
	# NAME=`echo $NAME | sed 's/[^a-zA-Z0-9]/_/g'`
	# echo "NAME $NAME"
	# Absolute path to this script. /home/user/bin/foo.sh
	ORIGNAME="$PWD/`basename $1`"
	# echo "ORIGNAME $ORIGNAME"
	# Strip /home/lang/home/
	ARGS=$@
	shift
	NAME=$ORIGNAME
	NAME=`echo $NAME | sed 's/^\/home\/blang1\/*//'`
	NAME=`echo $NAME $@`
	MYNAME=$NAME
	# echo "NAME $NAME"
	# Strip non-alphanumerics
	NAME=`echo $NAME | sed 's/[^a-zA-Z0-9]/_/g'`
	# echo "NAME $NAME"

	# Shorten to 200ish
	# NAMELEN=${#NAME}
	NAME=${NAME:0:200}
fi

if [[ $ORIGNAME == "" && $ARGS != "" ]]
then
	echo
	echo "--------------------------------------------------------------------------------";
	echo "> Error: '$ARGS' doesn't exist, cancelling";
	echo "--------------------------------------------------------------------------------";
	echo
	exit
fi

if [[ $NAME == "" && $ARGS != "" ]]
then
	echo
	echo "--------------------------------------------------------------------------------";
	echo "> Error: '$ORIGNAME' parsed to '$NAME', cancelling";
	echo "--------------------------------------------------------------------------------";
	echo
	exit
fi

echo "#!/bin/sh" > log-command-$NAME.txt
echo "export TERM=xterm-256color" >> log-command-$NAME.txt
echo "export START=\`date +%s\`" >> log-command-$NAME.txt
# echo "export TERM=xterm-color" >> log-command-$NAME.txt
echo "echo" >> log-command-$NAME.txt
echo "echo \"================================================================================\"" >> log-command-$NAME.txt
echo "echo \"> $ARGS (\`hostname\`)\"" >> log-command-$NAME.txt
echo "echo \"================================================================================\"" >> log-command-$NAME.txt
echo "echo" >> log-command-$NAME.txt
echo "source /home/blang1/.bash_profile" >> log-command-$NAME.txt
echo "echo" >> log-command-$NAME.txt
echo "$ARGS" >> log-command-$NAME.txt
echo "export STOP=\`date +%s\`" >> log-command-$NAME.txt
echo "export ELAPSED=\`expr \$STOP - \$START\`" >> log-command-$NAME.txt
echo "export DAYS=\`expr \$ELAPSED / 86400\`" >> log-command-$NAME.txt
echo "export HOURS=\`expr \$ELAPSED % 86400 / 3600\`" >> log-command-$NAME.txt
echo "export MINS=\`expr \$ELAPSED % 86400 % 3600 / 60\`" >> log-command-$NAME.txt
echo "export SECS=\`expr \$ELAPSED % 86400 % 3600 % 60\`" >> log-command-$NAME.txt
echo "echo" >> log-command-$NAME.txt
echo "echo \"================================================================================\"" >> log-command-$NAME.txt
echo "echo \"> Done! (Time elapsed: \$DAYS days \$HOURS h \$MINS min \$SECS sec)\"" >> log-command-$NAME.txt
echo "echo \"================================================================================\"" >> log-command-$NAME.txt

rm -f log-errors-$NAME.txt
rm -f log-output-$NAME.txt
echo 
echo
echo "--------------------------------------------------------------------------------";
echo "> Submitting job '$NAME'";
echo "--------------------------------------------------------------------------------";
echo
bsub -P idr -J $NAME -o /dev/null -L /bin/bash -q standard -n 1 -R "rusage[mem=15G]" "bash log-command-$NAME.txt > log-output-$NAME.txt 2> log-errors-$NAME.txt"
echo
