#
if [ $# -lt 1 ] ; then
    echo 'Usage: configure-run exp_name in /expscr directory'
    echo "`date` $0: Wrong usage." >>ERROR.out
    exit 8
fi
#
CURR_DIR=`pwd`
LIBS_DIR=`grep '^LIBS_DIR=' $SRCS_DIR/configure-src | \
                        sed 's/ *$//g' | cut -d'=' -f2`
machine=`grep '^#define ' $LIBS_DIR/machine.h | cut -d' ' -f2 | tr '[A-Z]' '[a-z]'`
march=`grep '^MARCH=' $LIBS_DIR/configure-lib | cut -d'=' -f2`
echo "set variables from $LIBS_DIR/opt_libs/options-$machine-$march"
        cp $LIBS_DIR/opt_libs/options-$machine-$march options
chmod a+x options
. ./options
CPP=`eval echo $CPP`
$CPP -P -I$SRCS_DIR $SRCS_DIR/def/cvar.F > cvar.i
sed '/^ *$/d;s/\//#/;s/\!/#/' cvar.i >cvar.env
rm -f cvar.i
. ./cvar.env
