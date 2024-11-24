#!/usr/bin/env bash

if [ -z "$1" ]; then
    echo "Error: No vendor directory supplied!"
    exit 1
fi

VENDOR_DIR=$1

if [ "$2" ]; then
    VERSION=$2
    echo "Using g++ version $VERSION"
    export CC=gcc-$VERSION
    export CXX=g++-$VERSION
    export CPP=g++-$VERSION
fi

if [ ! -d "$VENDOR_DIR" ]; then
    echo "Directory $VENDOR_DIR does not exists, please create it"
    exit 1
fi

NUM_PROCESSORS=8

#SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

pushd "$VENDOR_DIR" || exit
echo "------------------------------> ENTER PBSIM3 INSTALL <------------------------------"

if [ ! -d "pbsim3/" ]; then

    git clone https://github.com/yukiteruono/pbsim3.git
    pushd "pbsim3/" || exit
    # TODO: hardcode checkout commit
    git checkout f3e5f1cf3d0e8346b5e4598ac238b2b570b223e8
    ./configure
    make -j

    popd || exit # "pbsim3"
    echo "Finished installing PbSim3 simulator"

else
    echo "$VENDOR_DIR/pbsim3 already exists, skipping install!"
fi
echo "------------------------------> EXIT PBSIM3 INSTALL <------------------------------"

echo "------------------------------> ENTER LJA INSTALL <------------------------------"

if [ ! -d "LJA" ]; then

    git clone https://github.com/AntonBankevich/LJA.git
    pushd "LJA/" || exit
    git fetch
    git checkout -t origin/experimental_ml
    # git checkout c8cdeaf629c78adcde0603f3e4fa4241eedf9e2e
    #git apply "$SCRIPT_DIR/lja_eval.patch"
    cmake .
    make -j $NUM_PROCESSORS lja
    make -j $NUM_PROCESSORS mlgraph
    make -j $NUM_PROCESSORS jumboDBG
    make -j $NUM_PROCESSORS align_and_print
    popd || exit # "LJA"
    echo "Finished installing La Jolla assembler"

else
    echo "$VENDOR_DIR/LJA already exists, skipping install!"
fi
echo "------------------------------> EXIT LJA INSTALL <------------------------------"

popd || exit # $VENDOR_DIR

exit 0
