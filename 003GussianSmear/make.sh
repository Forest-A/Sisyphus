#!/bin/bash

compile() {
    code="$1"
    output="$2"

    echo "Compiling $code into $output"

    # Check if root-config is available
    if ! command -v root-config &> /dev/null; then
        echo "Error: root-config not found. Make sure ROOT is installed and the environment is set up correctly."
        exit 1
    fi

    g++ -g -O3 -Wall -Werror -std=c++17 -fPIC -o "$output" "$code" $(root-config --cflags --libs) -lTree -lMinuit -lASImage -L/home/motoko/TKI_Comparison/style -lstyle -Wl,-rpath=/home/motoko/TKI_Comparison/style



    if [ -e "$output" ]; then
        echo "Compilation successful: $output"
    else
        echo "Compilation failed"
        exit 1
    fi
}

# Check for the correct number of arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 source_file output_executable"
    exit 1
fi

# Call the compile function with the provided arguments
compile "$1" "$2"


