code=$1

if [ ${code}aa == aa ]
then
    echo help!

    exit
fi

mkdir -p outplot

codein=$code.C

comm="g++ -g -O3 -Wall -Werror -std=c++17 -fPIC -o $code $codein $(root-config --cflags --libs) -lTree -lMinuit -lASImage -L/home/motoko/TKI_Comparison/style -lstyle -Wl,-rpath=/home/motoko/TKI_Comparison/style"

echo $comm

eval $comm

./$code 
