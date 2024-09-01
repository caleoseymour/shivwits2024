## Downloads a bunch of information from KEGG using the REST api.

rm module
rm ko
rm brite


wget http://rest.kegg.jp/list/module
wget http://rest.kegg.jp/list/ko
wget http://rest.kegg.jp/list/brite

## Download KEGG module definitions and parse them into a file.
touch module_raw.txt
cut -f1 module | while read -r m1;
do
    read -r m2
    read -r m3
    read -r m4
    read -r m5
    read -r m6
    read -r m7
    read -r m8
    read -r m9
    read -r m10
    
    mlist=$(echo ${m1}+${m2}+${m3}+${m4}+${m5}+${m6}+${m7}+${m8}+${m9}+${m10} | sed 's/\+*$//g')
    wget -qO-  http://rest.kegg.jp/get/${mlist} >> module_raw.txt
done

grep "DEFINITION" module_raw.txt | sed "s/DEFINITION[[:space:]]*//g" > module_definitions.txt

paste module module_definitions.txt > module_data.tsv

# ## Header for XML files. We don't want these.
xmlhead='<?xml version="1.0" encoding="UTF-8"?>'

mkdir brite_files
touch brite_raw.txt
cut -f1 brite | while read -r m1;
do
    
    fname=$(echo ${m1} | sed 's/^br://g')
    echo $fname
    wget -qO- http://rest.kegg.jp/get/${m1} >> brite_files/$fname.txt
    
    if [ "$(head -1 brite_files/$fname.txt)" != "$xmlhead" ];
    then
        grep -v "#" brite_files/$fname.txt >> brite_raw.txt
    fi
done

