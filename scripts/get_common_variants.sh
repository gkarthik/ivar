#!/bin/bash

awk '
    FNR==1 && FNR == NR{
	FILES[1] = FILENAME
    }
    FNR==1 && FNR!=NR {
	FILES[length(FILES) +1] = FILENAME
    }
    {
	if(arr[$1,$2,$3]){
		arr[$1,$2,$3]++
	} else {
		arr[$1,$2,$3] = 1
	}
	line[FILENAME,$1,$2,$3] = $0
    }
    END {
        num_files = ARGC -1
	printf "POS\tREF\tALT\t"
	for(k in FILES){
	      print "AD_" FILES[k] "\tRAD_" FILES[k] "\tDP_" FILES[k] "\tQUAL_" FILES[k] "\t"
	}
        for ( key in arr ) {
            if ( arr[key] < num_files ) { continue }
	    printf "%s\t", line[FILES[1],key]
	    for(f = 2; f <= length(FILES); f++){
		# printf "KEY: %s\nFILE: %s\nLINE: %s\n", key, FILES[f], line[FILES[f],key]
		split(line[FILES[f],key], line_arr, "\t")
		for ( i = 4; i <= length( line_arr ); i++ ) {
                    printf "%s\t", line_arr[ i ]
		}
	    }
	    printf "\n"
        }
    }
' "$@"
