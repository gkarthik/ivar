#!/bin/bash

awk '
    FNR==1 && FNR == NR{
	FILES[1] = FILENAME
    }
    FNR==1 && FNR!=NR {
	FILES[length(FILES) +1] = FILENAME
    }
    {
	if(arr[$1,$2,$3,$4]){
		arr[$1,$2,$3,$4]++
	} else {
		arr[$1,$2,$3,$4] = 1
	}
	line[FILENAME,$1,$2,$3,$4] = $0
    }
    END {
        num_files = ARGC -1
	printf "REGION\tPOS\tREF\tALT\t"
	ORS="\t"
	for(k in FILES){
	      print "REF_DP_" FILES[k] "\tREF_RV_" FILES[k] "\tREF_QUAL_" FILES[k] "\tALT_DP_" FILES[k] "\tALT_RV_" FILES[k] "\tALT_QUAL_" FILES[k] "\tALT_FREQ_" FILES[k] "\tTOTAL_DP_" FILES[k] "\tPVAL_" FILES[k] "\tPASS_" FILES[k]
	}
	ORS="\n"
	printf "\n"
        for ( key in arr ) {
            if ( arr[key] < num_files ) { continue }
	    split(line[FILES[1],key], line_arr, "\t")
	    if (line_arr[ 1 ] == "REGION") { continue }
	    printf "%s\t", line[FILES[1],key]
	    for(f = 2; f <= length(FILES); f++){
		# printf "KEY: %s\nFILE: %s\nLINE: %s\n", key, FILES[f], line[FILES[f],key]
		split(line[FILES[f],key], line_arr, "\t")
		for ( i = 5; i <= length( line_arr ); i++ ) {
                    printf "%s\t", line_arr[ i ]
		}
	    }
	    printf "\n"
        }
    }
' "$@"
