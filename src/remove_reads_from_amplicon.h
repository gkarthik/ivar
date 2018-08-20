#include<iostream>
#include<stdint.h>

#ifndef removereads_from_amplicon
#define removereads_from_amplicon

int rmv_reads_from_amplicon(std::string bam, std::string region_, std::string bam_out, uint16_t* amplicon, int amp_n);

#endif
