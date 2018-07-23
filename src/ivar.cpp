#include <iostream>
#include <unistd.h>

#include "remove_reads_from_amplicon.h"
#include "call_consensus_pileup.h"
#include "call_variants.h"
#include "trim_primer_quality.h"
#include "get_masked_amplicons.h"
#include "suffix_tree.h"

struct args_t {
  std::string bam;		// -i
  std::string bed;		// -b
  std::string prefix;		// -p
  std::string ref;		// -r
  std::string region;		// -R
  uint8_t min_qual;		// -q
  uint8_t sliding_window;	// -s
  double min_threshold;	// -t
  int min_length;		// -m
  std::string f1;		// -1
  std::string f2;		// -2
  std::string adp_path;	// -a
} g_args;

void print_usage(){
  std::cout <<
    "Usage:	ivar [command <trim|callvariants|filtervariants|consensus|createbed|removereads|getmasked>]\n"
    "\n"
    "        Command       Description\n"
    "           trim       Trim reads in aligned bams\n"
    "       variants       Call Variants from aligned bam\n"
    " filtervariants       Filter variants across replicates\n"
    "      consensus       Call consensus from vcf file\n"
    "    removereads       Remove reads from aligned bam\n"
    "      getmasked       Get amplicons with primer mismatches\n"
    "      trimadapter     Trim adapter sequences from reads\n"
    "\n"
    "To view detailed usage for each command type `ivar <command>` \n";
}

void print_trim_usage(){
  std::cout <<
    "Usage: ivar trim [-i input] [-b <bedfile>] [-p prefix]\n\n"
    "Input Options    Description\n"
    "           -i    (Required) Input, aligned bam file to trim primers and quality\n"
    "           -b    (Required) Bed File, Bed file with primer sequences and positions.\n"
    "           -m    Minimum length of read\n"
    "           -q    Minimum quality threshold for sliding window to pass (Default: 20)\n"
    "           -s    Width of sliding window (Default: 4)\n\n"
    "Output Options   Description\n"
    "           -p    (Required) Prefix, prefix for the output file\n";
}

void print_variants_usage(){
  std::cout <<
    "Usage: samtools mpileup -A -d 300000 --reference <reference-fasta> -Q 0 -F 0 -i <input-bam> | ivar variants [-p prefix]\n\n"
    "Note : samtools mpileup output must be piped into ivar variants\n\n"
    "Input Options    Description\n"
    "           -r    (Required) Reference fasta file\n"
    "           -q    Minimum quality threshold to count base (Default: 20)\n"
    "           -t    Minimum Frequency Threshold(0 - 1) to call variants (Default: 0.03)."
    "Output Options   Description\n"
    "           -p    (Required) Prefix, prefix for the output tsv file\n";
}

void print_filtervariants_usage(){
  std::cout <<
    "Usage: ivar filtervariants [-p prefix] replicate-one.tsv replicate-two.csv ... \n\n"
    "Output Options   Description\n"
    "           -p    (Required) Prefix, prefix for the filtered tsv file\n";
}

void print_consensus_usage(){
  std::cout <<
    "Usage: samtools mpileup -A -d 300000 -Q 0 -F 0 [<input-bam>] | ivar consensus [-p prefix] \n\n"
    "Note : samtools mpileup output must be piped into `ivar consensus`\n\n"
    "Input Options    Description\n"
    "           -q    Minimum quality threshold to count base (Default: 20)\n"
    "           -t    Threshold(0 - 100) to call consensus (Default: 0)\n"
    "Output Options   Description\n"
    "           -p    (Required) Prefix, prefix for the output tsv file\n";
}

void print_removereads_usage(){
  std::cout <<
    "Usage: ivar removereads [-i <input-bam>] [-p prefix] primer-index-1 primer-index-2 primer-index-3 primer-index-4 ...  \n\n"
    "Input Options    Description\n"
    "           -i    (Required) Input BAM file run through `ivar trim` command whcih adds the primer number to BAM auxillary data\n"
    "Output Options   Description\n"
    "           -p    (Required) Prefix, prefix for the filtered BAM file\n";
}

void print_getmasked_usage(){
  std::cout <<
    "Usage: ivar getmasked [-i <input-filtered-tsv>] [-b <bed-file>]\n\n"
    "Input Options    Description\n"
    "           -i    (Required) Input filtered variants tsv generated from `ivar filtervariants`\n"
    "           -b    (Required) Bed file with primer indices\n";
}

void print_trimadapter_usage(){
  std::cout <<
    "Usage: ivar trimadapter [-f1 <input-fastq>] [-f2 <input-fastq-2>] [-p prefix] [-a <adapter-fasta-file>]\n\n"
    "Input Options    Description\n"
    "           -1    (Required) Input fastq file\n"
    "           -2    Input fastq file 2 (for pair ended reads)\n"
    "           -a    (Required) Adapter Fasta File\n"
    "           -p    (Required) Prefix of output fastq files\n";
}

static const char *trim_opt_str = "i:b:p:m::q::s::h?";
static const char *variants_opt_str = "p:t::q::h?";
static const char *consensus_opt_str = "p:q::t::h?";
static const char *removereads_opt_str = "i:p:h?";
static const char *filtervariants_opt_str = "p:h?";
static const char *getmasked_opt_str = "i:b:h?";
static const char *trimadapter_opt_str = "1:2::p:a:h?";

int main(int argc, char* argv[]){
  if(argc == 1){
    print_usage();
    return -1;
  }
  std::string cmd(argv[1]);
  int opt = 0, res = 0;
  // Sift arg by 1 for getopt
  argv[1] = argv[0];
  argv++;
  argc--;
  g_args.min_qual = 255;
  g_args.sliding_window = 255;
  g_args.min_threshold = -1;
  g_args.min_length = -1;
  if (cmd.compare("trim") == 0){
    opt = getopt( argc, argv, trim_opt_str);
    while( opt != -1 ) {
      switch( opt ) {
      case 'i':
	g_args.bam = optarg;
	break;
      case 'b':
	g_args.bed = optarg;
	break;
      case 'p':
	g_args.prefix = optarg;
	break;
      case 'm':
	g_args.min_length = atoi(optarg);
	break;
      case 'q':
	g_args.min_qual = atoi(optarg);
	break;
      case 's':
	g_args.sliding_window = atoi(optarg);
	break;
      case 'h':
      case '?':
	print_trim_usage();
	return -1;
	break;
      }
      opt = getopt( argc, argv, trim_opt_str);
    }
    if(g_args.bam.empty() || g_args.bed.empty() || g_args.prefix.empty()){
      print_trim_usage();
      return -1;
    }
    g_args.min_qual = (g_args.min_qual == 255) ? 20 : g_args.min_qual;
    g_args.sliding_window = (g_args.sliding_window == 255) ? 4 : g_args.sliding_window;
    g_args.min_length = (g_args.min_length == -1) ? 30 : g_args.min_length;
    res = trim_bam_qual_primer(g_args.bam, g_args.bed, g_args.prefix, g_args.region, g_args.min_qual, g_args.sliding_window, g_args.min_length);
  } else if (cmd.compare("variants") == 0){
    opt = getopt( argc, argv, variants_opt_str);
    while( opt != -1 ) {
      switch( opt ) {
      case 'p':
	g_args.prefix = optarg;
	break;
      case 't':
	g_args.min_threshold = atof(optarg);
	break;
      case 'q':
	g_args.min_qual = atoi(optarg);
	break;
      case 'h':
      case '?':
	print_variants_usage();
	return 0;
	break;
      }
      opt = getopt( argc, argv, variants_opt_str);
    }
    if(g_args.prefix.empty()){
      print_variants_usage();
      return -1;
    }
    g_args.min_threshold = (g_args.min_threshold == -1 || g_args.min_threshold < 0 || g_args.min_threshold >= 1) ? 0.03: g_args.min_threshold;
    g_args.min_qual = (g_args.min_qual == 255) ? 20 : g_args.min_qual;
    res = call_variants_from_plup(std::cin, g_args.prefix, g_args.min_qual, g_args.min_threshold);
  } else if (cmd.compare("consensus") == 0){
    opt = getopt( argc, argv, consensus_opt_str);
    g_args.min_threshold = 0;
    while( opt != -1 ) {
      switch( opt ) {
      case 't':
	g_args.min_threshold = atof(optarg);
	break;
      case 'p':
	g_args.prefix = optarg;
	break;
      case 'q':
	g_args.min_qual = atoi(optarg);
	break;
      case 'h':
      case '?':
	print_consensus_usage();
	return 0;
	break;
      }
      opt = getopt( argc, argv, consensus_opt_str);
    }
    if(g_args.prefix.empty()){
      print_consensus_usage();
      return -1;
    }
    g_args.min_qual = (g_args.min_qual == 255) ? 20 : g_args.min_qual;
    std::cout <<"Min Quality" << g_args.min_qual << std::endl;
    std::cout << "Threshold: " << g_args.min_threshold << std::endl;
    res = call_consensus_from_plup(std::cin, g_args.prefix, g_args.min_qual, g_args.min_threshold);
  } else if (cmd.compare("removereads") == 0){
    opt = getopt( argc, argv, removereads_opt_str);
    while( opt != -1 ) {
      switch( opt ) {
      case 'i':
	g_args.bam = optarg;
	break;
      case 'p':
	g_args.prefix = optarg;
	break;
      case 'h':
      case '?':
	print_removereads_usage();
	return 0;
	break;
      }
      opt = getopt( argc, argv, removereads_opt_str);
    }
    if(optind >= argc){
      print_removereads_usage();
      return -1;
    }
    uint16_t amp[100];
    for(int i = optind; i<argc;i++){
      amp[i] = atoi(argv[i]);
      std::cout << amp[i];
    }
    if(g_args.bam.empty() || g_args.prefix.empty()){
      print_removereads_usage();
      return -1;
    }
    res = rmv_reads_from_amplicon(g_args.bam, g_args.region, g_args.prefix, amp, argc - optind);
  } else if(cmd.compare("filtervariants") == 0){
    opt = getopt( argc, argv, filtervariants_opt_str);
    while( opt != -1 ) {
      switch( opt ) {
      case 'p':
	g_args.prefix = optarg;
	break;
      case 'h':
      case '?':
	print_filtervariants_usage();
	return 0;
	break;
      }
      opt = getopt( argc, argv, filtervariants_opt_str);
    }
    if(optind >= argc){
      print_filtervariants_usage();
      return -1;
    }
    if(g_args.prefix.empty()){
      print_filtervariants_usage();
      return -1;
    }
    std::string rep = "get_common_variants.sh ";
    for(int i = optind; i<argc;i++){
      rep += argv[i];
      rep += " ";
    }
    rep += " | sort -s -n -k 2 > "+g_args.prefix+".tsv";
    system(rep.c_str());
  } else if(cmd.compare("getmasked") == 0){
    opt = getopt( argc, argv, getmasked_opt_str);
    while( opt != -1 ) {
      switch( opt ) {
      case 'i':
	g_args.bam = optarg;
	break;
      case 'b':
	g_args.bed = optarg;
	break;
      case 'h':
      case '?':
	print_getmasked_usage();
	return 0;
	break;
      }
      opt = getopt( argc, argv, getmasked_opt_str);
    }
    if(g_args.bed.empty() || g_args.bam.empty()){
      print_getmasked_usage();
      return -1;
    }
    res = get_primers_with_mismatches(g_args.bed, g_args.bam);
  } else if(cmd.compare("trimadapter") == 0){
    opt = getopt( argc, argv, trimadapter_opt_str);
    while( opt != -1 ) {
      switch( opt ) {
      case '1':
	g_args.f1 = optarg;
	break;
      case '2':
	g_args.f2 = optarg;
	break;
      case 'p':
	g_args.prefix = optarg;
	break;
      case 'a':
	g_args.adp_path = optarg;
	break;
      case 'h':
      case '?':
	print_trimadapter_usage();
	return 0;
	break;
      }
      opt = getopt( argc, argv, trimadapter_opt_str);
    }
    if(g_args.f1.empty() || g_args.prefix.empty() || g_args.adp_path.empty()){
      print_trimadapter_usage();
      return -1;
    }
    res = trim_adapter(g_args.f1, g_args.f2, g_args.adp_path, g_args.prefix);
  } else {
    print_usage();
  }
  return res;
}
