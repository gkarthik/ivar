/*! \file */

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "remove_reads_from_amplicon.h"
#include "call_consensus_pileup.h"
#include "call_variants.h"
#include "trim_primer_quality.h"
#include "get_masked_amplicons.h"
#include "suffix_tree.h"

const std::string VERSION = "1.0";

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
    "Usage:	ivar [command <trim|callvariants|filtervariants|consensus|createbed|getmasked|removereads|version|help>]\n"
    "\n"
    "        Command       Description\n"
    "           trim       Trim reads in aligned BAM file\n"
    "       variants       Call variants from aligned BAM file\n"
    " filtervariants       Filter variants across replicates\n"
    "      consensus       Call consensus from aligned BAM file\n"
    "      getmasked       Get amplicons with primer mismatches\n"
    "    removereads       Remove reads from trimmed BAM file\n"
    "        version       Show version information\n"
    "    trimadapter       (EXPERIMENTAL) Trim adapter sequences from reads\n"
    "\n"
    "To view detailed usage for each command type `ivar <command>` \n";
}

void print_trim_usage(){
  std::cout <<
    "Usage: ivar trim -i <input.bam> -b <primers.bed> -p <prefix> [-m <min-length>] [-q <min-quality>] [-s <sliding-window-width>]\n\n"
    "Input Options    Description\n"
    "           -i    (Required) Indexed aligned bam file to trim primers and quality\n"
    "           -b    (Required) BED file with primer sequences and positions\n"
    "           -m    Minimum length of read to retain after trimming (Default: 30)\n"
    "           -q    Minimum quality threshold for sliding window to pass (Default: 20)\n"
    "           -s    Width of sliding window (Default: 4)\n\n"
    "Output Options   Description\n"
    "           -p    (Required) Prefix for the output BAM file\n";
}

void print_variants_usage(){
  std::cout <<
      "Usage: samtools mpileup -A -d 300000 --reference <reference-fasta> -B -Q 0 -F 0 <input.bam> | ivar variants -p <prefix> [-q <min-quality>] [-t <min-frequency-threshold>]\n\n"
    "Note : samtools mpileup output must be piped into ivar variants\n\n"
    "Input Options    Description\n"
    "           -q    Minimum quality score threshold to count base (Default: 20)\n"
    "           -t    Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)\n\n"
    "Output Options   Description\n"
    "           -p    (Required) Prefix for the output tsv variant file\n\n";
}

void print_filtervariants_usage(){
  std::cout <<
    "Usage: ivar filtervariants -p <prefix> replicate-one.tsv replicate-two.tsv ... \n\n"
    "Input: Variant tsv files for each replicate\n\n"
    "Output Options   Description\n"
    "           -p    (Required) Prefix for the output filtered tsv file\n";
}

void print_consensus_usage(){
  std::cout <<
    "Usage: samtools mpileup -A -d 300000 -Q 0 -F 0 <input.bam> | ivar consensus -p <prefix> \n\n"
    "Note : samtools mpileup output must be piped into `ivar consensus`\n\n"
    "Input Options    Description\n"
    "           -q    Minimum quality score threshold to count base (Default: 20)\n"
    "           -t    Minimum frequency threshold(0 - 1) to call consensus. (Default: 0)\n"
    "                 Frequently used thresholds | Description\n"
    "                 ---------------------------|------------\n"
    "                                          0 | Majority or most common base\n"
    "                                        0.2 | Bases that make up atleast 20% of the depth at a position\n"
    "                                        0.5 | Strict or bases that make up atleast 50% of the depth at a position\n"
    "                                        0.9 | Strict or bases that make up atleast 90% of the depth at a position\n"
    "                                          1 | Identical or bases that make up 100% of the depth at a position. Will have highest ambiguities\n"
    "Output Options   Description\n"
    "           -p    (Required) Prefix for the output fasta file and quality file\n";
}

void print_removereads_usage(){
  std::cout <<
    "Usage: ivar removereads -i <input.trimmed.bam> -p <prefix> primer-index-1 primer-index-2 primer-index-3 primer-index-4 ...  \n\n"
    "Input Options    Description\n"
    "           -i    (Required) Input BAM file  trimmed with ‘ivar trim’. Must be sorted and indexed, which can be done using sort_index_bam.sh\n"
    "           -t    (Required) Text file with primer indices separated by spaces\n"
    "Output Options   Description\n"
    "           -p    (Required) Prefix for the output filtered BAM file\n";
}

void print_getmasked_usage(){
  std::cout <<
    "Usage: ivar getmasked -i <input-filtered.tsv> -b <primers.bed> -p <prefix>\n\n"
    "Input Options    Description\n"
    "           -i    (Required) Input filtered variants tsv generated from `ivar filtervariants`\n"
    "           -b    (Required) BED file with primer sequences and positions\n"
    "Output Options   Description\n"
    "           -p    (Required) Prefix for the output text file\n";
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

void print_version_info(){
  std::cout << "iVar version " << VERSION << std::endl <<
    "\nPlease raise issues and bug reports at https://github.com/andersen-lab/ivar/\n\n";
}

static const char *trim_opt_str = "i:b:p:m::q::s::h?";
static const char *variants_opt_str = "p:t::q::h?";
static const char *consensus_opt_str = "p:q::t::h?";
static const char *removereads_opt_str = "i:p:t:h?";
static const char *filtervariants_opt_str = "p:h?";
static const char *getmasked_opt_str = "i:b:p:h?";
static const char *trimadapter_opt_str = "1:2::p:a:h?";

std::string get_filename_without_extension(std::string f, std::string ext){
  if(f.substr(f.length()-ext.length(), ext.length()).compare(ext) == 0){
    return f.substr(0,f.length()-ext.length());
  }
  return f;
}

/*!
 Main Function

 This is where the command line arguments into iVar are parsed.
 iVar first parses the first argument and depending on the command it either returns the version or uses GNU getopt to parse the remaining arguments.
 */

int main(int argc, char* argv[]){
  if(argc == 1){
    print_usage();
    return -1;
  }
  std::string cmd(argv[1]);
  if(cmd.compare("-v") == 0){
    print_version_info();
    return 0;
  }
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
    g_args.prefix = get_filename_without_extension(g_args.prefix,".bam");
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
    g_args.prefix = get_filename_without_extension(g_args.prefix,".tsv");
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
    g_args.prefix = get_filename_without_extension(g_args.prefix,".fa");
    g_args.prefix = get_filename_without_extension(g_args.prefix,".fasta");
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
      case 't':
	g_args.bed = optarg;
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
    if(g_args.bam.empty() || g_args.prefix.empty() || g_args.bed.empty()){
      print_removereads_usage();
      return -1;
    }
    uint16_t amp[150], i = 0;	// Max primer indices set to 150.
    std::string s;
    std::ifstream fin(g_args.bed);
    while(getline(fin, s, ' ' ) ){
      amp[i] = stoi(s);
      i++;
    }
    g_args.prefix = get_filename_without_extension(g_args.prefix,".bam");
    res = rmv_reads_from_amplicon(g_args.bam, g_args.region, g_args.prefix, amp, i);
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
    g_args.prefix = get_filename_without_extension(g_args.prefix,".tsv");
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
      case 'p':
	g_args.prefix = optarg;
	break;
      case 'h':
      case '?':
	print_getmasked_usage();
	return 0;
	break;
      }
      opt = getopt( argc, argv, getmasked_opt_str);
    }
    if(g_args.bed.empty() || g_args.bam.empty() || g_args.prefix.empty()){
      print_getmasked_usage();
      return -1;
    }
    g_args.prefix = get_filename_without_extension(g_args.prefix,".txt");
    res = get_primers_with_mismatches(g_args.bed, g_args.bam, g_args.prefix);
  } else if (cmd.compare("trimadapter") == 0){
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
  } else if(cmd.compare("version") == 0){
    print_version_info();
  } else {
    print_usage();
  }
  return res;
}
