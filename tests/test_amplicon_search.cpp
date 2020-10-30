#include <iostream>
#include <vector>

#include "../src/primer_bed.h"
#include "../src/interval_tree.h"

IntervalTree test_populate_amplicons(std::string pair_info_file, std::string bed_file)
{
    std::vector<primer> primers;
    primers = populate_from_file(bed_file);
    int amplicon_start = -1;
    int amplicon_end = -1;
    IntervalTree tree = IntervalTree();
    std::cout << "Total Number of primers: " << primers.size() << std::endl;
    populate_pair_indices(primers, pair_info_file);
    for (auto & p : primers) {
        if (p.get_strand() == '+')
        {
            if (p.get_pair_indice() != -1){
                amplicon_start = p.get_start();
                amplicon_end = primers[p.get_pair_indice()].get_end();
                std::cout << amplicon_start << ":" << amplicon_end << std::endl;
                tree.insert(Interval(amplicon_start, amplicon_end));

            }
        }
    }
    return tree;
}

int test_itree_overlap(IntervalTree tree, Interval queries[], int num_tests, bool expected[]){
    int result = 0;
    for (int i = 0; i < num_tests; i++)
    {
        result = tree.overlapSearch(queries[i]);
        if (result != expected[i])
        {
            std::cout << "Interval Tree overlap behavior incorrect for interval " << queries[i].low << ":" << queries[i].high 
            << " - " << "Expected: " << expected[i] << " Got: " << result << std::endl;
            return 1;
        }
    }
    return 0;
}


int main(){
    int result = 0;
    std::string pair_indice = "../data/pair_information.tsv";
    std::string bed = "../data/test.bed";
    IntervalTree tree;

    tree = test_populate_amplicons(pair_indice, bed);
    tree.inOrder();
    Interval queries[6] = {Interval(8, 380), Interval(230, 679), Interval(10, 250), Interval(10, 390), Interval(220, 600), Interval(240, 680)};
    bool expected[6] = {true, true, false, false, false, false};
    int num_tests = sizeof(queries) / sizeof(queries[0]);
    result = test_itree_overlap(tree, queries, num_tests, expected);
    return result;
}