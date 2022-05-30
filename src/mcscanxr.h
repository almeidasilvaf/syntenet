/*
 * Original Author: Haibao Tang <bao@uga.edu>, May 10, 2007
 *
 * Modified by Author: Yupeng Wang <wyp1125@uga.edu>, Mar 31, 2011
 *
 * Modified by Author: Kristian Ullrich <ullrich@evolbio.mpg.de>, May 29, 2022
 * 
 * Original files can be found here: https://github.com/wyp1125/MCScanX
 * 
 * Combined multiple header files into one header file: mcscanx.h
 * dagchainer.h, mcscan.h, msa.h, out_utils.h, permutation.h, read_data.h,
 * read_homology.h, struct.h
 * 
 * Original publication:
 * Wang, Yupeng, Haibao Tang, Jeremy D. DeBarry, Xu Tan, Jingping Li,
 * Xiyin Wang, Tae-ho Lee et al.
 * "MCScanX: a toolkit for detection and evolutionary analysis of gene synteny
 * and collinearity." Nucleic acids research 40, no. 7 (2012): e49-e49.
*/

#ifndef __MCSCANX_H
#define __MCSCANX_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <cstdio>
#include <cstring>
#include <cstdarg>
#include <cerrno>
#include <climits>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <getopt.h>
#include <string>
#include <vector>
#include <list>
#include <queue>
#include <set>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define sameString(a, b) (strcmp((a), (b))==0)
#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))
#define LABEL_LEN 200

#ifndef __GNUC__
#define __attribute__(x)
#endif

/***** STRUCT *****/

struct Blast_record{
    std::string gene1, gene2;
    std::string mol_pair;
    int pair_id;
    int node;
    double score;
};

struct more_feat{
    int tandem;
    int depth;
};

struct Gene_feat{
    std::vector<int> cursor;
    std::string name;
    std::string mol;
    int mid;
    int gene_id;
    bool operator < (const Gene_feat &g) const{
        return (mol == g.mol && mid < g.mid) || mol < g.mol || (mol == g.mol && mid == g.mid && name.compare(g.name)<0);
    }
};

struct geneCmp {bool operator() (const Gene_feat *a, const Gene_feat *b) const{
        return (a->mol == b->mol && a->mid < b->mid) || a->mol < b->mol || (a->mol == b->mol && a->mid == b->mid && a->name.compare(b->name)<0);
    }
};

struct Seg_feat{
    std::vector<int> pids;
    int index;
    Gene_feat *s1, *t1, *s2, *t2;
    double score, e_value;
    std::string mol_pair;
    bool sameStrand;
};

struct Cell_t{
    float raw;
    int score : 30;
    unsigned from : 2;
};

struct Score_t{
    int pairID;
    int x, y;
    float score;
    std::string gene1;
    std::string gene2;
    bool operator < (const Score_t & node) const{
        return  (x < node.x) || (x == node.x && y < node.y);
    }
};

struct Path_t{
    float score;
    int rc;
    int sub;
};

struct ortho_stat{
    int all_num;
    int syn_num;
};

struct New_endpoint{
    Gene_feat *n;
    int index;
    bool start;
    Gene_feat *e;
    bool operator <  (const New_endpoint &g) const
    {
        return (n->mol == g.n->mol && n->mid < g.n->mid) || n->mol < g.n->mol;
    }
};

typedef std::set<Gene_feat *, geneCmp> geneSet;

/***** Helper functions (Some from James Kent library) *****/

FILE *mustOpen(const char *fileName, const char *mode);

void read_blast(const std::string blast_infile);
void read_gff(const std::string gff_infile);
void feed_dag(const std::string &mol_pair);
void dag_main(std::vector<Score_t> &score, const std::string &mol_pair);

double ln_perm(int n, int r);
double ln_comb(int n, int k);

void print_params(FILE *fw);
void print_align(FILE *fw);

#endif
