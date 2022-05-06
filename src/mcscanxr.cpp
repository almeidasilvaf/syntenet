/*
 * Original Author: Haibao Tang <bao@uga.edu>, May 10, 2007
 *
 * Modified by Author: Yupeng Wang <wyp1125@uga.edu>, Mar 31, 2011
 *
 * Modified by Author: Kristian Ullrich <ullrich@evolbio.mpg.de>, May 5, 2022
 * 
 * Original files can be found here: https://github.com/wyp1125/MCScanX
 * 
 * Combined multiple cpp files into one cpp file: mcscanx.cpp
 * dagchainer.cc, mcscan.cc, msa.cc, out_utils.cc, out_homology.cc,
 * permutation.cc, read_data.cc, read_homology.cc, struct.cc
 * 
 * Original publication:
 * Wang, Yupeng, Haibao Tang, Jeremy D. DeBarry, Xu Tan, Jingping Li,
 * Xiyin Wang, Tae-ho Lee et al.
 * "MCScanX: a toolkit for detection and evolutionary analysis of gene synteny
 * and collinearity." Nucleic acids research 40, no. 7 (2012): e49-e49.
*/

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#include "mcscanxr.h"

/***** MAP *****/
map<string, ortho_stat> cmp_sp;
map<string, Gene_feat> gene_map;
map<string, int> mol_pairs;

/***** VECTOR *****/
vector<Blast_record> match_list;
vector<Seg_feat> seg_list;
vector<more_feat> gene_more;
vector<Score_t> score;
static vector<New_endpoint> endpoints;

/***** SET *****/

geneSet allg;

/***** CONSTANTS *****/

int MATCH_SCORE;
int MATCH_SIZE;
int GAP_PENALTY;
int GAP_SIZE;
int OVERLAP_WINDOW;
double E_VALUE;
int MAX_GAPS;
string PIVOT;
int EXTENSION_DIST;
int CUTOFF_SCORE;
int IN_SYNTENY;
int N_PROXIMAL;
int e_mode;
bool IS_PAIRWISE;
enum { DIAG, UP, LEFT, DEL };
int ali_ct = 0;
int Best_g = -1;
int Best_i;
int Best_j;
int Max_Y;
int max_level;

/***** VOID *****/

/* Print progress message */
void progress(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    vfprintf(stdout, format, args);
    fprintf(stdout, "\n");
    va_end(args);
}

/* Print error message but do not exit */
void err(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[Error] ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
}

/* Print error message but do not exit */
void warn(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[Warn] ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
}

void errAbort(const char *format, ...)
/* Print error message and exit */
{
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[Error] ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
    exit(1);
}

long clock1000()
/* A millisecond clock. */
{
    struct timeval tv;
    static long origSec;
    gettimeofday(&tv, NULL);
    if (origSec == 0) origSec = tv.tv_sec;
    return (tv.tv_sec-origSec)*1000 + tv.tv_usec / 1000;
}

void uglyTime(const char *label, ...)
/* Print label and how long it's been since last call.  Call with
 * a NULL label to initialize. */
{
    static long lastTime = 0;
    long time = clock1000();
    va_list args;
    va_start(args, label);
    if (label != NULL)
    {
        vfprintf(stdout, label, args);
        fprintf(stdout, " [%.3f seconds elapsed]\n", (time - lastTime)/1000.);
    }
    lastTime = time;
    va_end(args);
}

FILE *mustOpen(const char *fileName, const char *mode)
/* Open a file or die */
{
    FILE *f;
    /* gcc 4.2 has "deprecated conversion from string constant" problem"*/
    char *modeName = (char *)"";

    if (sameString(fileName, "stdin")) return stdin;
    if (sameString(fileName, "stdout")) return stdout;
    if ((f = fopen(fileName, mode)) == NULL)
    {
        if (mode)
        {
            if (mode[0] == 'r') modeName = (char *)" to read";
            else if (mode[0] == 'w') modeName = (char *)" to write";
            else if (mode[0] == 'a') modeName = (char *)" to append";
        }
        errAbort("Can't open %s%s: %s", fileName, modeName, strerror(errno));
    }
    return f;
}

// read_data.cpp

// incremental sorting y coord
static bool cmp_y (const Score_t& t1, const Score_t& t2){
    return t1.y < t2.y || (t1.y == t2.y && t1.x < t2.x);
}

// incremental sorting e-value
static bool cmp_ev (const Score_t& t1, const Score_t& t2){
    return t1.score < t2.score;
}

// filter the blast -m8 output by the following threshold:
// lexically sorted, gene #1 < gene #2
// non-self blast match
// both be present in the mcl output file and in the same group
//void read_blast(const char *prefix_fn, bool gff_flag=true){
void read_blast(const char *blast_fn, bool gff_flag=true){
    //char fn[LABEL_LEN], g1[LABEL_LEN], g2[LABEL_LEN];
    char fn[LABEL_LEN];
    //sprintf(fn,"%s.blast",prefix_fn);
    sprintf(fn,"%s",blast_fn);
    ifstream in(fn);
    int i;
    int total_num=0;
    string line,word,geneids,gene1,gene2;
    double evalue;
    map<string, double> blast_map;
    map<string, double>::iterator it;
    cout<<"Reading BLAST file and pre-processing"<<endl;
    while (!in.eof()){
        getline(in,line);
        if (line==""){
            break;
        }
        istringstream test(line);
        getline(test,gene1,'\t');
        getline(test,gene2,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        getline(test,word,'\t');
        istringstream double_iss(word);
        double_iss>>evalue;
        i=gene1.compare(gene2);
        if (i==0){
            continue;
        }else if (i<0){
            geneids=gene1+"&"+gene2;
        }else{
            geneids=gene2+"&"+gene1;
        }
        it = blast_map.find(geneids);
        if (it==blast_map.end()){
            blast_map[geneids]=evalue;
        }else{
            if (evalue<it->second){
                it->second=evalue;
            }
        }
        total_num++;
    }
    in.close();
    //double score;
    Blast_record br;
    int pair_id = 0;
    map<string, Gene_feat>::iterator it1, it2;
    Gene_feat *gf1, *gf2;
    cout<<"Generating BLAST list"<<endl;
    for (it=blast_map.begin(); it!=blast_map.end(); it++){
        istringstream test(it->first);
        getline(test,gene1,'&');
        getline(test,gene2,'&');
        it1 = gene_map.find(gene1);
        it2 = gene_map.find(gene2);
        if (it1==gene_map.end() || it2==gene_map.end()){
            continue;
        }
        gf1 = &(it1->second), gf2 = &(it2->second);
        if (gf1->mol.empty() || gf2->mol.empty()){
            continue;
        }
        if (IN_SYNTENY==1 && gf1->mol.substr(0,2) != gf2->mol.substr(0,2)){
            continue;
        }
/////////////bug here/////////////////////////////////////////////////////////////
        i=gf1->mol.compare(gf2->mol);
//////////////////////////////////////////////////////////////////////////////////
        if (i<0){
            br.gene1=gene1;
            br.gene2=gene2;
            br.mol_pair = gf1->mol+"&"+gf2->mol;
        }else if (i==0){
            if (gf1->mid<=gf2->mid){
                br.gene1=gene1;
                br.gene2=gene2;
            }else{
                br.gene1=gene2;
                br.gene2=gene1;
            }
            br.mol_pair = gf1->mol+"&"+gf2->mol;
        }else{
            br.gene1=gene2;
            br.gene2=gene1;
            br.mol_pair = gf2->mol+"&"+gf1->mol;
        }
//////////////////////////////////////////////////////////////////////////////////
        if(IN_SYNTENY!=2 || gf1->mol.substr(0,2)!=gf2->mol.substr(0,2)){
            mol_pairs[br.mol_pair]++;
        }
        br.pair_id = pair_id++;
        br.score = it->second;
        match_list.push_back(br);
    }
    int selected_num = match_list.size();
    progress("%d matches imported (%d discarded)",
             selected_num, total_num - selected_num);
}

//void read_gff(const char *prefix_fn){
void read_gff(const char *gff_fn){
    char fn[LABEL_LEN];
    string mol,gn,line,word;
    Gene_feat gf;
    //sprintf(fn, "%s.gff", gff_fn);
    sprintf(fn,"%s",gff_fn);
    ifstream in(fn);
    while (!in.eof()){
        getline(in,line);
        if (line==""){
            break;
        }
        istringstream test(line);
        getline(test,mol,'\t');
        gf.mol = mol;
        getline(test,gn,'\t');
        gf.name = gn;
        getline(test,word,'\t');
        gf.mid = atoi(word.c_str());
        gene_map[gf.name] = gf;
    }
    in.close();
}

static void filter_matches_x (){
    // match_bin is a list of records that are potentially repetitive
    vector<Score_t> match_bin, score_cpy;
    vector<Score_t>::const_iterator it, prev_rec;
    sort(score.begin(), score.end());
    prev_rec = it = score.begin();
    it++;
    match_bin.push_back(*(prev_rec));
    for (; it != score.end(); it++){
        // scan whether it has a linking window with previous one
        if ((prev_rec->x != it->x) || (it->y - prev_rec->y) > OVERLAP_WINDOW){
            // record last match_bin, take only least e-value
            score_cpy.push_back(*min_element(match_bin.begin(), match_bin.end(), cmp_ev));
            // start a new match_bin
            match_bin.clear();
        }
        match_bin.push_back(*it);
        prev_rec = it;
    }
    // don't forget the last match_bin
    score_cpy.push_back(*min_element(match_bin.begin(), match_bin.end(), cmp_ev));
    match_bin.clear();
    // copy into score
    score.clear();
    score = score_cpy;
    score_cpy.clear();
}

static void filter_matches_y (){
    // match_bin is a list of records that are potentially repetitive
    vector<Score_t> match_bin, score_cpy;
    vector<Score_t>::const_iterator it, prev_rec;
    sort(score.begin(), score.end(), cmp_y);
    prev_rec = it = score.begin();
    it++;
    match_bin.push_back(*(prev_rec));
    for (; it != score.end(); it++){
        // scan whether it has a linking window with previous one
        if ((prev_rec->y != it->y) || (it->x - prev_rec->x) > OVERLAP_WINDOW){
            // record last match_bin, take only least e-value
            score_cpy.push_back(*min_element(match_bin.begin(), match_bin.end(), cmp_ev));
            // start a new match_bin
            match_bin.clear();
        }
        match_bin.push_back(*it);
        prev_rec = it;
    }
    // don't forget the last match_bin
    score_cpy.push_back(*min_element(match_bin.begin(), match_bin.end(), cmp_ev));
    match_bin.clear();
    // copy into score
    score.clear();
    score = score_cpy;
    score_cpy.clear();
}

// feed into dagchainer
void feed_dag(const string &mol_pair){
    // two additional filters will be applied here
    // best hsp (least e-value)
    // non-repetitive in a window of 50kb region
    vector<Blast_record>::const_iterator it;
    Score_t cur_score;
    for (it = match_list.begin(); it < match_list.end(); it++){
        if (it->mol_pair != mol_pair){
            continue;
        }
        cur_score.pairID = it->pair_id;
        cur_score.x = gene_map[it->gene1].gene_id;
        cur_score.y = gene_map[it->gene2].gene_id;
        cur_score.gene1=it->gene1;
        cur_score.gene2=it->gene2;
        cur_score.score = MATCH_SCORE;
        score.push_back(cur_score);
    }
    // sort by both axis and remove redundant matches within
    // a given window length (default 50kb)
    filter_matches_x();
    filter_matches_y();
    dag_main(score, mol_pair);
}

// permutation.cpp

static double fact(int x)
/* returns factorial of x */
{
    double ans = 1;
    while (x > 1) ans *= x--;
    return ans;
}

static double ln_fact(int x)
/* ln(x!) using Stirling's formula, see Knuth I: 111 */
{
    double dx = x, invx, invx2, invx3, invx5, invx7, sum;

    if (x < 12) return log(fact(x));
    else
    {
        invx = 1 / dx;
        invx2 = invx * invx;
        invx3 = invx2 * invx;
        invx5 = invx3 * invx2;
        invx7 = invx5 * invx2;

        sum = ((dx + 0.5) * log(dx)) - dx;
        sum += log(2 * M_PI) / 2;
        sum += invx / 12 - invx3 / 360;
        sum += invx5/ 1260 - invx7 / 1680;

        return sum;
    }
}

double ln_perm(int n, int r)
/* natural log of permutation number */
{
    if (r > n || r <= 0) return 0;
    return ln_fact(n) - ln_fact(n-r);
}

double ln_comb(int n, int k)
/* natural log of combination number */
{
    if (k <= 0 || k >= n) return 0;
    return ln_fact(n) - ln_fact(k) - ln_fact(n-k);
}

// dagchainer.cpp

// check whether an alignment overlap (tandem alignment)
static bool check_overlap(vector<int>& xx, vector<int>& yy){
    int xmin = *min_element(xx.begin(), xx.end());
    int xmax = *max_element(xx.begin(), xx.end());
    int ymin = *min_element(yy.begin(), yy.end());
    int ymax = *max_element(yy.begin(), yy.end());
    return xmin <= ymax && ymin <= xmax;
}

static void retrieve_pos(int pid, int *pos1, int *pos2){
    /* returns pos1, pos2 for the blast pair */
    Blast_record *match_rec = &match_list[pid];
    *pos1 = gene_map[match_rec->gene1].mid;
    *pos2 = gene_map[match_rec->gene2].mid;
}

static bool is_significant(Seg_feat *sf, vector<Score_t>& score){
    /* test if a syntenic block is significant, see description in permutation.cc */
    /* see formula in permutation.cc, unknowns are m, N, L1, L2, l_1i l_2i*/
    int s1_a, s1_b, s2_a, s2_b, m, N=0, L1, L2, i;
    double l1, l2, summation=0;
    /* get the start and stop coordinates on each syntenic segment */
    s1_a = sf->s1->mid, s1_b = sf->t1->mid;
    s2_a = sf->s2->mid, s2_b = sf->t2->mid;
    /* calculate m, number of anchor points */
    m = sf->pids.size();
    /* calculate N, number of matches in the defined region*/
    vector<Score_t>::const_iterator it;
    for (it=score.begin(); it!=score.end(); it++){
        if (it->x >=s1_a && it->x <=s1_b && it->y >=s2_a && it->y <=s2_b){
            N++;
        }
    }
    /* calculate l1, l2, distance between successive anchor points */
    int l1_pos1, l1_pos2, l2_pos1, l2_pos2;
    retrieve_pos(sf->pids[0], &l1_pos1, &l2_pos1);
    for (i=1; i<m; i++){
        retrieve_pos(sf->pids[i], &l1_pos2, &l2_pos2);
        l1 = fabs(l1_pos2 - l1_pos1);
        l2 = fabs(l2_pos2 - l2_pos1);
        l1_pos1 = l1_pos2;
        l2_pos1 = l2_pos2;
        summation += log(l1)+log(l2);
    }
    /* calculate L1, L2, respective length of the matching region */
    L1 = s1_b - s1_a, L2 = s2_b - s2_a;
    /* this is the formula */
    sf->e_value = exp(M_LN2 + ln_perm(N, m) + summation - (m-1)*(log(L1)+log(L2)));
    return sf->e_value < E_VALUE;
}

static bool Descending_Score(const Path_t &a, const Path_t &b){
    return (a.score > b.score) || (a.score == b.score && a.rc > b.rc);
}

// whether a mol_pair is self comparison, e.g. "at1&at1"
static bool check_self (const string &s){
    int pos = s.find('&');
    return s.substr(0, pos) == s.substr(pos+1);
}

static void print_chains(vector<Score_t>& score, const string &mol_pair){
    /* Find and output highest scoring chains in score treating it as a DAG*/
    vector<float> path_score;
    vector<int> from, ans;
    vector<Path_t> high;
    vector<int> xx, yy;
    Path_t  p;
    bool done;
    int i, j, m, n, s, pid, num_gaps;
    int del_x, del_y;
    double x;
    bool is_self = check_self(mol_pair);
    sort(score.begin(), score.end());
    do{
        done = true;
        n = score.size();
        path_score.resize(n);
        from.resize(n);
        for (i=0; i<n; i++){
            path_score[i] = score[i].score;
            from[i] = -1;
        }
        for (j=1; j<n; j++){
            for (i=j-1; i>=0; i--){
                del_x = score[j].x - score[i].x - 1;
                del_y = score[j].y - score[i].y - 1;
                if  (del_x >= 0 && del_y >= 0){
                    // if (del_x > EXTENSION_DIST && del_y > EXTENSION_DIST)
                    //     break;
                    // if (del_x > EXTENSION_DIST || del_y > EXTENSION_DIST)
                    //     continue;
                    if (del_x > MAX_GAPS){
                        break;
                    }
                    if (del_y > MAX_GAPS){
                        continue;
                    }
                    //num_gaps = MAX(del_x, del_y)/UNIT_DIST;
                    num_gaps = MAX(del_x, del_y);
                    x = path_score[i] + score[j].score;
                    /* gap penalty */
                    if (num_gaps > 0){
                        x +=num_gaps*GAP_PENALTY;
                    }
                    if  (x > path_score[j]){
                        path_score[j] = x;
                        from[j] = i;
                    }
                }
            }
        }
        high.clear();
        for (i=0; i<n; i++){
            if (path_score[i] >= CUTOFF_SCORE){
                p.score = path_score[i];
                p.sub = i;
                p.rc = score[i].x + score[i].y;
                high.push_back(p);
            }
        }
        sort (high.begin(), high.end(), Descending_Score);
        m = high.size();
        for  (i=0; i<m; i++){
            if  (from[high[i].sub] != -2){
                ans.clear();
                for  (j=high[i].sub; from[j]>=0; j=from[j]){
                    ans.push_back(j);
                }
                ans.push_back(j);
                if (from[j] == -2){
                    done = false;
                    break;
                }else{
                    reverse(ans.begin(), ans.end());
                    s = ans.size();
                    if (is_self){
                        for  (j=0; j<s; j++){
                            from[ans[j]] = -2;
                            xx.push_back(score[ans[j]].x);
                            yy.push_back(score[ans[j]].y);
                        }
                    }
                    Seg_feat sf;
                    Blast_record *br;
                    if (!(is_self && check_overlap(xx, yy))){
                        sf.score = path_score[high[i].sub];
                        for (j=0; j<s; j++){
                            from[ans[j]] = -2;
                            pid = score[ans[j]].pairID;
                            br = &match_list[pid];
                            sf.pids.push_back(pid);
                        }
                        /* start and stop positions for two sub-segments */
                        br = &match_list[sf.pids.front()];
                        sf.s1 = &gene_map[br->gene1];
                        sf.s2 = &gene_map[br->gene2];
                        br = &match_list[sf.pids.back()];
                        sf.t1 = &gene_map[br->gene1];
                        sf.t2 = &gene_map[br->gene2];
                        /* determine the orientation of the alignment */
                        sf.sameStrand = *(sf.s2) < *(sf.t2);
                        if (!sf.sameStrand){
                            swap(sf.s2, sf.t2);
                        }
                        sf.mol_pair = mol_pair;
                        /* significance testing */
                        if (is_significant(&sf, score)){
                            seg_list.push_back(sf);
                        }
                    }
                    xx.clear(), yy.clear();
                }
            }
        }
        if (!done){
            for (i=j=0; i<n; i++){
                if (from[i] != -2){
                    if (i!=j){
                        score[j] = score[i];
                    }
                    j++;
                }
            }
            score.resize(j);
        }
    }
    while (!done);
}

void dag_main(vector<Score_t> &score, const string &mol_pair){
    int i, n=score.size();
    // should be sorted by y incremental
    Max_Y = score[n-1].y;
    // forward direction
    print_chains(score, mol_pair);
    // reverse complement the second coordinate set.
    for (i=0; i<n; i++){
        score[i].y = Max_Y - score[i].y + 1;
    }
    // reverse direction
    print_chains(score, mol_pair);
    score.clear();
}

// out_utils.cpp
void print_align(FILE* fw)
/* print pair-wise alignment */
{
    int i, j, pid;
    int nseg = seg_list.size(), nanchor;
    Seg_feat *s;

    print_params(fw);
//////////////////////////////////////////////////
    set<string> colgenes;
    for (i=0; i<nseg; i++)
    {
        s = &seg_list[i];
        nanchor = s->pids.size();
        for (j=0; j<nanchor; j++)
        {
            pid = s->pids[j];
            colgenes.insert(match_list[pid].gene1);
            colgenes.insert(match_list[pid].gene2);
        }
    }
    fprintf( fw, "############### Statistics ###############\n");
    double temp=100*(double)colgenes.size()/(double)gene_map.size();
    fprintf(fw,"# Number of collinear genes: %d, Percentage: %.2f\n",(int)colgenes.size(),temp);
    fprintf(fw,"# Number of all genes: %d\n", (int)gene_map.size());
    fprintf( fw, "##########################################\n");
//////////////////////////////////////////////////
    for (i=0; i<nseg; i++)
    {
        s = &seg_list[i];
        nanchor = s->pids.size();
        fprintf(fw, "## Alignment %d: score=%.1f e_value=%.2g N=%d %s %s\n",
                i, s->score, s->e_value, nanchor, s->mol_pair.c_str(),
                s->sameStrand?"plus":"minus");
        for (j=0; j<nanchor; j++)
        {
            pid = s->pids[j];
            fprintf(fw, "%3d-%3d:\t%s\t%s\t%7.1g\n",
                    i, j, match_list[pid].gene1.c_str(),
                    match_list[pid].gene2.c_str(), match_list[pid].score);
        }
    }
}

void print_params(FILE *fw)
/* print parameters */
{
    fprintf( fw, "############### Parameters ###############\n");
    fprintf( fw, "# MATCH_SCORE: %d\n", MATCH_SCORE );
    fprintf( fw, "# GAP_PENALTY: %d\n", GAP_PENALTY );
    fprintf( fw, "# MATCH_SIZE: %d\n", MATCH_SIZE );
    fprintf( fw, "# OVERLAP_WINDOW: %d\n", OVERLAP_WINDOW );
    fprintf( fw, "# E_VALUE: %lg\n", E_VALUE );
    fprintf( fw, "# MAX GAPS: %d\n", MAX_GAPS );
    fprintf( fw, "# IS_PAIRWISE: %d\n", IS_PAIRWISE );
    fprintf( fw, "# IN_SYNTENY: %d\n", IN_SYNTENY );
}

// out_homology.cpp

void print_align_homology(FILE* fw)
/* print pair-wise alignment */
{
    int i, j, pid;
    int nseg = seg_list.size(), nanchor;
    Seg_feat *s;
    string sp1,sp2,spc;

    print_params(fw);
    //////////////////////////////////////////////////
    set<string> colgenes;
    for (i=0; i<nseg; i++)
    {
        s = &seg_list[i];
        nanchor = s->pids.size();
        for (j=0; j<nanchor; j++)
        {
            pid = s->pids[j];
            colgenes.insert(match_list[pid].gene1);
            colgenes.insert(match_list[pid].gene2);
        }
    }
    fprintf( fw, "############### Statistics ###############\n");
    double temp=100*(double)colgenes.size()/(double)gene_map.size();
    fprintf(fw,"# Number of collinear genes: %d, Percentage: %.2f\n",(int)colgenes.size(),temp);
    fprintf(fw,"# Number of all genes: %d\n", (int)gene_map.size());
    fprintf( fw, "##########################################\n");
//////////////////////////////////////////////////


    for (i=0; i<nseg; i++)
    {
        s = &seg_list[i];
        nanchor = s->pids.size();
        fprintf(fw, "## Alignment %d: score=%.1f e_value=%.2g N=%d %s %s\n",
                i, s->score, s->e_value, nanchor, s->mol_pair.c_str(),
                s->sameStrand?"plus":"minus");
        sp1=(s->s1)->mol.substr(0,2);
        sp2=(s->s2)->mol.substr(0,2);
        spc=sp1+"&"+sp2; 
        cmp_sp[spc].syn_num+=nanchor;
        
        for (j=0; j<nanchor; j++)
        {
            pid = s->pids[j];
            fprintf(fw, "%3d-%3d:\t%s\t%s\t%7.1g\n",
                    i, j, match_list[pid].gene1.c_str(),
                    match_list[pid].gene2.c_str(), match_list[pid].score);
        }
    }
}

void print_params_homology(FILE *fw)
/* print parameters */
{
    fprintf( fw, "############### Parameters ###############\n");
    fprintf( fw, "# MATCH_SCORE: %d\n", MATCH_SCORE );
    fprintf( fw, "# GAP_PENALTY: %d\n", GAP_PENALTY );
    fprintf( fw, "# MATCH_SIZE: %d\n", MATCH_SIZE );
    fprintf( fw, "# OVERLAP_WINDOW: %d\n", OVERLAP_WINDOW );
    fprintf( fw, "# E_VALUE: %lg\n", E_VALUE );
    fprintf( fw, "# MAX GAPS: %d\n", MAX_GAPS );
    fprintf( fw, "# IS_PAIRWISE: %d\n", IS_PAIRWISE );
    fprintf( fw, "# IN_SYNTENY: %d\n", IN_SYNTENY );
}

// msa.cpp

void get_endpoints()
{
    int n=seg_list.size();
    int i;
    //int j;
    Seg_feat *s;
    for (i=0;i<n;i++)
    {
        s = &seg_list[i];
        s->index=i;
        New_endpoint ep;

        ep.n=s->s1;
        ep.index=2*i;
        ep.start=true;
        ep.e=s->t1;
        endpoints.push_back(ep);
        ep.n=s->t1;
        ep.index=2*i;
        ep.start=false;
        ep.e=s->s1;
        endpoints.push_back(ep);

        ep.n=s->s2;
        ep.index=2*i+1;
        ep.start=true;
        ep.e=s->t2;
        endpoints.push_back(ep);
        ep.n=s->t2;
        ep.index=2*i+1;
        ep.start=false;
        ep.e=s->s2;
        endpoints.push_back(ep);
    }
    sort(endpoints.begin(), endpoints.end());
}

void add_block(Gene_feat* s, Gene_feat* t, int level)
{
    int i,j;
    i=0;
    geneSet::iterator it,it1,it11;
    it=allg.find(s);
    it1=allg.find(t);
    it11=allg.end();
    it11--;
    for (;(*it)->mid<=(*it1)->mid&&(*it)->mol==(*it1)->mol;it++)
    {
        i=(*it)->cursor.size();
        if (i<level)
        {
            for (j=i+1;j<level;j++)
            {
                (*it)->cursor.push_back(0);
            }
            (*it)->cursor.push_back(1);
        }
        else
        {
            (*it)->cursor[level-1]=1;
        }
        if (it==it11)
        {
            break;
        }
    }
}

void add_matchpoints(int seg_index,int level)
{
    Seg_feat *s;
    int i;
    int j;
    //int k;
    map<string, Gene_feat>::iterator it5;
    i=(int)(seg_index/2);
    s=&seg_list[i];
    if (seg_index%2==0)
    {
        for (j=0;j<s->pids.size();j++)
        {
            it5=gene_map.find(match_list[s->pids[j]].gene1);
            if (it5->second.cursor.size()>=level)
                it5->second.cursor[level-1]=s->pids[j]+2;
        }
    }
    else
    {
        for (j=0;j<s->pids.size();j++)
        {
            it5=gene_map.find(match_list[s->pids[j]].gene2);
            if (it5->second.cursor.size()>=level)
                it5->second.cursor[level-1]=-(s->pids[j]+2);
        }
    }
}

void traverse()
{
    int i,j,k,lev;
    Gene_feat gf4;
    //char* temp;
    map<string, Gene_feat>::iterator it4;
    for (i=0;i<endpoints.size();i++)
    {
        if (endpoints[i].start==1)
        {
            it4=gene_map.find(endpoints[i].n->name);
            gf4=it4->second;
            k=gf4.cursor.size();
            if (k==0)
            {
                add_block(endpoints[i].n,endpoints[i].e,1);
                add_matchpoints(endpoints[i].index,1);
            }

            else
            {
                for (j=0;j<k;j++)
                {
                    if (gf4.cursor[j]==0)
                    {
                        lev=j+1;
                        break;
                    }
                }
                if (j==k)
                    lev=j+1;
                add_block(endpoints[i].n,endpoints[i].e,lev);
                add_matchpoints(endpoints[i].index,lev);
                if (lev>max_level)
                    max_level=lev;
            }
        }
        else
        {
            ;
        }
    }
}

void mark_tandem(const char *prefix_fn)
{
    int i,j;
    i=0;
    geneSet::const_iterator it7=allg.begin();
    more_feat mf;
    for (; it7!=allg.end(); it7++)
    {
        //(*it7)->gene_id=i;
        mf.depth=0;
        mf.tandem=0;
        for (j=0;j<(*it7)->cursor.size();j++)
        {
            if ((*it7)->cursor[j]!=0)
            {
                mf.depth++;
            }
        }
        gene_more.push_back(mf);
        i++;
    }
    vector<string>tpair1;
    vector<string>tpair2;
    map<string, Gene_feat>::iterator it8,it9;
    for (i=0;i<match_list.size();i++)
    {
        it8=gene_map.find(match_list[i].gene1);
        it9=gene_map.find(match_list[i].gene2);
        if (fabs(it8->second.gene_id-it9->second.gene_id)==1&&it8->second.mol==it9->second.mol)
        {
            gene_more[it8->second.gene_id].tandem=1;
            gene_more[it9->second.gene_id].tandem=1;
            tpair1.push_back(it8->second.name);
            tpair2.push_back(it9->second.name);
        }
    }
    if(tpair1.size()>0)
    {
    ofstream result;
    char fn[LABEL_LEN];
    sprintf(fn,"%s.tandem",prefix_fn);
    cout<<"Tandem pairs written to "<<fn<<endl;
    result.open(fn,ios::out);   
    for(i=0;i<tpair1.size();i++)
    {
    //cout<<tpair1[i]<<","<<tpair2[i]<<endl;
    result<<tpair1[i]<<","<<tpair2[i]<<endl;
    }
    result.close();
    }
}

void print_html()
{
    int i;
    int j;
    //int k;
    string color;
    ofstream result;
    string prev_mol="";
    Gene_feat *n;
    char result_dir[200];
    geneSet::iterator it6;
    i=0;
    for (it6=allg.begin();it6!=allg.end();it6++)
    {
        n=(*it6);
        if (n->mol!=prev_mol)
        {
            if (i>0)
            {
                result<<"</table></html>";
                result.close();
            }
//sprintf(result_dir,"%s.html/%s.html",prefix_fn,n->mol.c_str());
            sprintf(result_dir,"%s.html",n->mol.c_str());
            cout<<result_dir<<endl;
            result.open(result_dir,ios::out);
            result<<"<html><table cellspacing='0' cellpadding='0' align='left'>";
            result<<"<tr align='center'><td>Duplication depth</td><td>&nbsp;&nbsp;Reference chromosome</td><td align='left' colspan='"<<2*max_level<<"'>&nbsp;&nbsp;Collinear blocks</td></tr>"<<endl;
            prev_mol=n->mol;
//cout<<prev_mol<<endl;
            i++;
        }
        color="'#dddddd'";
        if (gene_more[n->gene_id].tandem)
            color="'#ee0000'";
//if(gene_more[n->gene_id].break_point)
//color="'#0000ee'";
        result<<"<tr align='center'><td>"<<gene_more[n->gene_id].depth<<"</td><td bgcolor="<<color<<">"<<n->name<<"</td>";
        for (j=0;j<n->cursor.size();j++)
        {
            result<<"<td>&nbsp;&nbsp;</td>";
            if (n->cursor[j]==0)
            {
                result<<"<td>&nbsp;</td>";
            }
            else if (n->cursor[j]==1)
            {
                result<<"<td>|&nbsp;|</td>";
            }
            else if (n->cursor[j]>1)
            {
                result<<"<td bgcolor='#ffff99'>"<<match_list[n->cursor[j]-2].gene2<<"</td>";
            }
            else
            {
                result<<"<td bgcolor='#ffff99'>"<<match_list[-n->cursor[j]-2].gene1<<"</td>";
            }
        }
        for (j=n->cursor.size();j<max_level;j++)
        {
            result<<"<td>&nbsp;</td>";
        }
        result<<"</tr>"<<endl;
    }
    result<<"<html><table>"<<endl;
    result.close();
}

void print_test()
{
    int j=0;
    geneSet::const_iterator i=allg.begin();
    for (; i!=allg.end(); i++)
    {
        cout<<(*i)->name.c_str()<<" "<<(*i)->cursor.size()<<endl;
        j++;
        if (j>1000)
            break;
    }

}

void msa_main(const char *prefix_fn)
{
    max_level=1;
//    fill_allg();
    get_endpoints();
    traverse();
    mark_tandem(prefix_fn);
    char html_fn[LABEL_LEN];
    printf("Writing multiple syntenic blocks to HTML files\n");
    sprintf(html_fn,"%s.html",prefix_fn);
    if (chdir(html_fn)<0){
        mkdir(html_fn,S_IRWXU|S_IRGRP|S_IXGRP);
        chdir(html_fn);
    }
    print_html();
}

void fill_allg()
{
    Gene_feat *gf1;
    map<string, Gene_feat>::iterator it;
    for (it=gene_map.begin(); it!=gene_map.end(); it++)
    {
        gf1 = &(it->second);
        allg.insert(gf1);
    }
    int i=0;
    geneSet::const_iterator it77=allg.begin();   
    for (; it77!=allg.end(); it77++)
    {
    (*it77)->gene_id=i;
    i++;
    }
}

//' @useDynLib syntenet, .registration = TRUE
//' @import Rcpp
//' @title rcpp_mcscanx_file
//' @name rcpp_mcscanx_file
//' @description MCSanX provides a clustering module for viewing the
//' relationship of colinear segments in multiple genomes (or heavily redundant
//' genomes). It takes the predicted pairwise segments from dynamic programming
//' (DAGchainer in particular) and then try to build consensus segments from a
//' set of related, overlapping segments.
//' @return list
//' @param blast_file blast input
//' @param gff_file gff input
//' @param prefix output prefix (default: out)
//' @param outdir output directory (default: "")
//' @param match_score match score (default: 50)
//' @param gap_penalty gap penalty (default: -1)
//' @param match_size match_size (default: 5)
//' @param e_value e_value (default: 1e-5)
//' @param max_gaps max gaps (default: 25)
//' @param overlap_window overlap window (default: 5)
//' @param is_pairwise specify if only pairwise blocks should be reported
//' (default: FALSE)
//' @param is_synteny specify patterns of collinear blocks.
//' 0: intra- and inter-species (default); 1: intra-species; 2: inter-species
//' @references Wang et al. (2012) MCScanX: a toolkit for detection and
//' evolutionary analysis of gene synteny and collinearity.
//' \emph{Nucleic acids research}. \bold{40.7}, e49-e49.
//' @references Haas et al. (2004) DAGchainer: a tool for mining segmental
//' genome duplications and synteny. \emph{Bioinformatics}. \bold{20.18}
//' 3643-3646.
//' @export rcpp_mcscanx_file
//' @author Kristian K Ullrich
// [[Rcpp::export]]
int rcpp_mcscanx_file(
    std::string blast_file,
    std::string gff_file,
    std::string prefix="out",
    std::string outdir="",
    int match_score=50,
    int gap_penalty=-1,
    int match_size=5,
    double e_value=1e-5,
    int max_gaps=25,
    int overlap_window=5,
    bool is_pairwise=false,
    int in_synteny=0){
    char curwd[256];
    getcwd(curwd, 256);
    MATCH_SCORE = match_score;
    GAP_PENALTY = gap_penalty;
    MATCH_SIZE = match_size;
    E_VALUE = e_value;
    MAX_GAPS = max_gaps;
    OVERLAP_WINDOW = overlap_window;
    IS_PAIRWISE = is_pairwise;
    IN_SYNTENY = in_synteny;
    CUTOFF_SCORE = MATCH_SCORE*MATCH_SIZE;
    /* Start the timer */
    uglyTime(NULL);
    char align_fn[LABEL_LEN];
    FILE *fw;
    const char *prefix_fn = prefix.c_str();
    const char *blast_fn = blast_file.c_str();
    const char *gff_fn = gff_file.c_str();
    read_gff(gff_fn);
    read_blast(blast_fn);
    if (outdir!=""){
        const char *outdir_fn = outdir.c_str();
        if (chdir(outdir_fn)<0){
            mkdir(outdir_fn,S_IRWXU|S_IRGRP|S_IXGRP);
            chdir(outdir_fn);
        }
    }
    sprintf(align_fn, "%s.collinearity", prefix_fn);
    fw = mustOpen(align_fn, "w");

    progress("%d pairwise comparisons", (int) mol_pairs.size());
    fill_allg();
    map<string, int>::const_iterator ip;
    for (ip=mol_pairs.begin(); ip!=mol_pairs.end(); ip++)
    {
        if (ip->second >= MATCH_SIZE) feed_dag(string(ip->first));
    }

    progress("%d alignments generated", (int) seg_list.size());
    print_align(fw);
    fclose(fw);
    uglyTime("Pairwise collinear blocks written to %s", align_fn);

    if (IS_PAIRWISE){
        seg_list.clear();
        match_list.clear();
        gene_more.clear();
        score.clear();
        endpoints.clear();
        chdir(curwd);
        /* End the timer */
        uglyTime("Done!");
        return 0;
    }
    msa_main(prefix_fn);
    seg_list.clear();
    match_list.clear();
    gene_more.clear();
    score.clear();
    endpoints.clear();
    chdir(curwd);
    /* End the timer */
    uglyTime("Done!");
    return 0;
}

