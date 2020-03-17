/****************************************************************************
 *
 *     LongAGE -- Ultra-long Alignment with Gap Excision
 *     Copyright (C)  Alexej Abyzov
 *                
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Creative Commons license
 * (Attribution-NonCommerical).
 * See terms at http://creativecommons.org/licenses/by-nc/2.5/legalcode
 *
 * Author: Alexej Abyzov
 */

#include <stdint.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <limits>

#include "LongAGE.h"
using namespace std;

LongAGE::LongAGE(LongAGEParameters paras) {
    scr = paras.scr;
    swap = false;

    if (paras.seq1.length() >= paras.seq2.length()) {
        seq1 = paras.seq1;
        seq2 = paras.seq2;
    } else {
        swap = true;
        seq1 = paras.seq2;
        seq2 = paras.seq1;
    }
}

pair<pair<long, long>,pair<long,long>> LongAGE::find_BPs(Scorer* scr, bool tdup_inv) {
    if (tdup_inv) {
        return find_inv_tdup_BPs(scr);
    } else {
        return find_indel_BPs();
    }    
}

pair<pair<long, long>,pair<long,long>> LongAGE::find_inv_tdup_BPs(Scorer* scr) {
    pair<pair<ScorePosition, ScorePosition>,pair<ScorePosition,ScorePosition> > c1, c2;
    pair<ScorePosition, ScorePosition> st_pair, ed_pair;
    string tmp = seq1;
    seq1 = seq2;
    seq2 = tmp;
    len1 = seq1.length();
    len2 = seq2.length();
    score_n2 = len2 + 2;
    int *M1;
    M1 = new int[2*score_n2];
    local_run(len1, 1, len2, M1);
    c1 = get_candidatate_BPs(scr, M1);
    long local_score_first = c1.first.second.score + c1.second.second.score;

    tmp = seq1;
    seq1 = seq2;
    seq2 = tmp;
    len1 = seq1.length();
    len2 = seq2.length();
    score_n2 = len2 + 2;
    int *M2;
    M2 = new int[2*score_n2];
    local_run(len1, 1, len2, M2);
    c2 = get_candidatate_BPs(scr, M2);
    long local_score_second = c2.first.second.score + c2.second.second.score;

    if (local_score_first > local_score_second) {
        _score = local_score_first;
        tmp = seq1;
        seq1 = seq2;
        seq2 = tmp;
        len1 = seq1.length();
        len2 = seq2.length();
        score_n2 = len2 + 2;
        st_pair = c1.first;
        ed_pair = c1.second;
    } else {
        _score = local_score_second;
        st_pair = c2.first;
        ed_pair = c2.second;
    }
    return make_pair(make_pair(st_pair.first.position, st_pair.second.position), 
                     make_pair(ed_pair.first.position, ed_pair.second.position));
}

pair<pair<long, long>,pair<long,long>> LongAGE::find_indel_BPs() {
    len1 = seq1.length();
    len2 = seq2.length();
    score_n2 = len2 + 2;
    int *M;
    M = new int[2*score_n2];

    pair<ScorePosition, pair<ScorePosition, long> > first_results = local_run(len1, 1, len2, M);
    ScorePosition left_end = first_results.first;
    ScorePosition right_end = first_results.second.first;

    long right_start = first_results.second.second + score_n2 + 1;
    _score = right_end.score;

    long right_start_n = split(right_start).n;
    long right_start_m = split(right_start).m;
    long left_end_n = split(left_end.position).n;
    long left_end_m = split(left_end.position).m;

    if (left_end_n >= right_start_n || left_end_m >= right_start_m) {
        pair<ScorePosition, pair<ScorePosition, long> > contained_results = local_run(right_start_n - 1, 1, right_start_m - 1, M);
        left_end = contained_results.first;
    }

    ScorePosition left_start = get_fisrt_BP(left_end);

    pair<long, long> left = make_pair(left_start.position, left_end.position);
    pair<long, long> right = make_pair(right_start, right_end.position);
    return make_pair(left, right);
}

void LongAGE::findAlignment(LongAGEParameters paras) {

    pair<pair<long, long>, pair<long, long> > bps = find_BPs(scr, paras.tdup_inv);

    long l_st = bps.first.first;
    long l_ed = bps.first.second;
    long r_st = bps.second.first;
    long r_ed = bps.second.second;
    Location l_st_s = split(l_st);
    Location l_ed_s = split(l_ed);
    Location r_st_s = split(r_st);
    Location r_ed_s = split(r_ed);

    long left_s1_st = l_st_s.n - 1;
    long left_s1_ed = l_ed_s.n - 1;

    long left_s2_st = l_st_s.m - 1;
    long left_s2_ed = l_ed_s.m - 1;

    long right_s1_st = r_st_s.n - 1;
    long right_s1_ed = r_ed_s.n - 1;

    long right_s2_st = r_st_s.m - 1;
    long right_s2_ed = r_ed_s.m - 1;

    string left_s1 = seq1.substr(left_s1_st, left_s1_ed-left_s1_st+1);
    string left_s2 = seq2.substr(left_s2_st, left_s2_ed-left_s2_st+1);
    string right_s1 = seq1.substr(right_s1_st, right_s1_ed-right_s1_st+1);
    string right_s2 = seq2.substr(right_s2_st, right_s2_ed-right_s2_st+1);

    pair<string, string> left_result = get_aligned_seqs(left_s1, left_s2);
    pair<string, string> right_result = get_aligned_seqs(right_s1, right_s2);

    AliFragment *f1, *f2;
    if (swap) {
        f1 = new AliFragment(left_result.second, left_result.first,
                paras.start2 + paras.inc2*(left_s2_st), paras.start1 + paras.inc1*(left_s1_ed),
                paras.start2 + paras.inc2*(left_s2_ed), paras.start1 + paras.inc1*(left_s1_ed));
        f2 = new AliFragment(right_result.second, right_result.first,
                paras.start2 + paras.inc2*(right_s2_st),  paras.start1 + paras.inc1*(right_s1_st),
                paras.start2 + paras.inc2*(right_s2_ed), paras.start1 + paras.inc1*(right_s1_ed));
    } else {
        f1 = new AliFragment(left_result.first, left_result.second,
                paras.start1 + paras.inc1*(left_s1_st), paras.start2 + paras.inc2*(left_s2_st),
                paras.start1 + paras.inc1*(left_s1_ed), paras.start2 + paras.inc2*(left_s2_ed));
        f2 = new AliFragment(right_result.first, right_result.second,
                paras.start1 + paras.inc1*(right_s1_st), paras.start2 + paras.inc2*(right_s2_st),
                paras.start1 + paras.inc1*(right_s1_ed), paras.start2 + paras.inc2*(right_s2_ed));
    }

    f1->next(f2);
    _frags = f1;
}


const int MIN =  5000 - std::numeric_limits<std::int32_t>::max();
const int MIN_EPS = 10000 - std::numeric_limits<std::int32_t>::max();;

inline int Max(int a, int b) {
    int local_max = a > b ? a : b;
    return local_max < MIN_EPS ? MIN : local_max;
}

inline int Max(int a, int b, int c) {
    return Max(Max(a, b), c);
}

inline int Max(int a, int b, int c, int d) {
    return Max(Max(a, b, c), d);
}

inline int Max(int a, int b, int c, int d, int e) {
    return Max(Max(a, b, c, d), e);
}

int *res;
int nres;
char *seqc0, *seqc1;
int n0, n1, nd;

static int nmax=0;
static int g, hh, m;

#define gap(k)  ((k) <= 0 ? 0 : g+hh*(k))

static int *sapp;
static int  last;

inline void DEL(int k) { (last < 0) ? last = sapp[-1] -= (k) : last = *sapp++ = -(k); }

inline void INS(int k) { (last < 0) ? ({ sapp[-1] = (k); *sapp++ = last; }) : last = *sapp++ = (k); }

inline void REP() { last = *sapp++ = 0; }

static int *CC, *DD;
static int *RR, *SS;

int compute_consensus(char *aa0, int n0, char *aa1, int n1, int *res)
{
    int i0, i1;
    int op, nc;
    char *sp0, *sp1;
    int *rp;

    sp0 = seqc0;
    sp1 = seqc1;
    rp = res;
    nc = nd = i0 = i1 = op = 0;

    while (i0 < n0 || i1 < n1) {
        if (op == 0 && *rp == 0) {
            op = *rp++;
            *sp0 = aa0[i0++];
            *sp1 = aa1[i1++];
            nc++;
            if (*sp0++ == *sp1++) nd++;
        }
        else {
            if (op==0)
                op = *rp++;
            if (op>0) {
                *sp0++ = '-';
                *sp1++ = aa1[i1++];
                op--;
                nc++;
            }
            else {
                *sp0++ = aa0[i0++];
                *sp1++ = '-';
                op++;
                nc++;
            }
        }
    }
    return nc;
}

void create_aligned_seq(int seqsiz)
{
    res = (int *) calloc((size_t) seqsiz, sizeof(int));
    seqc0 = (char*) calloc((size_t) seqsiz, sizeof(char));
    seqc1 = (char*)calloc((size_t) seqsiz, sizeof(char));
    if (res == NULL || seqc0 == NULL || seqc1 == NULL) {
        fprintf(stderr, "Cannot allocate consensus arrays %d\n", seqsiz);
        exit(1);
    }
}

static int align(char *A, char*B,
                 int M, int N,
                 int tb, int te,
                 int match, int mismatch)
{
    int midi, midj, type;
    int midc;
    int c1, c2;
    {
        int   i, j;
        int c, e, d, s;
        int t;

        if (N <= 0)
        {
            if (M > 0) DEL(M);
            return -gap(M);
        }
        if (M <= 1)
        {
            if (M <= 0)
            {
                INS(N);
                return -gap(N);
            }

            if (tb < te) tb = te;
            midc = (tb-hh) - gap(N);
            midj = 0;

            for (j = 1; j <= N; j++)
            {
                if (B[j] == A[1]) {
                    c = -gap(j-1) + match - gap(N-j);
                } else {
                    c = -gap(j-1) + mismatch - gap(N-j);
                }
                
                if (c > midc)
                {
                    midc = c;
                    midj = j;
                }
            }
            if (midj == 0)
            {
                INS(N);
                DEL(1);
            }
            else
            {
                if (midj > 1) INS(midj-1);
                REP();
                if (midj < N) INS(N-midj);
            }
            return midc;
        }

        midi = M/2;
        CC[0] = 0;
        t = -g;
        for (j = 1; j <= N; j++)
        {
            CC[j] = t = t-hh;
            DD[j] = t-g;
        }
        t = tb;
        for (i = 1; i <= midi; i++)
        {
            s = CC[0];
            CC[0] = c = t = t-hh;
            e = t-g;
            for (j = 1; j <= N; j++)
            {
                if ((c =   c   - m) > (e =   e   - hh)) e = c;
                if ((c = CC[j] - m) > (d = DD[j] - hh)) d = c;
                c = B[j] == A[i] ? s + match : s + mismatch;
                if (e > c) c = e;
                if (d > c) c = d;
                s = CC[j];
                CC[j] = c;
                DD[j] = d;
            }
        }
        DD[0] = CC[0];

        RR[N] = 0;
        t = -g;
        for (j = N-1; j >= 0; j--)
        {
            RR[j] = t = t-hh;
            SS[j] = t-g;
        }
        t = te;
        for (i = M-1; i >= midi; i--)
        {
            s = RR[N];
            RR[N] = c = t = t-hh;
            e = t-g;
            for (j = N-1; j >= 0; j--)
            {
                if ((c =   c   - m) > (e =   e   - hh)) e = c;
                if ((c = RR[j] - m) > (d = SS[j] - hh)) d = c;
                c = B[j+1] == A[i+1] ? s + match : s + mismatch;
                if (e > c) c = e;
                if (d > c) c = d;
                s = RR[j];
                RR[j] = c;
                SS[j] = d;
            }
        }
        SS[N] = RR[N];

        midc = CC[0]+RR[0];
        midj = 0;
        type = 1;
        for (j = 0; j <= N; j++)
            if ((c = CC[j] + RR[j]) >= midc)
            if (c > midc || CC[j] != DD[j] && RR[j] == SS[j])
            {
                midc = c;
                midj = j;
            }
        for (j = N; j >= 0; j--)
            if ((c = DD[j] + SS[j] + g) > midc)
            {
                midc = c;
                midj = j;
                type = 2;
            }
    }

    if (type == 1)
    {
        c1 = align(A, B, midi, midj, tb, -g, match, mismatch);
        c2 = align(A+midi, B+midj, M-midi, N-midj, -g, te, match, mismatch);
    }
    else
    {
        align(A, B, midi-1, midj, tb, 0, match, mismatch);
        DEL(2);
        align(A+midi+1, B+midj, M-midi-1, N-midj, 0, te, match, mismatch);
    }
    return midc;
}

void allocate(int M, int N, int match, int mismatch, int G, int H, int S[], int* NC)
{
    int c, ck;
    g = G;
    hh = H;
    m = g+hh;
    sapp = S;
    last = 0;

    if (CC==NULL) {
        nmax = N;
        CC=(int *)calloc((size_t)(nmax+1),sizeof(int));
        DD=(int *)calloc((size_t)(nmax+1),sizeof(int));
        RR=(int *)calloc((size_t)(nmax+1),sizeof(int));
        SS=(int *)calloc((size_t)(nmax+1),sizeof(int));
    }
    else if (N > nmax ) {
        nmax = N;
        CC=(int *)realloc(CC,(size_t)(nmax+1)*sizeof(int));
        DD=(int *)realloc(DD,(size_t)(nmax+1)*sizeof(int));
        RR=(int *)realloc(RR,(size_t)(nmax+1)*sizeof(int));
        SS=(int *)realloc(SS,(size_t)(nmax+1)*sizeof(int));
    }

    if (CC==NULL || DD==NULL || RR==NULL || SS==NULL) {
        fprintf(stderr," cannot allocate llmax arrays\n");
        exit(1);
    }
}

void run_align(char **rs, int len1, int len2,
               char *p1, char *p2,
               int gap_open, int gap_extend, int match, int mismatch) {
    n0 = len1; n1 = len2;
    create_aligned_seq(n0+n1);
    allocate(n0, n1, match, mismatch, -gap_open, -gap_extend, res, &nres);
    align(&p1[-1], &p2[-1], n0, n1,-g,-g,match,mismatch);
    compute_consensus(p1, n0, p2, n1, res);
    rs[0] = seqc0; rs[1] = seqc1;
}

LongAGEParameters::LongAGEParameters() {}

pair<string, string> LongAGE::get_aligned_seqs(string s1, string s2) {
    char* rs[2];
    char* p1 = new char[s1.length()];
    char* p2 = new char[s2.length()];

    strcpy(p1, s1.c_str());
    strcpy(p2, s2.c_str());
    for (long i=0; i<s1.length(); i++) {
        p1[i] = toupper(p1[i]);
    }
    for (long i=0; i<s2.length(); i++) {
        p2[i] = toupper(p2[i]);
    }
    run_align(rs, s1.length(), s2.length(), p1, p2,
              scr->getGapOpen(), scr->getGapExtend(),
              scr->getMatch(), scr->getMismatch());
    string res1(rs[0]);
    string res2(rs[1]);
    return make_pair(res1, res2);
}

long LongAGE::Score() {
    return _score;
}

AliFragment* LongAGE::Fragments() {
    return _frags;
}

LongAGE::~LongAGE() {}

pair< pair<ScorePosition, ScorePosition>,pair<ScorePosition,ScorePosition>> 
LongAGE::get_candidatate_BPs(Scorer* scr, int* M) {

    string s1_rev = seq1;
    string s2_rev = seq2;

    reverse(s1_rev.begin(), s1_rev.end());
    reverse(s2_rev.begin(), s2_rev.end());

    string hold1 = seq1;
    string hold2 = seq2;
    seq1 = s1_rev;
    seq2 = s2_rev;
    int M_rev[2*score_n2];
    local_run(len1, 1, len2, M_rev);

    seq1 = hold1;
    seq2 = hold2;
    int M_back[2*score_n2];
    for (long i=0, j=score_n2-1; i<score_n2 && j>=0; i++, j--) {
        M_back[i] = M_rev[j];
    }
    long max_m = -1;
    long max_m_score = -1;
    for (long m=1; m<score_n2-2; m++) {
        long check = M[m] + M_back[m+1];
        if (check > max_m_score) {
            max_m_score = check;
            max_m = m;
        }
    }

    pair<ScorePosition, pair<ScorePosition, long> > top = local_run(len1, 1, max_m, M);
    pair<ScorePosition, pair<ScorePosition, long> > bottom = local_run(len1, max_m+1, len2, M);
    ScorePosition top_end = top.first;
    ScorePosition bottom_end = bottom.first;
    ScorePosition top_start = get_fisrt_BP(top_end);
    ScorePosition bottom_start = get_fisrt_BP(bottom_end);

    pair<ScorePosition, ScorePosition> top_pair = make_pair(top_start, top_end);
    pair<ScorePosition, ScorePosition> bottom_pair = make_pair(bottom_start, bottom_end);
    return make_pair(top_pair, bottom_pair);
}


Location LongAGE::split(long len_n_m) {

    long n = len_n_m / score_n2;
    long m = len_n_m - (n*score_n2);
    return Location(n,m);
}

pair<ScorePosition, pair<ScorePosition, long> > LongAGE::local_run(long n_stop, long m_start, long m_stop, int* M) {

    ScorePosition left, right;
    long last_m = -1;

    int *S1, *X1, *Y1;
    S1 = new int[2*score_n2];
    X1 = new int[2*score_n2];
    Y1 = new int[2*score_n2];

    long *R_S, *R_X, *R_Y;
    R_S = new long[2*score_n2];
    R_X = new long[2*score_n2];
    R_Y = new long[2*score_n2];

    int *S2, *X2, *Y2;
    S2 = new int[2*score_n2];
    X2 = new int[2*score_n2];
    Y2 = new int[2*score_n2];

    for (long i=0; i<2*score_n2; i++) {
        R_S[i] = -1;
        R_X[i] = -1;
        R_Y[i] = -1;
        M[i] = 0;
        S1[i] = 0;
        S2[i] = 0;
        X1[i] = 0;
        Y1[i] = 0;
        X2[i] = 0;
        Y2[i] = 0;
    }

    for (long i=0; i<score_n2; i++) {
        X1[i] = MIN;
        Y1[i] = MIN;
        X2[i] = MIN;
        Y2[i] = MIN;
    }

    S1[score_n2] = 0;
    X1[score_n2] = MIN;
    Y1[score_n2] = MIN;

    M[score_n2] = 0;
    S2[score_n2] = 0;
    X2[score_n2] = MIN;
    Y2[score_n2] = MIN;

    S1[2*score_n2-1] = 0;
    X1[2*score_n2-1] = MIN;
    Y1[2*score_n2-1] = MIN;

    M[2*score_n2-1] = 0;
    S2[2*score_n2-1] = 0;
    X2[2*score_n2-1] = MIN;
    Y2[2*score_n2-1] = MIN;

    for (long n=1; n<=n_stop; n++) {
        char c = toupper(seq1[n-1]);
        //char c = seq1[n-1];
        short* scs = scr->getScores(c);
        for (long m=m_start; m<=m_stop; m++) {

            long n_m = m + score_n2;
            long n_minus_1_m = m;
            long n_m_minus_1 = n_m - 1;
            long n_minus_1_m_minus_1 = (n_m - 1) - score_n2;

            int h = scr->getGapOpen();
            int g = scr->getGapExtend();
            char d = toupper(seq2[m-1]);
            int match_n_m = scs[d];

            S1[n_m] = Max(S1[n_minus_1_m_minus_1] + match_n_m,
                          X1[n_minus_1_m_minus_1] + match_n_m,
                          Y1[n_minus_1_m_minus_1] + match_n_m,
                          0);

            X1[n_m] = Max(S1[n_m_minus_1] + h + g,
                          X1[n_m_minus_1] + g);

            Y1[n_m] = Max(S1[n_minus_1_m] + h + g,
                          Y1[n_minus_1_m] + g);

            if (S1[n_m] > left.score) {
                left.score = S1[n_m];
                left.position = score_n2 * n + m;
            }

            M[n_m] = Max(M[n_minus_1_m], M[n_m_minus_1], S1[n_m]);

            S2[n_m] =  Max(S2[n_minus_1_m_minus_1] + match_n_m,
                          X2[n_minus_1_m_minus_1] + match_n_m,
                          Y2[n_minus_1_m_minus_1] + match_n_m,
                          M[n_minus_1_m_minus_1] + match_n_m,
                          0);

            X2[n_m] = Max(S2[n_m_minus_1] + h + g,
                          X2[n_m_minus_1] + g);

            Y2[n_m] = Max(S2[n_minus_1_m] + h + g,
                          Y2[n_minus_1_m] + g);

            if (S2[n_m] == S2[n_minus_1_m_minus_1] + match_n_m && R_S[n_minus_1_m_minus_1] != -1) {
                R_S[n_m] =  R_S[n_minus_1_m_minus_1];
            }
            if (S2[n_m] == X2[n_minus_1_m_minus_1] + match_n_m && R_X[n_minus_1_m_minus_1] != -1) {
                R_S[n_m] =  R_X[n_minus_1_m_minus_1];
            }
            if (S2[n_m] == Y2[n_minus_1_m_minus_1] + match_n_m && R_Y[n_minus_1_m_minus_1] != -1) {
                R_S[n_m] =  R_Y[n_minus_1_m_minus_1];
            }
            if (X2[n_m] == S2[n_m_minus_1] + h + g && R_S[n_m_minus_1] != -1) {
                R_X[n_m] = R_S[n_m_minus_1];
            }
            if (X2[n_m] == X2[n_m_minus_1] + g && R_X[n_m_minus_1] != -1) {
                R_X[n_m] = R_X[n_m_minus_1];
            }
            if (Y2[n_m] == S2[n_minus_1_m] + h + g && R_S[n_minus_1_m] != -1) {
                R_Y[n_m] = R_S[n_minus_1_m];
            }
            if (Y2[n_m] == Y2[n_minus_1_m] + g && R_Y[n_minus_1_m] != -1) {
                R_Y[n_m] = R_Y[n_minus_1_m];
            }

            if (S2[n_m] == M[n_minus_1_m_minus_1] + match_n_m) {
                R_S[n_m] = score_n2 * (n-1) + (m-1);
            }

            if (S2[n_m] >= right.score) {
                right.score = S2[n_m];
                right.position = score_n2 * n + m;
                last_m = R_S[n_m];
            }
        }
        for (long m=1; m<=len2; m++) {
            S1[m] = S1[m + score_n2];
            X1[m] = X1[m + score_n2];
            Y1[m] = Y1[m + score_n2];
            M[m] = M[m + score_n2];
            S2[m] = S2[m + score_n2];
            X2[m] = X2[m + score_n2];
            Y2[m] = Y2[m + score_n2];
            R_S[m] = R_S[m + score_n2];
            R_X[m] = R_X[m + score_n2];
            R_Y[m] = R_Y[m + score_n2];
        }
        S1[score_n2] = 0;
        X1[score_n2] = MIN;
        Y1[score_n2] = MIN;

        S2[score_n2] = 0;
        X2[score_n2] = MIN;
        Y2[score_n2] = MIN;

        S1[2*score_n2-1] = 0;
        X1[2*score_n2-1] = MIN;
        Y1[2*score_n2-1] = MIN;

        S2[2*score_n2-1] = 0;
        X2[2*score_n2-1] = MIN;
        Y2[2*score_n2-1] = MIN;
    }

    return make_pair(left, make_pair(right, last_m));
}

ScorePosition LongAGE::get_fisrt_BP(ScorePosition end) {

    ScorePosition st;
    bool exit = false;

    int S[2*score_n2];
    int X[2*score_n2];
    int Y[2*score_n2];

    for (long i=0; i<2*score_n2; i++) {
        S[i] = MIN;
        X[i] = MIN;
        Y[i] = MIN;
    }

    long total_neg = 6*score_n2;
    total_neg--;

    long end_n = split(end.position).n;
    long end_m = split(end.position).m;

    S[end_m] = scr->getMatch();

    for (long n = end_n; n > 0 && !exit; n--) {
        char c = toupper(seq1[n - 1]);
        //char c = seq1[n - 1];
        short* scs = scr->getScores(c);
        for (long m = end_m; m > 0; m--) {

            if (n * score_n2 + m == end.position) {
                continue;
            }

            long n_m = m;
            long n_plus_1_m = m + score_n2;
            long n_m_plus_1 = n_m + 1;
            long n_plus_1_m_plus_1 = m + score_n2 + 1;

            int h = scr->getGapOpen();
            int g = scr->getGapExtend();
            char d = toupper(seq2[m-1]);
            int match_n_m = scs[d];

            bool prev_s_inf = S[n_m] == MIN;
            bool prev_x_inf = X[n_m] == MIN;
            bool prev_y_inf = Y[n_m] == MIN;

            S[n_m] = Max(S[n_plus_1_m_plus_1] + match_n_m,
                         X[n_plus_1_m_plus_1] + match_n_m,
                         Y[n_plus_1_m_plus_1] + match_n_m,
                         0);

            X[n_m] = Max(S[n_m_plus_1] + h + g,
                         X[n_m_plus_1] + g);

            Y[n_m] = Max(S[n_plus_1_m] + h + g,
                         Y[n_plus_1_m] + g);

            if (S[n_m] <= 0) {
                S[n_m] = MIN;
            }
            if (X[n_m] <= 0) {
                X[n_m] = MIN;
            }
            if (Y[n_m] <= 0) {
                Y[n_m] = MIN;
            }
            if (S[n_m] == MIN && !prev_s_inf) {
                total_neg++;
            }
            if (X[n_m] == MIN && !prev_x_inf) {
                total_neg++;
            }
            if (Y[n_m] == MIN && !prev_y_inf) {
                total_neg++;
            }
            if (S[n_m] > 0 && prev_s_inf) {
                total_neg--;
            }
            if (X[n_m] > 0 && prev_x_inf) {
                total_neg--;
            }
            if (Y[n_m] > 0 && prev_y_inf) {
                total_neg--;
            }

            if (S[n_m] > st.score) {
                st.score = S[n_m];
                st.position = n * score_n2 + m;
            }
        }
        if (total_neg == 6*score_n2) {
            exit = true;
            break;
        }
        for (long m=score_n2+1; m<2*score_n2-1 && !exit; m++) {
            S[m] = S[m - score_n2];
            X[m] = X[m - score_n2];
            Y[m] = Y[m - score_n2];
        }
    }
    return st;
}
