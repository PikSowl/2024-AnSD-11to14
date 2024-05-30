#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstring>
#include <unordered_map>

#define intLimit 2147483647
#define matrixSize 10

using std::cin;
using std::cout;
using std::string;
using std::memset;
using std::swap;
using std::endl;
using std::vector;
using std::unordered_map;

void textInput(string& text){
    string temp;
    std::ifstream ifs("Babylon1.txt");
    while (!ifs.eof())
    {
        if (ifs.eof()) break;
        ifs >> temp;
        text += temp;
    }
    ifs.close();
}



//
///*
// 11 Aho-Corasick algorithm

const int abc_Size=28;

struct bohr_vrtx{
    int next_vrtx[abc_Size],pat_num,suff_link,auto_move[abc_Size],par,suff_flink;
    bool flag;
    char symb;
};

vector<bohr_vrtx> bohr;
vector<string> pattern;
bohr_vrtx make_bohr_vrtx(int p,char c){
    bohr_vrtx v;
    memset(v.next_vrtx, 255, sizeof(v.next_vrtx));
    memset(v.auto_move, 255, sizeof(v.auto_move));
    v.flag=false;
    v.suff_link=-1;
    v.par=p;
    v.symb=c;
    v.suff_flink=-1;
    return v;
}

void add_string_to_bohr(const string& s){
    int num=0;
    for (int i=0; i<s.length(); i++){
        char ch=s[i]-'a';
        if (bohr[num].next_vrtx[ch]==-1){
            bohr.push_back(make_bohr_vrtx(num,ch));
            bohr[num].next_vrtx[ch]=bohr.size()-1;
        }
        num=bohr[num].next_vrtx[ch];
    }
    bohr[num].flag=true;
    pattern.push_back(s);
    bohr[num].pat_num=pattern.size()-1;
}

bool is_string_in_bohr(const string& s){
    int num=0;
    for (int i=0; i<s.length(); i++){
        char ch=s[i]-'a';
        if (bohr[num].next_vrtx[ch]==-1){
            return false;
        }
        num=bohr[num].next_vrtx[ch];
    }
    return true;
}

int get_auto_move(int v, char ch);

int get_suff_link(int v){
    if (bohr[v].suff_link==-1)
        if (v==0||bohr[v].par==0)
            bohr[v].suff_link=0;
        else
            bohr[v].suff_link=get_auto_move(get_suff_link(bohr[v].par), bohr[v].symb);
    return bohr[v].suff_link;
}

int get_auto_move(int v, char ch){
    if (bohr[v].auto_move[ch]==-1)
        if (bohr[v].next_vrtx[ch]!=-1)
            bohr[v].auto_move[ch]=bohr[v].next_vrtx[ch];
        else
        if (v==0)
            bohr[v].auto_move[ch]=0;
        else
            bohr[v].auto_move[ch]=get_auto_move(get_suff_link(v), ch);
    return bohr[v].auto_move[ch];
}

int get_suff_flink(int v){
    if (bohr[v].suff_flink==-1){
        int u=get_suff_link(v);
        if (u==0)
            bohr[v].suff_flink=0;
        else
            bohr[v].suff_flink=(bohr[u].flag)?u:get_suff_flink(u);
    }
    return bohr[v].suff_flink;
}

void check(int v,int i){
    for(int u=v;u!=0;u=get_suff_flink(u)){
        if (bohr[u].flag) cout << i-pattern[bohr[u].pat_num].length()+1 << ':' << pattern[bohr[u].pat_num] << endl;
    }
}

void find_all_pos(const string& s){
    int u=0;
    for(int i=0;i<s.length();i++){
        u=get_auto_move(u,s[i]-'a');
        check(u,i+1);
    }
}

void lab_11() { //Aho–Corasick
    bohr.push_back(make_bohr_vrtx(0,'$'));
    add_string_to_bohr("eyes");
    add_string_to_bohr("farm");
    add_string_to_bohr("ores");
    add_string_to_bohr("blue");
    add_string_to_bohr("nail");

    string text;
    textInput(text);
    find_all_pos(text);
}

//
//*/
///*
// 12 Knuth–Morris–Pratt algorithm

vector<int> prefix_function(const string& s) {
    vector<int> pi(s.length(), 0);
    for (int i = 1; i < s.length(); i++) {
        int j = pi[i - 1];

        while (j > 0 && s[i] != s[j])
            j = pi[j - 1];

        if (s[i] == s[j])  pi[i] = j + 1;
        else pi[i] = j;
    }
    return pi;
}

void lab_12(){
    string search_for, text;
    cout << "Input line to search for:";
    cin >> search_for;
    textInput(text);

    vector<int> pi = prefix_function(search_for + '#' + text);

    int sf_length = search_for.length();

    for (int i = 0; i < text.length(); i++) {
        if (pi[sf_length + 1 + i] == sf_length) {
            cout << i - sf_length + 1 << ':' << search_for << endl;
        }
    }
}

//
//*/
///*
// 13 Boyer–Moore algorithm
vector<int> prefix_func(const string &s) {
    vector<int> p(s.length());

    int k = 0;
    p[0] = 0;
    for (int i = 1; i < s.length(); ++i) {
        while (k > 0 && s[k] != s[i]) {
            k = p[k - 1];
        }
        if (s[k] == s[i]) {
            ++k;
        }
        p[i] = k;
    }
    return p;
}

int find(string &text, string &search_for) {
    if (text.length() < search_for.length()) {
        return -1;
    }

    typedef unordered_map<char, int> TStopTable;
    typedef unordered_map<int, int> TSufficsTable;
    TStopTable stop_table;
    TSufficsTable suffics_table;

    for (int i = 0; i < search_for.length(); ++i) {
        stop_table[search_for[i]] = i;
    }

    string rt(search_for.rbegin(), search_for.rend());
    vector<int> p = prefix_func(search_for), pr = prefix_func(rt);
    for (int i = 0; i < search_for.length() + 1; ++i) {
        suffics_table[i] = search_for.length() - p.back();
    }

    for (int i = 1; i < search_for.length(); ++i) {
        int j = pr[i];
        suffics_table[j] = std::min(suffics_table[j], i - pr[i] + 1);
    }

    for (int shift = 0; shift <= text.length() - search_for.length();) {
        int pos = search_for.length() - 1;

        while (search_for[pos] == text[pos + shift]) {
            if (pos == 0) return shift;
            --pos;
        }

        if (pos == search_for.length() - 1) {
            TStopTable::const_iterator stop_symbol = stop_table.find(text[pos + shift]);
            int stop_symbol_additional = pos - (stop_symbol != stop_table.end() ? stop_symbol->second : -1);
            shift += stop_symbol_additional;
        } else {
            shift += suffics_table[search_for.length() - pos - 1];
        }
    }
    return -1;
}

void lab_13(){
    string search_for, text;
    cout << "Input line to search for:";
    cin >> search_for;
    textInput(text);
    cout << find(text, search_for) << ':' << search_for << endl;
}
//
//*/
///*
// 14 Rabin–Karp algorithm
#define d 10

void rabinKarp(string& text, string& search_for, int q) {
    int m = search_for.length();
    int n = text.length();
    int i, j;
    int p = 0;
    int t = 0;
    int h = 1;

    for (i = 0; i < m - 1; i++)
        h = (h * d) % q;

    for (i = 0; i < m; i++) {
        p = (d * p + search_for[i]) % q;
        t = (d * t + text[i]) % q;
    }

    for (i = 0; i <= n - m; i++) {
        if (p == t) {
            for (j = 0; j < m; j++) {
                if (text[i + j] != search_for[j])
                    break;
            }
            if (j == m)
                cout << i << ':' << search_for << endl;
        }
        if (i < n - m) {
            t = (d * (t - text[i] * h) + text[i + m]) % q;
            if (t < 0)
                t = (t + q);
        }
    }
}

void lab_14(){
    string text, search_for;
    cout << "Input line to search for:";
    cin >> search_for;
    textInput(text);
    rabinKarp(text,search_for, 26);
}

int main() {
    int lab_num;
    cout << "Chose lab number (options: 11, 12, 13, 14):";
    cin >> lab_num;
    while (true) {
        switch (lab_num) {
            case 11:
                lab_11();
                return 0;
            case 12:
                lab_12();
                return 0;
            case 13:
                lab_13();
                return 0;
            case 14:
                lab_14();
                return 0;
            default:
                cout << "wrong input" << endl << "Chose lab number:";
                cin >> lab_num;
                break;
        }
    }
}

