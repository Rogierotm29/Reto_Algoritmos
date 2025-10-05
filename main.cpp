#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <limits>
#include <unordered_map>
#include <array>

using namespace std;

// =============== helpers de texto ===============
static inline string trim(const string& s){
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a == string::npos) return "";
    return s.substr(a, b - a + 1);
}
static inline char toU(char c){ return (char)toupper((unsigned char)c); }
static inline string only_AA(const string& s){
    string t; t.reserve(s.size());
    for(char c: s){
        char u = toU(c);
        if ((u>='A' && u<='Z') || u=='*') t.push_back(u);
    }
    return t;
}
static inline string upper_copy(string s){
    for(char& c: s) c = toU(c);
    return s;
}
static inline bool is_digits(const string& s){
    if (s.empty()) return false;
    for(char c: s) if(!isdigit((unsigned char)c)) return false;
    return true;
}

// =============== estructuras ===============
struct ProteinItem { string name; string aa; };

struct Data {
    string genome_wuhan;
    string genome_texas;
    string gene_M;
    string gene_ORF1ab;
    string gene_S;
    vector<ProteinItem> proteins;
} DB;

// =============== paths/utilidades ===============
static string path_join(const string& base, const string& file){
    std::filesystem::path p(base); p /= file; return p.string();
}
static string existing_file(const string& base, const string& nameNoExt){
    string p1 = path_join(base, nameNoExt);
    if (std::filesystem::exists(p1)) return p1;
    string p2 = path_join(base, nameNoExt + ".txt");
    if (std::filesystem::exists(p2)) return p2;
    return p1; // deja que el lector lance error informativo
}

// =============== lectores FASTA ===============
string read_fasta_dna_clean(const string& path){
    ifstream in(path);
    if(!in) throw runtime_error("No pude abrir: " + path);
    string line, seq; seq.reserve(400000);
    while(getline(in, line)){
        if(!line.empty() && line[0]=='>') continue;
        for(char c: line){
            char u = toU(c);
            if(u=='A'||u=='C'||u=='G'||u=='T') seq.push_back(u);
        }
    }
    return seq;
}
vector<ProteinItem> read_protein_fasta(const string& path){
    ifstream in(path);
    if(!in) throw runtime_error("No pude abrir: " + path);
    vector<ProteinItem> items;
    string line, current_name, current_seq;
    auto push_current = [&](){
        if(!current_name.empty() && !current_seq.empty()){
            items.push_back({current_name, current_seq});
        }
        current_name.clear(); current_seq.clear();
    };
    while(getline(in, line)){
        string s = trim(line);
        if(s.empty()) continue;
        if(s[0]=='>'){
            size_t i=0; while(i<s.size() && s[i]=='>') ++i;
            string header = trim(s.substr(i));
            if(header.empty()){
                header = "unnamed_" + to_string(items.size()+1);
            }else{
                size_t sp = header.find_first_of(" \t");
                if(sp != string::npos) header = header.substr(0, sp);
            }
            push_current();
            current_name = header;
        }else{
            for(char c: s){
                char u = toU(c);
                if((u>='A' && u<='Z') || u=='*') current_seq.push_back(u);
            }
        }
    }
    push_current();
    return items;
}

// =============== carga completa ===============
void load_all(const string& BASE){
    DB.genome_wuhan = read_fasta_dna_clean(existing_file(BASE, "SARS-COV-2-MN908947.3"));
    DB.genome_texas = read_fasta_dna_clean(existing_file(BASE, "SARS-COV-2-MT106054.1"));
    DB.gene_M       = read_fasta_dna_clean(existing_file(BASE, "gen-M"));
    DB.gene_ORF1ab  = read_fasta_dna_clean(existing_file(BASE, "gen-ORF1AB"));
    DB.gene_S       = read_fasta_dna_clean(existing_file(BASE, "gen-S"));
    DB.proteins     = read_protein_fasta(existing_file(BASE, "seq-proteins"));
}

//============= Búsqueda Naive (propia) =============
inline bool match_at (const string& text, size_t i, const string& pattern){
    if (i + pattern.size() > text.size()) return false;
    for(size_t k=0; k<pattern.size(); ++k){
        if(text[i + k] != pattern[k]) return false;
    }
    return true;
}
vector<size_t> naive_search_all(const string& text, const string& pattern){
    vector<size_t> pos;
    if(pattern.empty() || text.empty() || pattern.size() > text.size()) return pos;
    for(size_t i=0; i + pattern.size() <= text.size(); ++i){
        if(match_at(text, i, pattern)) pos.push_back(i);
    }
    return pos;
}

//============= ADN básico (comp, revcomp) ==========
inline char comp(char b){
    switch(b){
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        default:  return 'N';
    }
}
string revcomp(const string& s){
    string rc; rc.resize(s.size());
    for(size_t i=0, n=s.size(); i<n; ++i) rc[n-1-i] = comp(s[i]);
    return rc;
}

//============= Palíndromo rev-comp (centros pares) =========
pair<size_t,size_t> longest_revcomp_palindrome(const string& s){
    auto eq_comp = [&](char a, char b){ return comp(a) == b; };
    size_t bestL = 0, bestR = 0;
    if(s.size() < 2) return {0,0};
    for(size_t i=0; i+1<s.size(); ++i){
        long long L = (long long)i, R = (long long)i + 1;
        while(L>=0 && (size_t)R<s.size() && eq_comp(s[(size_t)L], s[(size_t)R])){
            if((size_t)(R-L) > (bestR-bestL)){ bestL=(size_t)L; bestR=(size_t)R; }
            --L; ++R;
        }
    }
    return {bestL, bestR};
}

// Palíndromo "directo" (
pair<size_t,size_t> longest_direct_palindrome(const string& s){
    if (s.empty()) return {0,0};
    size_t bestL = 0, bestR = 0;

    auto expand = [&](long long L, long long R){
        while (L >= 0 && (size_t)R < s.size() && s[(size_t)L] == s[(size_t)R]){
            if ((size_t)(R - L) > (bestR - bestL)){ bestL = (size_t)L; bestR = (size_t)R; }
            --L; ++R;
        }
    };

    // centros impares y pares
    for (size_t i = 0; i < s.size(); ++i){
        expand((long long)i, (long long)i);         // impar
        if (i + 1 < s.size()) expand((long long)i, (long long)i + 1); // par
    }
    return {bestL, bestR};
}


//============= Código genético y traducción =========
static const unordered_map<string, char> GENETIC_CODE = {
    {"TTT",'F'},{"TTC",'F'},{"TTA",'L'},{"TTG",'L'},
    {"TCT",'S'},{"TCC",'S'},{"TCA",'S'},{"TCG",'S'},
    {"TAT",'Y'},{"TAC",'Y'},{"TAA",'*'},{"TAG",'*'},
    {"TGT",'C'},{"TGC",'C'},{"TGA",'*'},{"TGG",'W'},
    {"CTT",'L'},{"CTC",'L'},{"CTA",'L'},{"CTG",'L'},
    {"CCT",'P'},{"CCC",'P'},{"CCA",'P'},{"CCG",'P'},
    {"CAT",'H'},{"CAC",'H'},{"CAA",'Q'},{"CAG",'Q'},
    {"CGT",'R'},{"CGC",'R'},{"CGA",'R'},{"CGG",'R'},
    {"ATT",'I'},{"ATC",'I'},{"ATA",'I'},{"ATG",'M'},
    {"ACT",'T'},{"ACC",'T'},{"ACA",'T'},{"ACG",'T'},
    {"AAT",'N'},{"AAC",'N'},{"AAA",'K'},{"AAG",'K'},
    {"AGT",'S'},{"AGC",'S'},{"AGA",'R'},{"AGG",'R'},
    {"GTT",'V'},{"GTC",'V'},{"GTA",'V'},{"GTG",'V'},
    {"GCT",'A'},{"GCC",'A'},{"GCA",'A'},{"GCG",'A'},
    {"GAT",'D'},{"GAC",'D'},{"GAA",'E'},{"GAG",'E'},
    {"GGT",'G'},{"GGC",'G'},{"GGA",'G'},{"GGG",'G'}
};
string translate_frame(const string& dna, size_t frame, bool until_stop=false){
    string aa;
    if (frame>2 || dna.size() < frame+3) return aa;
    aa.reserve(dna.size()/3);
    for(size_t i=frame; i+2<dna.size(); i+=3){
        string codon; codon.push_back(dna[i]); codon.push_back(dna[i+1]); codon.push_back(dna[i+2]);
        auto it = GENETIC_CODE.find(codon);
        char c = (it==GENETIC_CODE.end()) ? 'X' : it->second;
        if (until_stop && c=='*') break;
        aa.push_back(c);
    }
    return aa;
}
array<string,3> translate_3_frames_forward(const string& dna){
    return { translate_frame(dna,0,false),
             translate_frame(dna,1,false),
             translate_frame(dna,2,false) };
}
array<string,3> translate_3_frames_reverse(const string& dna){
    string rc = revcomp(dna);
    return { translate_frame(rc,0,false),
             translate_frame(rc,1,false),
             translate_frame(rc,2,false) };
}

//============= Buscar AA en el genoma (6 marcos) =========
struct HitAA{
    bool reverse;  // 0=fwd, 1=rev
    int frame;     // 0,1,2
    size_t aa_pos; // índice AA en la proteína traducida
    size_t nt_start, nt_end; // mapeo en nt
};
vector<HitAA> find_protein_in_genome(const string& genome, const string& aa_pattern){
    vector<HitAA> hits;
    if (aa_pattern.empty()) return hits;

    // forward
    auto fwd = translate_3_frames_forward(genome);
    for(int f=0; f<3; ++f){
        auto occ = naive_search_all(fwd[(size_t)f], aa_pattern);
        for(size_t p: occ){
            size_t nt_start = p*3 + (size_t)f;
            size_t nt_end   = nt_start + aa_pattern.size()*3 - 1;
            if (nt_end < genome.size()){
                hits.push_back({false, f, p, nt_start, nt_end});
            }
        }
    }
    // reverse (reusa una vez el rc)
    string rc = revcomp(genome);
    auto rev = array<string,3>{
        translate_frame(rc,0,false),
        translate_frame(rc,1,false),
        translate_frame(rc,2,false)
    };
    for(int f=0; f<3; ++f){
        auto occ = naive_search_all(rev[(size_t)f], aa_pattern);
        for(size_t p: occ){
            size_t start_rc = p*3 + (size_t)f;
            size_t end_rc   = start_rc + aa_pattern.size()*3 - 1;
            if (end_rc >= rc.size()) continue;
            size_t N = genome.size();
            size_t orig_start = N - 1 - end_rc;
            size_t orig_end   = N - 1 - start_rc;
            hits.push_back({true, f, p, orig_start, orig_end});
        }
    }
    sort(hits.begin(), hits.end(), [](const HitAA& a, const HitAA& b){
        return a.nt_start < b.nt_start;
    });
    return hits;
}

//============= Diff Wuhan vs Texas (SNP/INDEL greedy) ========
struct DiffEvent {
    enum Kind { SNP, DEL, INS } kind;
    size_t posA;   // índice en Wuhan
    size_t posB;   // índice en Texas
    string a;      // en Wuhan (1 base para SNP, k para DEL)
    string b;      // en Texas (1 base para SNP, k para INS)
};
static inline bool kmatch(const string& A, size_t ai, const string& B, size_t bj, size_t k){
    if (ai + k > A.size() || bj + k > B.size()) return false;
    for (size_t t=0; t<k; ++t) if (A[ai+t] != B[bj+t]) return false;
    return true;
}
static vector<DiffEvent> compare_genomes_anchored(const string& A, const string& B,
                                                  size_t k = 12, size_t L = 200){
    vector<DiffEvent> ev;
    size_t i = 0, j = 0;
    while (i < A.size() && j < B.size()){
        if (A[i] == B[j]) { ++i; ++j; continue; }
        if (i+1 < A.size() && j+1 < B.size() && A[i+1] == B[j+1]){
            ev.push_back({DiffEvent::SNP, i, j, string(1, A[i]), string(1, B[j])});
            ++i; ++j; continue;
        }
        size_t best_del = SIZE_MAX, best_ins = SIZE_MAX;
        for (size_t a = 1; a <= L && i + a + k <= A.size(); ++a){
            if (kmatch(A, i + a, B, j, k)){ best_del = a; break; }
        }
        for (size_t b = 1; b <= L && j + b + k <= B.size(); ++b){
            if (kmatch(A, i, B, j + b, k)){ best_ins = b; break; }
        }
        if (best_del == SIZE_MAX && best_ins == SIZE_MAX){
            ev.push_back({DiffEvent::SNP, i, j, string(1, A[i]), string(1, B[j])});
            ++i; ++j;
        } else if (best_ins == SIZE_MAX || (best_del < best_ins)){
            ev.push_back({DiffEvent::DEL, i, j, A.substr(i, best_del), ""});
            i += best_del;
        } else {
            ev.push_back({DiffEvent::INS, i, j, "", B.substr(j, best_ins)});
            j += best_ins;
        }
    }
    if (i < A.size()) ev.push_back({DiffEvent::DEL, i, j, A.substr(i), ""});
    if (j < B.size()) ev.push_back({DiffEvent::INS, i, j, "", B.substr(j)});
    return ev;
}

//============= Impacto de SNP en AA + frameshift ========
struct AAImpact {
    string gene;
    size_t nt_pos;      // posición nt en Wuhan
    size_t codon_index; // índice de codón relativo al inicio del gen
    string ref_codon, alt_codon;
    char ref_aa, alt_aa;
};
static inline char aa_from_codon(const string& codon){
    auto it = GENETIC_CODE.find(codon);
    return (it==GENETIC_CODE.end()) ? 'X' : it->second;
}
vector<AAImpact> aa_impacts_for_gene(const string& wuhan,
                                     const vector<DiffEvent>& ev,
                                     const string& gene_seq,
                                     const string& gene_name){
    vector<AAImpact> out;
    auto occ = naive_search_all(wuhan, gene_seq);
    if (occ.empty()) return out;
    size_t s = occ[0]; // inicio del gen en Wuhan
    string mutated = gene_seq;

    // aplica SOLO SNPs dentro del gen (INDEL se trata aparte)
    for (const auto& e : ev){
        if (e.kind != DiffEvent::SNP) continue;
        if (e.posA < s || e.posA >= s + gene_seq.size()) continue;
        mutated[e.posA - s] = e.b[0];
    }
    size_t codons = min(gene_seq.size(), mutated.size()) / 3;
    for (size_t i = 0; i < codons; ++i){
        string c_ref = gene_seq.substr(3*i,3);
        string c_alt = mutated.substr(3*i,3);
        if (c_ref == c_alt) continue;
        char a_ref = aa_from_codon(c_ref);
        char a_alt = aa_from_codon(c_alt);
        if (a_ref != a_alt){
            out.push_back({gene_name, s + 3*i, i, c_ref, c_alt, a_ref, a_alt});
        }
    }
    return out;
}
void warn_frameshift_in_gene(const vector<DiffEvent>& ev,
                             size_t gene_start, size_t gene_end,
                             const string& gname){
    for (const auto& d : ev){
        if (d.posA < gene_start || d.posA > gene_end) continue;
        if (d.kind == DiffEvent::DEL){
            if (d.a.size() % 3 != 0)
                cout << "   [Frameshift] " << gname << " DEL no múltiplo de 3 en nt@" << d.posA << "\n";
        } else if (d.kind == DiffEvent::INS){
            if (d.b.size() % 3 != 0)
                cout << "   [Frameshift] " << gname << " INS no múltiplo de 3 en nt@" << d.posA << "\n";
        }
    }
}

//=============== Menús ===============
void menu_find_gene(){
    cout << "\n[1] Índices de M/ORF1ab/S en Wuhan y Texas (primeros 12 nt)\n";

    auto print_for = [&](const string& label, const string& genome){
        cout << "== " << label << " ==\n";
        if (genome.empty()){
            cout << "   Error: genoma " << label << " no está cargado.\n";
            return;
        }

        struct GeneRow { string name; const string* seq; };
        vector<GeneRow> genes = {
            {"M",      &DB.gene_M},
            {"ORF1ab", &DB.gene_ORF1ab},
            {"S",      &DB.gene_S}
        };

        for (const auto& g : genes){
            if (g.seq->empty()){
                cout << "Gen " << g.name << ": secuencia del gen vacía.\n";
                continue;
            }
            auto occ = naive_search_all(genome, *g.seq);
            if (occ.empty()){
                cout << "Gen " << g.name << ": no encontrado en " << label << ".\n";
                continue;
            }
            size_t idx = occ[0];
            string first12 = (idx + 12 <= genome.size()) ? genome.substr(idx, 12) : "";
            cout << "Gen " << g.name << ":\n";
            cout << "  Índice: " << idx
                 << ", Primeros 12 nt: " << first12 << "\n";
        }
    };

    print_for("Wuhan 2019", DB.genome_wuhan);
    print_for("Texas 2020", DB.genome_texas);
}


void menu_longest_palindrome(){
    cout << "\n[2] Palíndromo más largo en un gen (elige rev-comp o directo)\n";
    cout << "   Gen: 1) M  2) ORF1ab  3) S  -> " << flush;
    int op; 
    if(!(cin >> op)){ cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n'); cout<<"   Entrada inválida.\n"; return; }

    const string* gene = nullptr; string gname;
    if (op==1){ gene=&DB.gene_M; gname="M"; }
    else if (op==2){ gene=&DB.gene_ORF1ab; gname="ORF1ab"; }
    else if (op==3){ gene=&DB.gene_S; gname="S"; }
    else { cout << "   Opción inválida.\n"; return; }

    cout << "   Tipo: 1) Rev-Comp  2) Directo  -> " << flush;
    int tipo;
    if(!(cin >> tipo)){ cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n'); cout<<"   Entrada inválida.\n"; return; }

    pair<size_t,size_t> LR;
    string label;

    if (tipo == 2){
        label = "palíndromo DIRECTO";
        LR = longest_direct_palindrome(*gene);
        size_t L = LR.first, R = LR.second;
        if (R <= L){ cout << "   No se detectó palíndromo directo no trivial.\n"; return; }
        if (R - L + 1 < 2){ cout << "   No se detectó palíndromo directo no trivial.\n"; return; }
        cout << "   Gen " << gname << " → mejor " << label << "\n";
        cout << "   L=" << L << "  R=" << R << "  (longitud=" << (R - L + 1) << ")\n";
        string pal = gene->substr(L, R - L + 1);
        for(size_t i=0; i<pal.size(); i+=60)
            cout << "   " << pal.substr(i, min<size_t>(60, pal.size()-i)) << "\n";
    } else {
        label = "palíndromo REV-COMP";
        LR = longest_revcomp_palindrome(*gene);
        size_t L = LR.first, R = LR.second;
        if (R <= L){ cout << "   No se detectó palíndromo rev-comp no trivial.\n"; return; }
        cout << "   Gen " << gname << " → mejor " << label << "\n";
        cout << "   L=" << L << "  R=" << R << "  (longitud=" << (R - L + 1) << ")\n";
        string pal = gene->substr(L, R - L + 1);
        for(size_t i=0; i<pal.size(); i+=60)
            cout << "   " << pal.substr(i, min<size_t>(60, pal.size()-i)) << "\n";
    }
}


void menu_translate_gene(){
    cout << "\n[3] Traducir un gen a proteína (6 marcos)\n";
    cout << "   Gen: 1) M  2) ORF1ab  3) S  -> " << flush;
    int op; if(!(cin >> op)){ cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n'); cout<<"   Entrada inválida.\n"; return; }
    const string* gene = nullptr; string gname;
    if (op==1){ gene=&DB.gene_M; gname="M"; }
    else if (op==2){ gene=&DB.gene_ORF1ab; gname="ORF1ab"; }
    else if (op==3){ gene=&DB.gene_S; gname="S"; }
    else { cout << "   Opción inválida.\n"; return; }

    auto fwd = translate_3_frames_forward(*gene);
    auto rev = translate_3_frames_reverse(*gene);
    auto printAA = [&](const string& label, const string& aa){
        cout << "   " << label << " | len=" << aa.size() << "\n";
        for(size_t i=0; i<aa.size(); i+=60)
            cout << "   " << aa.substr(i, min<size_t>(60, aa.size()-i)) << "\n";
    };
    cout << "   Proteínas candidatas por marco:\n";
    printAA("FW frame 0", fwd[0]);
    printAA("FW frame 1", fwd[1]);
    printAA("FW frame 2", fwd[2]);
    printAA("RV frame 0", rev[0]);
    printAA("RV frame 1", rev[1]);
    printAA("RV frame 2", rev[2]);
}

void menu_find_protein_in_genome(){
    cout << "\n[4] Buscar una proteína (AA) en el GENOMA (6 marcos)\n";
    int mode = 1;
    if (!DB.proteins.empty()){
        cout << "   ¿Cómo quieres ingresar la proteína?\n";
        cout << "   1) Pegar SECUENCIA de AA\n";
        cout << "   2) Elegir por NOMBRE desde seq-proteins (" << DB.proteins.size() << " disponibles)\n";
        cout << "   -> " << flush;
        if(!(cin >> mode)){
            cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "   Entrada inválida.\n"; return;
        }
    }
    string aa; string chosen_name;
    if (mode == 2 && !DB.proteins.empty()){
        cout << "   Mostrar primeros 20? (y/n): "; char yn; cin >> yn;
        if (yn=='y' || yn=='Y'){
            size_t lim = min<size_t>(20, DB.proteins.size());
            for(size_t i=0;i<lim;++i)
                cout << "   [" << i << "] " << DB.proteins[i].name
                     << "  (lenAA=" << DB.proteins[i].aa.size() << ")\n";
        }
        cout << "   Escribe un INDICE o parte del nombre: ";
        string token; cin >> token;
        bool ok = false;
        if (is_digits(token)){
            size_t idx = (size_t)stoull(token);
            if (idx < DB.proteins.size()){
                aa = DB.proteins[idx].aa; chosen_name = DB.proteins[idx].name; ok = true;
            }
        }
        if (!ok){
            string q = upper_copy(token);
            vector<size_t> matches;
            for(size_t i=0;i<DB.proteins.size();++i){
                string up = upper_copy(DB.proteins[i].name);
                if (up.find(q) != string::npos) matches.push_back(i);
            }
            if (matches.empty()){ cout << "   Nombre no encontrado.\n"; return; }
            if (matches.size() > 20){ cout << "   " << matches.size() << " coincidencias, muestro primeras 20:\n"; matches.resize(20); }
            for(size_t k=0;k<matches.size();++k){
                size_t i = matches[k];
                cout << "   [" << i << "] " << DB.proteins[i].name
                     << "  (lenAA=" << DB.proteins[i].aa.size() << ")\n";
            }
            cout << "   Elige INDICE exacto: ";
            string idxs; cin >> idxs;
            if (!is_digits(idxs)){ cout << "   Entrada inválida.\n"; return; }
            size_t idx = (size_t)stoull(idxs);
            if (idx >= DB.proteins.size()){ cout << "   Índice fuera de rango.\n"; return; }
            aa = DB.proteins[idx].aa; chosen_name = DB.proteins[idx].name;
        }
    } else {
        cout << "   Escribe la SECUENCIA de AA (A-Z, '*' si quieres):\n";
        string raw; cin >> raw;
        aa = only_AA(raw);
        if (aa.empty()){ cout << "   Secuencia vacía o inválida.\n"; return; }
        if (aa.size() < 6) cout << "   Aviso: secuencia muy corta, puede haber muchos falsos negativos.\n";
    }
    cout << "   Genoma: 1) Wuhan 2019  2) Texas 2020  -> " << flush;
    int gg; 
    if(!(cin >> gg)){ cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n'); cout<<"   Entrada inválida.\n"; return; }
    const string* genome = (gg==2) ? &DB.genome_texas : &DB.genome_wuhan;
    string glabel = (gg==2) ? "Texas 2020" : "Wuhan 2019";

    auto hits = find_protein_in_genome(*genome, aa);
    if (hits.empty()){
        cout << "   " << (chosen_name.empty()? "(secuencia pegada)" : "(" + chosen_name + ")")
             << " no se encontró en " << glabel << ".\n";
        return;
    }
    cout << "   " << (chosen_name.empty()? "Secuencia" : "("+chosen_name+")")
         << " hallada en " << glabel << ": " << hits.size() << " ocurrencias\n";
    string aa4 = (aa.size()>=4 ? aa.substr(0,4) : aa);
    for(const auto& h: hits){
        string nt12 = (h.nt_start+12<=genome->size()? genome->substr(h.nt_start,12) : "");
        cout << "   - " << (h.reverse ? "REV" : "FWD")
             << " frame " << h.frame
             << " | nt[" << h.nt_start << "," << h.nt_end << "]"
             << " | aa_pos=" << h.aa_pos
             << " | 4AA=" << aa4
             << " | nt12=" << nt12
             << "\n";
    }
}

void menu_find_protein_by_name(){
    if (DB.proteins.empty()){
        cout << "\n[5] Buscar proteína por nombre\n";
        cout << "   No cargué seq-proteins.txt, así que esta opción no tiene data.\n";
        return;
    }
    cout << "\n[5] Buscar proteína por NOMBRE en el genoma (6 marcos)\n";
    cout << "   Mostrar primeros 20? (y/n): ";
    char show; cin >> show;
    if (show=='y' || show=='Y'){
        size_t lim = min<size_t>(20, DB.proteins.size());
        for(size_t i=0;i<lim;++i)
            cout << "   [" << i << "] " << DB.proteins[i].name 
                 << "  (lenAA=" << DB.proteins[i].aa.size() << ")\n";
    }
    cout << "   Escribe INDICE o parte del nombre: ";
    string token; cin >> token;

    string aa; string chosen; bool ok=false;
    if (is_digits(token)){
        size_t idx = (size_t)stoull(token);
        if (idx < DB.proteins.size()){ aa=DB.proteins[idx].aa; chosen=DB.proteins[idx].name; ok=true; }
    }
    if(!ok){
        string q = upper_copy(token);
        vector<size_t> matches;
        for(size_t i=0;i<DB.proteins.size();++i){
            string up = upper_copy(DB.proteins[i].name);
            if (up.find(q) != string::npos) matches.push_back(i);
        }
        if(matches.empty()){ cout << "   Nombre no encontrado.\n"; return; }
        if(matches.size()>20){ cout << "   " << matches.size() << " coincidencias, muestro primeras 20:\n"; matches.resize(20); }
        for(size_t k=0;k<matches.size();++k){
            size_t i=matches[k];
            cout << "   [" << i << "] " << DB.proteins[i].name 
                 << "  (lenAA=" << DB.proteins[i].aa.size() << ")\n";
        }
        cout << "   Elige INDICE exacto: ";
        string idxs; cin >> idxs;
        if(!is_digits(idxs)){ cout << "   Entrada inválida.\n"; return; }
        size_t idx = (size_t)stoull(idxs);
        if(idx>=DB.proteins.size()){ cout << "   Índice fuera de rango.\n"; return; }
        aa = DB.proteins[idx].aa; chosen = DB.proteins[idx].name;
    }

    cout << "   Genoma: 1) Wuhan 2019  2) Texas 2020  -> " << flush;
    int gg; if(!(cin >> gg)){ cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n'); cout<<"   Entrada inválida.\n"; return; }
    const string* genome = (gg==2)? &DB.genome_texas : &DB.genome_wuhan;
    string glabel = (gg==2)? "Texas 2020" : "Wuhan 2019";

    auto hits = find_protein_in_genome(*genome, aa);
    if (hits.empty()){
        cout << "   (" << chosen << ") no se encontró en " << glabel << ".\n";
        return;
    }
    cout << "   (" << chosen << ") hallada en " << glabel << ": " << hits.size() << " ocurrencias\n";
    string aa4 = (aa.size()>=4 ? aa.substr(0,4) : aa);
    for (const auto& h : hits) {
        string nt12 = (h.nt_start+12<=genome->size()? genome->substr(h.nt_start,12) : "");
        cout << "   - " << (h.reverse ? "REV" : "FWD") 
             << " frame " << h.frame 
             << " | nt[" << h.nt_start << "," << h.nt_end << "]"
             << " | aa_pos=" << h.aa_pos
             << " | 4AA=" << aa4
             << " | nt12=" << nt12
             << "\n";
    }
}

void menu_compare_genomes_snps(){
    cout << "\n[6] Comparar Wuhan 2019 vs Texas 2020 (SNPs / INDELs)\n";
    const string& A = DB.genome_wuhan; // Wuhan
    const string& B = DB.genome_texas; // Texas
    if (A.empty() || B.empty()){ cout << "   Error: algún genoma no está cargado.\n"; return; }

    auto ev = compare_genomes_anchored(A, B, /*k=*/12, /*L=*/200);

    size_t snp = 0, del_events = 0, ins_events = 0, del_bases = 0, ins_bases = 0;
    for (const auto& e : ev){
        if (e.kind == DiffEvent::SNP) snp++;
        else if (e.kind == DiffEvent::DEL){ del_events++; del_bases += e.a.size(); }
        else if (e.kind == DiffEvent::INS){ ins_events++; ins_bases += e.b.size(); }
    }
    cout << "   Wuhan len=" << A.size() << " | Texas len=" << B.size() << "\n";
    cout << "   SNPs: " << snp
         << " | DEL eventos: " << del_events << " (bases=" << del_bases << ")"
         << " | INS eventos: " << ins_events << " (bases=" << ins_bases << ")\n";

    size_t show = min<size_t>(50, ev.size());
    if (show == 0){ cout << "   No se encontraron diferencias.\n"; return; }
    cout << "   Primeros " << show << " eventos:\n";
    for (size_t k = 0; k < show; ++k){
        const auto& e = ev[k];
        if (e.kind == DiffEvent::SNP){
            cout << "   SNP  @Wuhan:" << e.posA << " | " << e.a << " -> " << e.b << "\n";
        } else if (e.kind == DiffEvent::DEL){
            cout << "   DEL  @Wuhan:" << e.posA << " | '" << e.a << "' (len=" << e.a.size() << ")\n";
        } else {
            cout << "   INS  @Wuhan:" << e.posA << " | Texas tiene '" << e.b
                 << "' (len=" << e.b.size() << ")\n";
        }
    }

    // Impacto AA por SNP dentro de M / ORF1ab / S
    auto impM   = aa_impacts_for_gene(A, ev, DB.gene_M,      "M");
    auto impORF = aa_impacts_for_gene(A, ev, DB.gene_ORF1ab, "ORF1ab");
    auto impS   = aa_impacts_for_gene(A, ev, DB.gene_S,      "S");
    auto printImp = [&](const vector<AAImpact>& v){
        for (const auto& x: v){
            cout << "   " << x.gene
                 << " nt@" << x.nt_pos
                 << " codon#" << x.codon_index
                 << " " << x.ref_codon << "(" << x.ref_aa << ") -> "
                 << x.alt_codon << "(" << x.alt_aa << ")\n";
        }
    };
    size_t totalAA = impM.size() + impORF.size() + impS.size();
    cout << "   Cambios de aminoácido por SNP dentro de genes (aprox): "
         << totalAA << "\n";
    printImp(impM); printImp(impORF); printImp(impS);

    // Avisos frameshift por INDEL no múltiplo de 3
    auto occM   = naive_search_all(A, DB.gene_M);
    auto occORF = naive_search_all(A, DB.gene_ORF1ab);
    auto occS   = naive_search_all(A, DB.gene_S);
    if (!occM.empty())   warn_frameshift_in_gene(ev, occM[0],   occM[0] + DB.gene_M.size() - 1,      "M");
    if (!occORF.empty()) warn_frameshift_in_gene(ev, occORF[0], occORF[0]+DB.gene_ORF1ab.size()-1,   "ORF1ab");
    if (!occS.empty())   warn_frameshift_in_gene(ev, occS[0],   occS[0] + DB.gene_S.size() - 1,      "S");
}

// =============== MAIN ===============
constexpr const char* DEFAULT_BASE =
R"(C:\Users\rogie\OneDrive\Escritorio\Reto_Algoritmo)";

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(&cout);

    try{
        string BASE = (argc >= 2) ? string(argv[1]) : string(DEFAULT_BASE);
        cout << "Cargando archivos desde: " << BASE << "\n";
        load_all(BASE);

        cout << "\n[OK] CARGA COMPLETA\n";
        cout << "Genomas: Wuhan=" << DB.genome_wuhan.size()
             << " | Texas=" << DB.genome_texas.size() << " bases\n";
        cout << "Genes:   M=" << DB.gene_M.size()
             << " | ORF1ab=" << DB.gene_ORF1ab.size()
             << " | S=" << DB.gene_S.size() << " bases\n";
        cout << "Proteínas cargadas: " << DB.proteins.size() << "\n";
        for(size_t i=0;i<DB.proteins.size() && i<3;++i){
            cout << "  [" << i << "] " << DB.proteins[i].name
                 << " (lenAA=" << DB.proteins[i].aa.size() << ")\n";
            cout << "      " << DB.proteins[i].aa.substr(0, min<size_t>(60, DB.proteins[i].aa.size())) << "\n";
        }

        while(true){
            cout << "\n=========== MENU ===========\n";
            cout << "1) Índices de M/ORF1ab/S en Wuhan y Texas (primeros 12 nt)\n";
            cout << "2) Palíndromo (rev-comp o directo) más largo en un gen\n";
            cout << "3) Traducir un gen a proteína (6 marcos)\n";
            cout << "4) Buscar proteína (AA) en el genoma (6 marcos)\n";
            cout << "5) Buscar proteína por NOMBRE (seq-proteins.txt)\n";
            cout << "6) Comparar Wuhan vs Texas (SNPs / INDELs)\n";
            cout << "0) Salir\n";
            cout << "Elige opción: " << flush;

            int op;
            if(!(cin >> op)){
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
                cout << "   Entrada inválida.\n";
                continue;
            }
            if(op==0) return 0;
            else if(op==1) menu_find_gene();
            else if(op==2) menu_longest_palindrome();
            else if(op==3) menu_translate_gene();
            else if(op==4) menu_find_protein_in_genome();
            else if(op==5) menu_find_protein_by_name();
            else if(op==6) menu_compare_genomes_snps();
            else cout << "   Opción inválida.\n";
        }
    }catch(const exception& e){
        cerr << "[ERROR] " << e.what() << "\n";
        return 1;
    }
}
