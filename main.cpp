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

// --- helpers extra para paso 4/5 ---
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
struct ProteinItem {
    string name;   // p.ej. "Spike_S1"
    string aa;     // secuencia AA (A-Z y *)
};

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
    std::filesystem::path p(base);
    p /= file;
    return p.string();
}

// Si el archivo no existe como "nombre", probamos "nombre.txt".
static string existing_file(const string& base, const string& nameNoExt){
    string p1 = path_join(base, nameNoExt);
    if (std::filesystem::exists(p1)) return p1;
    string p2 = path_join(base, nameNoExt + ".txt");
    if (std::filesystem::exists(p2)) return p2;
    // Devolvemos p1 para que el lector lance un error informativo
    return p1;
}

// =============== lectores FASTA ===============
// DNA: ignora cabeceras '>' y limpia a A/C/G/T en mayúsculas
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

// Proteínas: múltiples entradas con '>' nombre y secuencia multi-línea
vector<ProteinItem> read_protein_fasta(const string& path){
    ifstream in(path);
    if(!in) throw runtime_error("No pude abrir: " + path);

    vector<ProteinItem> items;
    string line, current_name, current_seq;

    auto push_current = [&](){
        if(!current_name.empty() && !current_seq.empty()){
            items.push_back({current_name, current_seq});
        }
        current_name.clear();
        current_seq.clear();
    };

    while(getline(in, line)){
        string s = trim(line);
        if(s.empty()) continue;

        if(s[0]=='>'){
            // corta uno o más '>' y toma hasta el primer espacio como nombre
            size_t i = 0; while(i<s.size() && s[i]=='>') ++i;
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
    // Usa nombres SIN extensión; existing_file probará con y sin ".txt"
    DB.genome_wuhan = read_fasta_dna_clean(existing_file(BASE, "SARS-COV-2-MN908947.3"));
    DB.genome_texas = read_fasta_dna_clean(existing_file(BASE, "SARS-COV-2-MT106054.1"));
    DB.gene_M       = read_fasta_dna_clean(existing_file(BASE, "gen-M"));
    DB.gene_ORF1ab  = read_fasta_dna_clean(existing_file(BASE, "gen-ORF1AB"));
    DB.gene_S       = read_fasta_dna_clean(existing_file(BASE, "gen-S"));
    DB.proteins     = read_protein_fasta(existing_file(BASE, "seq-proteins"));
}


//=============Busqueda Naive 
inline bool match_at (const string& text, size_t i, const string& pattern){
    if (i + pattern.size() > text.size()){
        return false;
    }
    for(size_t k{0}; k < pattern.size(); ++k){
        if(text[i + k] != pattern[k]){
            return false;
        }
    }
    return true;
}

vector<size_t> naive_search_all(const string& text, const string& pattern){
    vector<size_t> pos;
    if(pattern.empty() || text.empty() || pattern.size() > text.size()){
        return pos;
    }
    for(size_t i{0}; i + pattern.size() <= text.size(); i++){
        if(match_at(text, i, pattern)){
            pos.push_back(i);
        } 
    }
    return pos;
}

//ADN basico
inline char comp(char b){
    switch(b){
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N';
    }
}

string revcomp(const string& s){
    string rc;
    rc.resize(s.size());
    for(size_t i {0}, n = s.size(); i < n; ++i){
        rc[n - 1 - i] = comp(s[i]);
    }
    return rc;
}

//Palindromo rev-comp mas largo: solo centros pares (i, i+1)
pair<size_t,size_t> longest_revcomp_palindrome(const string& s){
    auto eq_comp = [&](char a, char b){
        return comp(a) == b;
    };
    size_t bestL = 0, bestR = 0;
    if(s.size() < 2){
        return {0,0};
    }
    for(size_t i{0}; i + 1 < s.size(); ++i){
        long long L = (long long)i, R = (long long) i + 1;
        while(L >= 0 && (size_t)R < s.size() && eq_comp(s[(size_t)L], s[(size_t)R])){
            if((size_t)(R - L) > (bestR - bestL)) {
                bestL = (size_t)L;
                bestR = (size_t)R;
            }
            --L;
            ++R; // expandimos hacia afuera
        }
    }
    return {bestL, bestR}; //si no hay, devuelve [0,0]
}

//=====Codigo genetico y traduccion

// 1) Tabla del código genético
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

// 2) Traduce un marco (definida ANTES de usarla abajo)
string translate_frame(const string& dna, size_t frame, bool until_stop=false){
    string aa;
    if (frame>2 || dna.size() < frame+3) return aa;
    aa.reserve(dna.size()/3);
    for(size_t i=frame; i+2<dna.size(); i+=3){
        string codon; 
        codon.push_back(dna[i]); 
        codon.push_back(dna[i+1]); 
        codon.push_back(dna[i+2]);
        auto it = GENETIC_CODE.find(codon);
        char c = (it==GENETIC_CODE.end()) ? 'X' : it->second;
        if (until_stop && c=='*') break;
        aa.push_back(c);
    }
    return aa;
}

// 3) Usa translate_frame para dar los 3 marcos forward y reverse
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

// buscar AA en el genoma (6 marcos) y mapear a nt
struct HitAA{
    bool reverse; 
    int frame;
    size_t aa_pos; 
    size_t nt_start; 
    size_t nt_end; // inclusive
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
    // reverse
    string rc = revcomp(genome);
    auto rev = translate_3_frames_reverse(genome); // usa rc internamente
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

//===============Menus
void menu_find_gene(){
    cout << "\n[1] Buscar gen en el genoma (búsqueda naive propia)\n";
    cout << "   Gen: 1) M  2) ORF1ab  3) S  -> " << flush;
    int op; 
    if(!(cin >> op)){ cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n'); cout << "   Entrada inválida.\n"; return; }

    const string* gene = nullptr;
    string gname;
    if (op == 1)      { gene = &DB.gene_M;      gname = "M"; }
    else if (op == 2) { gene = &DB.gene_ORF1ab; gname = "ORF1ab"; }
    else if (op == 3) { gene = &DB.gene_S;      gname = "S"; }
    else { cout << "   Opción inválida.\n"; return; }

    cout << "   Genoma: 1) Wuhan 2019  2) Texas 2020  -> " << flush;
    int gg; 
    if(!(cin >> gg)){ cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n'); cout << "   Entrada inválida.\n"; return; }

    const string* genome = (gg==2) ? &DB.genome_texas : &DB.genome_wuhan;
    string glabel = (gg==2) ? "Texas 2020" : "Wuhan 2019";

    auto occ = naive_search_all(*genome, *gene);
    if (occ.empty()){
        cout << "   No se encontró el gen " << gname << " en " << glabel << ".\n";
        return;
    }
    cout << "   Ocurrencias del gen " << gname << " en " << glabel << " (índices 0-based):\n";
    for(size_t i=0;i<occ.size();++i){
        cout << "   - inicio: " << occ[i]
             << "  fin: " << (occ[i] + gene->size() - 1) << "\n";
    }
}

void menu_longest_palindrome(){
    cout << "\n[2] Palíndromo (complemento inverso) más largo dentro de un gen\n";
    cout << "   Gen: 1) M  2) ORF1ab  3) S  -> " << flush;
    int op; if(!(cin >> op)){ cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n'); cout<<"   Entrada inválida.\n"; return; }

    const string* gene = nullptr; string gname;
    if (op==1){ gene=&DB.gene_M; gname="M"; }
    else if (op==2){ gene=&DB.gene_ORF1ab; gname="ORF1ab"; }
    else if (op==3){ gene=&DB.gene_S; gname="S"; }
    else { cout << "   Opción inválida.\n"; return; }

    auto LR = longest_revcomp_palindrome(*gene);
    size_t L = LR.first, R = LR.second;
    if (R <= L){ cout << "   No se detectó palíndromo rev-comp no trivial.\n"; return; }

    cout << "   Gen " << gname << " → mejor palíndromo (rev-comp)\n";
    cout << "   L=" << L << "  R=" << R << "  (longitud=" << (R - L + 1) << ")\n";

    string pal = gene->substr(L, R - L + 1);
    for(size_t i=0; i<pal.size(); i+=60)
        cout << "   " << pal.substr(i, min<size_t>(60, pal.size()-i)) << "\n";
}

// --- Traducción de un gen a proteína (6 marcos)
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

    cout << "   Tip: el marco “real” tiende a comenzar en ATG y tener pocos '*'.\n";
}

// [4] Buscar proteína AA en el genoma (6 marcos), con ingreso por secuencia o por nombre
void menu_find_protein_in_genome(){
    cout << "\n[4] Buscar una proteína (AA) en el GENOMA (6 marcos)\n";

    // Si hay base de proteínas, damos opción de pegar AA o elegir por nombre
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

    string aa;           // secuencia AA elegida
    string chosen_name;  // nombre (si se eligió por lista)

    if (mode == 2 && !DB.proteins.empty()){
        // --- flujo por NOMBRE ---
        cout << "   Mostrar primeros 20? (y/n): ";
        char yn; cin >> yn;
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
                aa = DB.proteins[idx].aa;
                chosen_name = DB.proteins[idx].name;
                ok = true;
            }
        }
        if (!ok){
            // búsqueda por substring (case-insensitive)
            string q = upper_copy(token);
            vector<size_t> matches;
            for(size_t i=0;i<DB.proteins.size();++i){
                string up = upper_copy(DB.proteins[i].name);
                if (up.find(q) != string::npos) matches.push_back(i);
            }
            if (matches.empty()){ cout << "   Nombre no encontrado.\n"; return; }

            if (matches.size() > 20){
                cout << "   " << matches.size() << " coincidencias, muestro primeras 20:\n";
                matches.resize(20);
            }
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
            aa = DB.proteins[idx].aa;
            chosen_name = DB.proteins[idx].name;
        }
    } else {
        // --- flujo pegando SEC. AA ---
        cout << "   Escribe la SECUENCIA de AA (A-Z, '*' si quieres):\n";
        string raw; cin >> raw;
        aa = only_AA(raw);
        if (aa.empty()){ cout << "   Secuencia vacía o inválida.\n"; return; }
        if (aa.size() < 6)
            cout << "   Aviso: secuencia muy corta, puede haber muchos falsos negativos.\n";
    }

    cout << "   Genoma: 1) Wuhan 2019  2) Texas 2020  -> " << flush;
    int gg; 
    if(!(cin >> gg)){ 
        cin.clear(); cin.ignore(numeric_limits<streamsize>::max(), '\n'); 
        cout<<"   Entrada inválida.\n"; return; 
    }
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
    for(const auto& h: hits){
        cout << "   - " << (h.reverse ? "REV" : "FWD")
             << " frame " << h.frame
             << " | nt[" << h.nt_start << "," << h.nt_end << "]"
             << " | aa_pos=" << h.aa_pos << "\n";
    }
}

// [5] Buscar proteína por NOMBRE (usa seq-proteins.txt)
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
    for (const auto& h : hits) {
        cout << "   - " << (h.reverse ? "REV" : "FWD") 
             << " frame " << h.frame 
             << " | nt[" << h.nt_start << "," << h.nt_end << "]"
             << " | aa_pos=" << h.aa_pos << "\n";
    }
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
            cout << "1) Buscar un gen en el genoma\n";
            cout << "2) Palíndromo rev-comp más largo en un gen\n";
            cout << "3) Traducir un gen a proteína (6 marcos)\n";
            cout << "4) Buscar proteína (AA) en el genoma (6 marcos)\n";
            cout << "5) Buscar proteína por NOMBRE (seq-proteins.txt)\n";
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
            else cout << "   Opción inválida.\n";
        }

    }catch(const exception& e){
        cerr << "[ERROR] " << e.what() << "\n";
        return 1;
    }
}
