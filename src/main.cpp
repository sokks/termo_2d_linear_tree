#include "proc/proc.h"

using std::string;
using std::cout;
using std::endl;


bool WRITE_LAYERS = true;
int write_freq = 100;
string baseFolderTemp = "data/temp/";

string from_num(int n) {
    string a = "";
    while (n > 0) {
        a = string(1, char(n%10) + '0') + a;
        n = n / 10;
    }

    return a;
}

char * gen_filename(string baseFolder, int n) {
    string num = from_num(n);
    int max = 6;
    int additional = max - num.length();
    string res = string(additional, '0') + num;
    string fname = baseFolder + res + ".out";

    char *fmane_c = new char[fname.size()];
    strcpy(fname_c, fname.c_str())
    return fname_c;
}


int main(int argc, char **argv) {
    if (argc < 6) {
        std::cout << "usage: prog <base_lvl> <max_lvl> <offsets_file> <grid_file> <time_steps>\n";
        return 0;
    }
    int    base_level = atoi(argv[1]);
    int    max_level  = atoi(argv[2]);
    char * offsets_file = argv[3];
    char * grid_file = argv[4];
    int    ts_n       = std::atoi(argv[5]);

    GridInit(base_level, max_level);
    SolverInit(ts_n);   

    Proc p;
    p.MPIInit(argc, argv);

    cout << "base_level=" << base_level << 
            " max_level=" << max_level <<
            " offsets_file=" << offsets_file << 
            " grid_file=" << grid_file << 
            " ts_n=" << ts_n << endl;
    
    p.InitMesh(offsets_file, grid_file);

    // usleep(10000000);

    p.BuildGhosts();
    p.BuildNeighs();

    cout << "PROC_MEMORY=" << p.GetProcAllocMem() << endl;

    // usleep(10000000);

    for (int k = 0; k < ts_n; k++) {
        cout << k  << endl;
        if (WRITE_LAYERS && (k%write_freq ==0)) {
            p.WriteT(gen_filename(baseFolderTemp, k));
        }
        p.MakeStep();
    }

    // p.WriteT(gen_filename(baseFolderTemp, ts_n));
    // p.WriteStat("data/stat.out");

    p.MPIFinalize();
}