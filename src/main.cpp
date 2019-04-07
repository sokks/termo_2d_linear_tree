#include "proc/proc.h"

using std::string;


bool WRITE_LAYERS = true;
int write_freq = 1000;
string baseFolderTemp = "data/temp/";

int refine_freq = 1000;
string baseFolderRefine = "data/refine/";

string gen_filename(string baseFolder, int n) {
    string num = std::to_string(n);
    int max = 6;
    int additional = max - num.length();
    string res = string(additional, '0') + num;
    string fname = baseFolder + res + ".out";
    return fname;
}


int main(int argc, char **argv) {
    if (argc < 3) {
        std::cout << "usage: prog <BaseSize> <TimeSteps>\n";
        return 0;
    }
    int base_sz, ts_n;
    base_sz = std::atoi(argv[1]);
    ts_n    = std::atoi(argv[2]);

    Init(base_sz, ts_n);    

    Proc p;
    p.MPIInit(argc, argv);
    
    p.InitMesh();
    p.FillStart(&Area::T0);
    p.MarkToRefine();
    p.Refine();
    p.LoadBalance();

    p.BuildGhosts();
    for (int k = 0; k < ts_n; k++) {
        if (WRITE_LAYERS && (k%write_freq ==0)) {
            p.WriteT(gen_filename(baseFolderTemp, k));
        }
        p.MakeStep();
    }

    p.WriteT(gen_filename(baseFolderTemp, ts_n));
    p.WriteStat("data/stat.out");

    p.MPIFinalize();
}