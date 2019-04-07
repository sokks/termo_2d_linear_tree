#include "grid/grid.h"

int main(int argc, char **argv) {
    int base_level = 0;
    int max_level  = 10;
    string filename = "data/refine/start_grid.dat";
    
    if (argc >= 4) {
        base_level = atoi(argv[1]);
        max_level  = atoi(argv[2]);
        filename   = argv[3]; 
    }
    GridInit(base_level, max_level);
    LinearTree grid = LinearTree(&Area::T0);
    while (grid.MarkToRefine()) {
        grid.DoRefine();
    }
    // grid.Balance21();

    grid.Write(filename);
}