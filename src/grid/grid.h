#include <vector>
#include <mpi.h>
#include <math.h>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>

#include "../area/area.h"

using std::vector;
using std::map;
using std::string;

/* * * * * * * * * * РАБОТА С ИНДЕКСАМИ * * * * * * * * * */

typedef long long int GlobalNumber_t;

GlobalNumber_t merge_ints(int a, int b);
void split_ints(GlobalNumber_t c, int *a, int *b);


enum Neigh   { DOWN, UP, LEFT, RIGHT };
enum CornerNeigh { LU, RU, LD, RD };
enum Child       { cLU, cRU, cLD, cRD };

struct CellIndex {
    int lvl, i, j;

    CellIndex() {}
    CellIndex(int _lvl, int _i, int _j): lvl(_lvl), i(_i), j(_j) {}
    CellIndex(const CellIndex& c): lvl(c.lvl), i(c.i), j(c.j) {}
    CellIndex(int _lvl, GlobalNumber_t globalNumber);
    CellIndex& operator=(const CellIndex& c) { lvl = c.lvl; i = c.i; j = c.j; return *this; }
    // CellIndex operator=(CellIndex c) { lvl = c.lvl; i = c.i; j = c.j; return *this; }
    // ~CellIndex() {}

    GlobalNumber_t get_global_number();

    bool is_left_border();
    bool is_right_border();
    bool is_upper_border();
    bool is_down_border();
    bool is_left_upper_corner();
    bool is_right_upper_corner();
    bool is_left_down_corner();
    bool is_right_down_corner();

    CellIndex get_child(Child c);
    CellIndex get_parent();
    Child get_child_pos();

    CellIndex get_face_neighbor(Neigh n);
    CellIndex get_corner_neighbor(CornerNeigh n);
    vector<CellIndex>  get_larger_possible_face_neighbour(Neigh n);
    vector<CellIndex>& get_larger_possible_face_neighbour_optimized(vector<CellIndex>&, Neigh n);
    vector<CellIndex>  get_halfsize_possible_face_neighbours(Neigh n);
    vector<CellIndex>& get_halfsize_possible_face_neighbours_optimized(vector<CellIndex>&, Neigh n);
    
    vector<CellIndex>  get_larger_possible_corner_neighbour(CornerNeigh n);
    vector<CellIndex>  get_halfsize_possible_corner_neighbours(CornerNeigh n);

    vector<CellIndex>  get_all_halfsize_possible_neighs();
    vector<CellIndex>& get_all_halfsize_possible_neighs_optimized(vector<CellIndex>&);
    vector<CellIndex> get_all_samesize_possible_neighs();
    vector<CellIndex>& get_all_samesize_possible_neighs_optimized(vector<CellIndex>&);
    vector<CellIndex> get_all_larger_possible_neighs();
    vector<CellIndex>& get_all_larger_possible_neighs_optimized(vector<CellIndex>&);

    vector<GlobalNumber_t> get_all_possible_neighbours_ids(MpiTimer& timer);

    bool is_border();
};



/* * * * * * * * РАБОТА С ЯЧЕЙКАМИ И ДЕРЕВОМ * * * * * * * */

double dist(double x1, double y1, double x2, double y2);

struct Cell: public CellIndex {
    char refine_mark = 0;
    double temp[2]; // for cur and next

    Cell() {}
    Cell(int _lvl, GlobalNumber_t globalNumber): CellIndex(_lvl, globalNumber) {}
    Cell(const Cell& c): CellIndex(c.lvl, c.i, c.j) {temp[0] = c.temp[0]; temp[1] = c.temp[1]; }
    Cell& operator=(const Cell& c) { lvl = c.lvl; i = c.i; j = c.j; temp[0] = c.temp[0]; temp[1] = c.temp[1]; return *this; }

    double *get_temp() { return &temp[0]; }
    void    get_spacial_coords(double *x, double *y);
    void    get_border_cond(char *cond_type, double (**cond_func)(double, double, double));
    double  get_S() {
        double lvl_dx = min_dx * pow(2, max_lvl - lvl);
        return lvl_dx * lvl_dx;
    }
};


struct LinearTree {

    int    base_sz;
    double base_dx;
    double min_dx;
    int    max_lvl;
    int    base_lvl;

    vector<Cell> cells;

    LinearTree() {}
    LinearTree(string filename);
    LinearTree(int base_size, int max_level, double (*Temp_func)(double, double));

    int FindCell(GlobalNumber_t target, Cell *cell);

    void MarkToRefine();
    void DoRefine();

    void Write(string filename);
};