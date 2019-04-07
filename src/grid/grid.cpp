#include "grid.h"

using std::vector;
using std::map;
using std::string;
using std::sort;
using std::cout;
using std::endl;

int    base_sz;
double base_dx;
double min_dx;
int    max_lvl;
int    base_lvl;

void GridInit(int base_level, int max_level) {
    base_lvl = base_level;
    max_lvl  = max_level;

    base_sz  = pow(2, base_lvl);
    base_dx  = (Area::x_end - Area::x_start) / base_sz;

    int max_sz = pow(2, max_lvl);
    min_dx  = (Area::x_end - Area::x_start) / max_sz;

    cout << "GRID PARAMS INITED" <<
            "\nmax_lvl=" << max_lvl <<
            "\nmin_dx=" << min_dx <<
            "\nbase_lvl=" << base_sz <<
            "\nbase_dx=" << base_dx <<
            "\nbase_sz=" << base_sz << endl;
}


/* * * * * * * * * * РАБОТА С ИНДЕКСАМИ * * * * * * * * * */

GlobalNumber_t merge_ints(int a, int b) {
    GlobalNumber_t c = 0;

    int c_pos = max_lvl * 2;
    int a_pos = max_lvl;
    int b_pos = max_lvl;

    for (int pos = max_lvl - 1; pos >= 0; pos--) {
        int a1 = a & (1 << pos);
        GlobalNumber_t a2 = a1 << (pos + 1);
        c = c | a2;
        
        int b1 = b & (1 << pos);
        GlobalNumber_t b2 = b1 << pos;
    }

    return c;
}

void split_ints(GlobalNumber_t c, int *a, int *b) {
    int aa = 0;
    int bb = 0;

    for (int pos = max_lvl - 1; pos >= 0; pos--) {
        GlobalNumber_t c1 = c & (1 << (pos*2+1));
        c1 = c1 >> (pos + 1);
        aa = aa | c1;

        GlobalNumber_t c2 = c & (1 << (pos*2));
        c2 = c2 >> (pos);
        bb = bb | c2;
    }

    *a = aa;
    *b = bb;
}

CellIndex::CellIndex(int _lvl, GlobalNumber_t globalNumber): lvl(_lvl) {
    std::cout << "new CellIndex(" << lvl << ", " << globalNumber << ")\n";
    split_ints(globalNumber, &i, &j);
}

GlobalNumber_t CellIndex::get_global_number() {
    return merge_ints(i, j);
}

CellIndex CellIndex::get_child(Child c) {
    CellIndex child;
    child.lvl = lvl+1;
    int h = 1 << (max_lvl - child.lvl);

    if (c == Child::cLD) {        // (0,0)
        child.i = i;
        child.j = j;
    } else if (c == Child::cRD) { // (0,1)
        child.i = i;
        child.j = j + h;
    } else if (c == Child::cLU) { // (1,0)
        child.i = i + h;
        child.j = j;
    } else {                     // (1,1)
        child.i = i + h;
        child.j = j + h;
    }

    return child;
}

CellIndex CellIndex::get_parent() {
    CellIndex parent;
    parent.lvl = lvl - 1;
    parent.i = i & ~(1 << (max_lvl - lvl));
    parent.j = j & ~(1 << (max_lvl - lvl));
    return parent;
}

Child CellIndex::get_child_pos() {
    int hi = i & (1 << (max_lvl - lvl));
    int hj = j & (1 << (max_lvl - lvl));
    if ((hi == 0) && (hj == 0)) {
        return Child::cLD;
    }
    if ((hi == 0) && (hj != 0)) {
        return Child::cRD;
    }
    if ((hi != 0) && (hj == 0)) {
        return Child::cLU;
    }
    if ((hi != 0) && (hj != 0)) {
        return Child::cRU;
    }
    return Child::cLD;
}

CellIndex CellIndex::get_face_neighbor(Neigh n) {
    int h = 1 << (max_lvl - lvl - 1); // 2^(b-l)
    CellIndex c;
    c.lvl = lvl;
    c.i   = i + ((n == Neigh::DOWN) ? -h : (n == Neigh::UP) ? h : 0);
    c.j   = j + ((n == Neigh::LEFT) ? -h : (n == Neigh::RIGHT) ? h : 0);
    return c;
}

CellIndex CellIndex::get_corner_neighbor(CornerNeigh n) {
    int h = 1 << (max_lvl - lvl - 1); // 2^(b-l)
    CellIndex c;
    c.lvl = lvl;
    c.i   = i + ( ((n == CornerNeigh::LU) || (n == CornerNeigh::LD)) ? -h : h);
    c.j   = j + ( ((n == CornerNeigh::LD) || (n == CornerNeigh::RD)) ? -h : h);
    return c;
}

vector<CellIndex> CellIndex::get_larger_possible_face_neighbour(Neigh n) {
    vector<CellIndex> ret;

    Child my_pos = get_child_pos();
    if (n == Neigh::RIGHT) {
        if ((my_pos == Child::cLD) || (my_pos == Child::cLU)) {
            return ret;
        }
        CellIndex c = get_parent();
        ret.push_back(c.get_face_neighbor(n));
        return ret;
    }

    if (n == Neigh::LEFT) {
        if ((my_pos == Child::cRD) || (my_pos == Child::cRU)) {
            return ret;
        }
        CellIndex c = get_parent();
        ret.push_back(c.get_face_neighbor(n));
        return ret;
    }

    if (n == Neigh::UP) {
        if ((my_pos == Child::cLD) || (my_pos == Child::cRD)) {
            return ret;
        }
        CellIndex c = get_parent();
        ret.push_back(c.get_face_neighbor(n));
        return ret;
    }
    
    if (n == Neigh::DOWN) {
        if ((my_pos == Child::cLU) || (my_pos == Child::cRU)) {
            return ret;
        }
        CellIndex c = get_parent();
        ret.push_back(c.get_face_neighbor(n));
        return ret;
    }

    return ret;
}

vector<CellIndex>& CellIndex::get_larger_possible_face_neighbour_optimized(vector<CellIndex>& buf, Neigh n) {
    Child my_pos = get_child_pos();
    if (n == Neigh::RIGHT) {
        if ((my_pos == Child::cLD) || (my_pos == Child::cLU)) {
            return buf;
        }
        CellIndex c = get_parent();
        buf.push_back(c.get_face_neighbor(n));
        return buf;
    }

    if (n == Neigh::LEFT) {
        if ((my_pos == Child::cRD) || (my_pos == Child::cRU)) {
            return buf;
        }
        CellIndex c = get_parent();
        buf.push_back(c.get_face_neighbor(n));
        return buf;
    }

    if (n == Neigh::UP) {
        if ((my_pos == Child::cLD) || (my_pos == Child::cRD)) {
            return buf;
        }
        CellIndex c = get_parent();
        buf.push_back(c.get_face_neighbor(n));
        return buf;
    }
    
    if (n == Neigh::DOWN) {
        if ((my_pos == Child::cLU) || (my_pos == Child::cRU)) {
            return buf;
        }
        CellIndex c = get_parent();
        buf.push_back(c.get_face_neighbor(n));
        return buf;
    }

    return buf;
}

vector<CellIndex> CellIndex::get_larger_possible_corner_neighbour(CornerNeigh n) {
    vector<CellIndex> ret;

    Child my_pos = get_child_pos();
    if (n == CornerNeigh::LD) {
        if (my_pos != Child::cRU) {
            return ret;
        }
    }
    if (n == CornerNeigh::LU) {
        if (my_pos != Child::cRD) {
            return ret;
        }
    }
    if (n == CornerNeigh::LD) {
        if (my_pos != Child::cRU) {
            return ret;
        }
    }
    if (n == CornerNeigh::LD) {
        if (my_pos != Child::cRU) {
            return ret;
        }
    }

    CellIndex c = get_parent();
    ret.push_back(c.get_corner_neighbor(n));
    return ret;
}

vector<CellIndex> CellIndex::get_halfsize_possible_face_neighbours(Neigh n) {
    CellIndex fullSizeNeigh = get_face_neighbor(n);
    vector<CellIndex> halfSizeNeighs;
    halfSizeNeighs.reserve( 2 );
    if (n == Neigh::DOWN) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cLU));
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cRU));
        return halfSizeNeighs;
    }
    if (n == Neigh::UP) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cLD));
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cRD));
        return halfSizeNeighs;
    }
    if (n == Neigh::LEFT) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cRU));
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cRD));
        return halfSizeNeighs;
    }
    if (n == Neigh::RIGHT) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cLU));
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cLD));
        return halfSizeNeighs;
    }
    return halfSizeNeighs;
}

vector<CellIndex>& CellIndex::get_halfsize_possible_face_neighbours_optimized(vector<CellIndex>& buf, Neigh n) {
    CellIndex fullSizeNeigh = get_face_neighbor(n);
    if (n == Neigh::DOWN) {
        buf.push_back(fullSizeNeigh.get_child(Child::cLU));
        buf.push_back(fullSizeNeigh.get_child(Child::cRU));
        return buf;
    }
    if (n == Neigh::UP) {
        buf.push_back(fullSizeNeigh.get_child(Child::cLD));
        buf.push_back(fullSizeNeigh.get_child(Child::cRD));
        return buf;
    }
    if (n == Neigh::LEFT) {
        buf.push_back(fullSizeNeigh.get_child(Child::cRU));
        buf.push_back(fullSizeNeigh.get_child(Child::cRD));
        return buf;
    }
    if (n == Neigh::RIGHT) {
        buf.push_back(fullSizeNeigh.get_child(Child::cLU));
        buf.push_back(fullSizeNeigh.get_child(Child::cLD));
        return buf;
    }

    return buf;
}

vector<CellIndex> CellIndex::get_halfsize_possible_corner_neighbours(CornerNeigh n) {
    CellIndex fullSizeNeigh = get_corner_neighbor(n);
    vector<CellIndex> halfSizeNeighs;
    if (n == CornerNeigh::LU) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cRD));
        return halfSizeNeighs;
    }
    if (n == CornerNeigh::RU) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cLD));
        return halfSizeNeighs;
    }
    if (n == CornerNeigh::LD) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cRU));
        return halfSizeNeighs;
    }
    if (n == CornerNeigh::RD) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cLU));
        return halfSizeNeighs;
    }
    return halfSizeNeighs;
}

vector<CellIndex> CellIndex::get_all_halfsize_possible_neighs() {
    vector<CellIndex> res;
    res.reserve( 8 );

    res = get_halfsize_possible_face_neighbours_optimized(res, Neigh::LEFT);
    res = get_halfsize_possible_face_neighbours_optimized(res, Neigh::RIGHT);
    res = get_halfsize_possible_face_neighbours_optimized(res, Neigh::UP);
    res = get_halfsize_possible_face_neighbours_optimized(res, Neigh::DOWN);

    return res;
}

vector<CellIndex>& CellIndex::get_all_halfsize_possible_neighs_optimized(vector<CellIndex>& buf) {

    buf = get_halfsize_possible_face_neighbours_optimized(buf, Neigh::LEFT);
    buf = get_halfsize_possible_face_neighbours_optimized(buf, Neigh::RIGHT);
    buf = get_halfsize_possible_face_neighbours_optimized(buf, Neigh::UP);
    buf = get_halfsize_possible_face_neighbours_optimized(buf, Neigh::DOWN);

    return buf;
}

vector<CellIndex> CellIndex::get_all_samesize_possible_neighs() {
    vector<CellIndex> res;
    res.reserve( 4 );

    res.push_back(get_face_neighbor(Neigh::LEFT));
    res.push_back(get_face_neighbor(Neigh::RIGHT));
    res.push_back(get_face_neighbor(Neigh::DOWN));
    res.push_back(get_face_neighbor(Neigh::UP));

    return res;
}

vector<CellIndex>& CellIndex::get_all_samesize_possible_neighs_optimized(vector<CellIndex>& buf) {

    buf.push_back(get_face_neighbor(Neigh::LEFT));
    buf.push_back(get_face_neighbor(Neigh::RIGHT));
    buf.push_back(get_face_neighbor(Neigh::DOWN));
    buf.push_back(get_face_neighbor(Neigh::UP));

    return buf;
}

vector<CellIndex> CellIndex::get_all_larger_possible_neighs() {
    vector<CellIndex> res;
    res.reserve(4);

    res = get_larger_possible_face_neighbour_optimized(res, Neigh::LEFT);
    res = get_larger_possible_face_neighbour_optimized(res, Neigh::RIGHT);
    res = get_larger_possible_face_neighbour_optimized(res, Neigh::UP);
    res = get_larger_possible_face_neighbour_optimized(res, Neigh::DOWN);

    return res;
}

vector<CellIndex>& CellIndex::get_all_larger_possible_neighs_optimized(vector<CellIndex>& buf) {

    buf = get_larger_possible_face_neighbour_optimized(buf, Neigh::LEFT);
    buf = get_larger_possible_face_neighbour_optimized(buf, Neigh::RIGHT);
    buf = get_larger_possible_face_neighbour_optimized(buf, Neigh::UP);
    buf = get_larger_possible_face_neighbour_optimized(buf, Neigh::DOWN);

    return buf;
}

// TODO optimize copying
vector<GlobalNumber_t> CellIndex::get_all_possible_neighbours_ids() {
    vector<CellIndex> all_neighs;
    all_neighs.reserve(20);

    all_neighs = get_all_halfsize_possible_neighs_optimized(all_neighs);
    all_neighs = get_all_samesize_possible_neighs_optimized(all_neighs);
    all_neighs = get_all_larger_possible_neighs_optimized(all_neighs);
    
    vector<GlobalNumber_t> all_neighs_ids;
    all_neighs_ids.reserve(all_neighs.size());
    for (CellIndex c: all_neighs) {
        all_neighs_ids.push_back(c.get_global_number());
    }

    sort(all_neighs_ids.begin(), all_neighs_ids.end());

    vector<GlobalNumber_t> all_neighs_ids_without_duplicates;
    all_neighs_ids_without_duplicates.reserve(all_neighs_ids.size());
    GlobalNumber_t prev = all_neighs_ids[0];
    all_neighs_ids_without_duplicates.push_back(prev);
    for (auto it = ++all_neighs_ids.begin(); it < all_neighs_ids.end(); it++) {
        if (*it != prev) {
            all_neighs_ids_without_duplicates.push_back(*it);
        }
        prev = *it;
    }
    return all_neighs_ids_without_duplicates;
}

bool CellIndex::is_left_border() {
    CellIndex c = get_face_neighbor(Neigh::LEFT);
    return (c.j < 0);
}

bool CellIndex::is_right_border() {
    CellIndex c = get_face_neighbor(Neigh::RIGHT);
    return ((c.j & (-1 << max_lvl)) != 0);
}

bool CellIndex::is_upper_border() {
    CellIndex c = get_face_neighbor(Neigh::UP);
    return ((c.i & (-1 << max_lvl)) != 0);
}

bool CellIndex::is_down_border() {
    CellIndex c = get_face_neighbor(Neigh::DOWN);
    return (c.i < 0);
}

bool CellIndex::is_left_upper_corner() {
    return is_left_border() && is_upper_border();
}

bool CellIndex::is_right_upper_corner() {
    return is_right_border() && is_upper_border();
}

bool CellIndex::is_left_down_corner() {
    return is_left_border() && is_down_border();
}

bool CellIndex::is_right_down_corner() {
    return is_right_border() && is_down_border();
}

bool CellIndex::is_border() {
    return is_left_border() || is_right_border() || is_upper_border() || is_down_border();
            // is_left_upper_corner() || is_right_upper_corner() || is_left_down_corner() || is_right_down_corner();
}


/* * * * * * * * РАБОТА С ЯЧЕЙКАМИ И ДЕРЕВОМ * * * * * * * */

double dist(double x1, double y1, double x2, double y2) {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

void Cell::get_spacial_coords(double *x, double *y) {
    double lvl_dx = min_dx * pow(2, max_lvl - lvl);
    *x = min_dx * i + lvl_dx/2; 
    *y = min_dx * j + lvl_dx/2; 
}

void Cell::get_border_cond(char *cond_type, double (**cond_func)(double, double, double)) {
    if (is_left_border()) { 
        Area::get_border_cond(Area::Border::LEFT, cond_type, cond_func);
        return;
    }

    if (is_right_border()) { 
        Area::get_border_cond(Area::Border::RIGHT, cond_type, cond_func);
        return;
    }

    if (is_down_border()) { 
        Area::get_border_cond(Area::Border::UP, cond_type, cond_func);
        return;
    }

    if (is_upper_border()) { 
        Area::get_border_cond(Area::Border::DOWN, cond_type, cond_func);
        return;
    }

    *cond_type = -1;
}

vector<Cell> Cell::split() {
    vector<Cell> children;
    
    CellIndex ci00 = ((CellIndex)(*this)).get_child(Child::cLD);
    Cell c00 = Cell(ci00, 0.0);
    double x, y;
    c00.get_spacial_coords(&x, &y);
    c00.temp[0] = Area::T0(x, y);
    children.push_back(c00);

    CellIndex ci01 = ((CellIndex)(*this)).get_child(Child::cRD);
    Cell c01 = Cell(ci01, 0.0);
    c01.get_spacial_coords(&x, &y);
    c01.temp[0] = Area::T0(x, y);
    children.push_back(c01);

    CellIndex ci10 = ((CellIndex)(*this)).get_child(Child::cLU);
    Cell c10 = Cell(ci10, 0.0);
    c10.get_spacial_coords(&x, &y);
    c10.temp[0] = Area::T0(x, y);
    children.push_back(c10);

    CellIndex ci11 = ((CellIndex)(*this)).get_child(Child::cRU);
    Cell c11 = Cell(ci11, 0.0);
    c11.get_spacial_coords(&x, &y);
    c11.temp[0] = Area::T0(x, y);
    children.push_back(c11);

    return children;
}

LinearTree::LinearTree(string filename) {
    // todo
}

LinearTree::LinearTree(double (*Temp_func)(double, double)) {
    
    max_present_lvl = base_lvl;

    int global_i_start = 0;
    int global_i_stop  = 1 << (max_lvl*2);
    int global_i_step  = 1 << (2*max_lvl - 2*base_lvl);

    for (int i = global_i_start; i < global_i_stop; i += global_i_step) {
        Cell c = Cell(base_lvl, i);
        double x, y;
        c.get_spacial_coords(&x, &y);
        c.temp[0] = Temp_func(x, y);
        cells.push_back(c);
    }
}

int LinearTree::FindCell(GlobalNumber_t target, Cell *cell) {
    // binary search
    int left = 0;
	int right = cells.size();
	while (left < right) {
		int midi = (right + left) / 2;
		Cell mid = cells[midi];
        GlobalNumber_t val = mid.get_global_number();
		if (val == target) {
			// return midi;
            *cell = mid;
            return midi;
		}
        if (val < target) {
			left = midi + 1;
		} else {
			right = midi;
		}
	}
	return -1;
}

int LinearTree::MarkToRefine() {
    cout << "MarkToRefine start\n";
    
    if (max_present_lvl == max_lvl) {
        return 0;
    }

    double refine_thresholds[] = {0, 0, 0, 0, 0, 10, 20, 30, 40, 50};

    int n_of_marks = 0;

    for (int i = 0; i < cells.size(); i++) {
        Cell c = cells[i];

        if (c.is_border()) {
            continue;
        }

        // vector<GlobalNumber_t> neighs = c.get_all_possible_neighbours_ids();
        // vector<double> vals = get_square_vals(neighs);
        double grad = 0.0;
        // double grad = get_grad(vals[0],  vals[1],  vals[2],
        //                        vals[3], c.temp[0], vals[4],
        //                        vals[5],  vals[6],  vals[7]);

        // if (grad > refine_thresholds[c.lvl]) {
        //     cells[i].refine_mark = 1;
        //     n_of_marks++;
        // }
        double x, y;
        c.get_spacial_coords(&x, &y);

        if ((max_present_lvl == base_lvl) && (Area::Refine1(x, y))) {
            cells[i].mark_to_refine();
            cells[i].temp[0] = 100;
            cout << int(cells[i].refine_mark)  << " ";
            n_of_marks++;
        }

        if ((max_present_lvl == base_lvl + 1) && (Area::Refine2(x, y))) {
            cells[i].refine_mark = char(1);
            n_of_marks++;
        }

        if ((max_present_lvl == base_lvl + 2) && (Area::Refine3(x, y))) {
            cells[i].refine_mark = 1;
            n_of_marks++;
        }

    }

    cout << "N_OF_MARKS=" << n_of_marks << endl;
    cout << "MarkToRefine end\n";
    return (n_of_marks > 0);
}

// по методу Базарова (из анализа изображений)
double H1[3][3] = {
    { 1,  2,  1},
    { 0,  0,  0},
    {-1, -2, -1}};

double H2[3][3] = {
    {-1, 0, 1},
    {-2, 0, 2},
    {-1, 0, 1}};

double get_grad(double w00, double w01, double w02,
                double w10, double w11, double w12,
                double w20, double w21, double w22, double dx) {
    // H1 * W
    double S1 = (w00 + 2*w01 + w02) - (w20 + 2*w21 + w22);
    // H2 * W
    double S2 =  (w02 + 2*w12 + w22) - (w00 + 2*w10 + w20);
    // модуль градиента
    double G = 1 / (8*dx) * sqrt(S1*S1 + S2*S2);

    return G;
}

void LinearTree::DoRefine() {
    int i = 0;
    cout << "DoRefine start\n";
    while (i < cells.size()) {
        Cell c = cells[i];
        // cout << "i=" << i << " Cell(" << c.lvl << ", " << c.i << "," << c.j << ", " << c.temp[0] << ")";
        if (c.refine_mark) {
            cells.erase(cells.begin()+i);

            vector<Cell> new_cells = c.split();
            // if (i < 200) {
            //     cout << "Cell(" << c.lvl << ", " << c.i << "," << c.j << ", " << c.temp[0] << ") --> ";
            //     cout << "Cell(" << new_cells[0].lvl << ", " << new_cells[0].i << "," << new_cells[0].j << ", " << new_cells[0].temp[0] << ")  ";
            //     cout << "Cell(" << new_cells[1].lvl << ", " << new_cells[1].i << "," << new_cells[1].j << ", " << new_cells[1].temp[0] << ")  ";
            //     cout << "Cell(" << new_cells[2].lvl << ", " << new_cells[2].i << "," << new_cells[2].j << ", " << new_cells[2].temp[0] << ")  ";
            //     cout << "Cell(" << new_cells[3].lvl << ", " << new_cells[3].i << "," << new_cells[3].j << ", " << new_cells[3].temp[0] << ")  ";
            // }

            cells.insert(cells.begin()+(i), new_cells.begin(), new_cells.end());
            i += 4;
            
        } else {
            i++;
        }
    }

    max_present_lvl++;
    cout << "DoRefine end\n";
}

void LinearTree::Balance21() {

}

void LinearTree::Write(string filename) {
    vector<char> buf = GenWriteStruct();
    int len = buf.size();

    std::ofstream fout(filename, std::ios::out | std::ios::binary);
    fout.write(&buf[0], len);
    fout.close();
}

vector<char> LinearTree::GenWriteStruct() {
    vector<char> buf;
    std::cout << " gen structs for write start\n";
    for (Cell c: cells) {
        char *tmp = (char *)(&c.lvl);
        for (int i = 0; i < sizeof(int); i++) {
            buf.push_back(tmp[i]);
        }
        tmp = (char *)(&c.i);
        for (int i = 0; i < sizeof(int); i++) {
            buf.push_back(tmp[i]);
        }
        tmp = (char *)(&c.j);
        for (int i = 0; i < sizeof(int); i++) {
            buf.push_back(tmp[i]);
        }
        if (c.temp[0] > 1) {
            // std::cout << (c.temp[0] > 5) << std::endl;
        }
        
        tmp = (char *)(&c.temp[0]);
        for (int i = 0; i < sizeof(double); i++) {
            buf.push_back(tmp[i]);
        }
    }
    std::cout << " gen structs for write finished\n";
    return buf;
}
