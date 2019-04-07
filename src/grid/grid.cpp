#include "grid.h"

using std::vector;
using std::map;
using std::string;
using std::sort;


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
vector<GlobalNumber_t> CellIndex::get_all_possible_neighbours_ids(MpiTimer& timer) {
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

    timer.Start();
    sort(all_neighs_ids.begin(), all_neighs_ids.end());
    timer.Stop();

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
    // можно сразу задать min_dx и использовать его
    double min_dx = base_dx/pow(2,max_lvl-lvl);
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


LinearTree::LinearTree(string filename) {
    // todo
}

LinearTree::LinearTree(int base_level, int max_level, double (*Temp_func)(double, double)) {
    
    base_lvl = base_level;
    max_lvl  = max_level;

    base_sz  = pow(2, base_lvl);
    base_dx  = (Area::x_end - Area::x_start) / base_sz;

    int max_sz = pow(2, max_lvl);
    min_dx  = (Area::x_end - Area::x_start) / max_sz;

    int global_i_start = 0;
    int global_i_stop  = 1 << (max_lvl*2);
    int global_i_step  = 1 << (2*max_lvl - 2*base_lvl);

    // int my_i_start = offsetG * global_i_step;
    // int my_i_stop = my_i_start + procG * (global_i_step);
    int my_i_start = global_i_start;
    int my_i_stop  = global_i_stop;

    /* std::cout << mpiInfo.comm_rank <<  
        " global_i_start=" << global_i_start <<
        " global_i_stop=" << global_i_stop <<
        " global_i_step=" << global_i_step <<
        " my_i_start=" << my_i_start <<
        " my_i_stop=" << my_i_stop << std::endl; */

    
    for (int i = my_i_start; i < my_i_stop; i += global_i_step) {
        tree.cells.push_back(Cell(base_lvl, i));
    }

    // fill metaInfo (пока вообще обмен не нужен, но если загружать сетку извне, то будет нужно)
    // meta = new MetaInfo[mpiInfo.comm_size];
    // meta[mpiInfo.comm_rank] = MetaInfo{ my_i_start, my_i_stop - global_i_step, procG };
    // MPI_Allgather(&meta[mpiInfo.comm_rank], sizeof(MetaInfo), MPI_CHAR, (void *)meta, sizeof(MetaInfo), MPI_CHAR, mpiInfo.comm);

    // std::cout << mpiInfo.comm_rank <<  " exchanged meta  ";
    // for (int i = 0; i < mpiInfo.comm_size; i++) {
    //     std::cout << meta[i].toString() << " | ";
    // }
}



// int LinearTree::FindCell(int lvl, int i, int j) {
//     GlobalNumber_t target = merge_ints(i, j);
//     return FindCell(target);
// }

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