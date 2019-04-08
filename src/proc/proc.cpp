#include "proc.h"

using std::vector;
using std::map;
using std::string;
using std::sort;

using std::cout;
using std::endl;

double base_tau = 0.001;
double      tau = 0.001;
int    time_steps = 100;

void Init(int gs_x, int ts_n) {
    base_sz = gs_x;
    base_lvl = (int) floor(log2(base_sz));

    base_dx = (Area::x_end - Area::x_start) / base_sz;
    base_tau = (base_dx * base_dx) / 4;
    time_steps = ts_n;

    // todo set tau

    std::cout << "a=" << Area::a << std::endl;
    std::cout << "base_sz=" << base_sz << " base_lvl=" << base_lvl << " max_lvl=" << max_lvl << std::endl;
    std::cout << "base_dx=" << base_dx << " base_tau=" << base_tau << " t_end=" << base_tau * time_steps << std::endl;
}


/* * * * * * * * * * * СТАТИСТИКА * * * * * * * * * * */

string MpiInfo::toString() {
    std::ostringstream stringStream;
    stringStream << "{ comm: " << comm << ", ";
    stringStream << "comm_size: " << comm_size << ", ";
    stringStream << "comm_rank: " << comm_rank << " }";

    return stringStream.str();
}

string Stat::toString() {
    std::ostringstream stringStream;
    for (auto t: timers) {
        stringStream << t.first << ": " << t.second.FullDur() << std::endl;
    }
    return stringStream.str();
}

string MetaInfo::toString() {
    std::ostringstream stringStream;
    stringStream << "{ (first: " << procStart << " ";
    stringStream << "last: " << procEnd << "), ";
    stringStream << "len: " << procG << " }";

    return stringStream.str();
}


/* * * * * * * * * * * РАБОТА ПРОЦЕССА * * * * * * * * * * */

Proc::Proc() {
    stat.timers["total"] = MpiTimer();
    stat.timers["total"].Start();
    stat.timers["io"] = MpiTimer();
    stat.timers["build_ghosts"] = MpiTimer();
    stat.timers["communication"] = MpiTimer();
    stat.timers["find_cell"] = MpiTimer();
    stat.timers["step"] = MpiTimer();
    stat.timers["get_border_cond"] = MpiTimer();
    stat.timers["get_possible_neighs"] = MpiTimer();
    stat.timers["sort_neighs"] = MpiTimer();
}

Proc::~Proc() {

    if (ghosts_in != nullptr) {
        delete[] ghosts_in;
    }
    if (ghosts_out != nullptr) {
        delete[] ghosts_out;
    }

    stat.timers["total"].Stop();
    std::cout << mpiInfo.comm_rank << " ";
    std::cout << stat.toString() << std::endl;
}

int Proc::MPIInit(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    mpiInfo.comm = MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiInfo.comm_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiInfo.comm_size);

    return 0;
}
int Proc::MPIFinalize() {
    MPI_Finalize();
    return 0;
}

int Proc::InitMesh() {
    stat.timers["init_mesh"] = MpiTimer();
    stat.timers["init_mesh"].Start();

    totalG = base_sz * base_sz;
    procG  = totalG / mpiInfo.comm_size;
    offsetG = procG * mpiInfo.comm_rank;

    int global_i_start = 0;
    int global_i_stop = 1 << (max_lvl*2);
    int global_i_step = 1 << (2*max_lvl - 2*base_lvl);

    int my_i_start = offsetG * global_i_step;
    int my_i_stop = my_i_start + procG * (global_i_step);

    std::cout << mpiInfo.comm_rank <<  
        " global_i_start=" << global_i_start <<
        " global_i_stop=" << global_i_stop <<
        " global_i_step=" << global_i_step <<
        " my_i_start=" << my_i_start <<
        " my_i_stop=" << my_i_stop << std::endl;

    for (int i = my_i_start; i < my_i_stop; i += global_i_step) {
        mesh.cells.push_back(Cell(base_lvl, i));
    }

    // fill metaInfo (пока вообще обмен не нужен, но если загружать сетку извне, то будет нужно)
    

    // todo set tau using min dx

    stat.timers["init_mesh"].Stop();
    std::cout << "Mesh inited\n";
    return 0;
}

int Proc::InitMesh(string offsets_filename, string cells_filename) {
    stat.timers["init_mesh"] = MpiTimer();
    stat.timers["init_mesh"].Start();

    // [0] - offset, [1] - len
    int range[2];
    int one_sz = 3 * sizeof(int) + sizeof(double);

    MPI_File fh;
    MPI_File_open( mpiInfo.comm, offsets_filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_at(fh, mpiInfo.comm_rank * sizeof(int), &range, 2, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    std::cout << mpiInfo.comm_rank << " read range[0]=" << range[0] << " range[1]=" << range[1] << std::endl;

    vector<char> buffer(range[1], 1);
    std::cout << "will read cells file\n";

    MPI_File_open( mpiInfo.comm, cells_filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_at(fh, range[0], &buffer[0], range[1]-1, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    std::cout << "read cells file buffer.size()=" << buffer.size() << std::endl;
    cout << mpiInfo.comm_rank << " buf[0]=" << int(buffer[0]) << endl;

    mesh.GenFromWriteStruct(buffer);

    stat.timers["init_mesh"].Stop();
    std::cout << mpiInfo.comm_rank <<  " MESH INITED   " <<
                 "n_of_cells=" << mesh.cells.size() <<
                " cells_offset=" << range[0] / one_sz << std::endl;
    

    MetaInfo my_meta;
    my_meta.procStart = mesh.cells[0].get_global_number();
    my_meta.procEnd   = mesh.cells[mesh.cells.size()-1].get_global_number();
    my_meta.procG     = mesh.cells.size();

    cout << mpiInfo.comm_rank << " my_meta={" << my_meta.procStart << "," << my_meta.procEnd << "," << my_meta.procG << "}\n";

    meta = FullMeta(my_meta, mpiInfo.comm_size, mpiInfo.comm_rank);


    // cout << mpiInfo.comm_rank << " sizeof(my_meta)=" << sizeof(MetaInfo) << " sizeof(meta)=" << meta.size() * sizeof(meta) << endl;
    MPI_Allgather(&meta.metas[mpiInfo.comm_rank*3], meta.one_meta_len, MPI_CHAR, &meta.metas[0], meta.one_meta_len, MPI_CHAR, mpiInfo.comm);

    std::cout << mpiInfo.comm_rank <<  " exchanged meta  ";
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        std::cout << meta.GetMetaOfProc(i).toString() << " | ";
    }
    
    return 0;
}


void Proc::MarkToRefine() {
    // not implemented
}

int Proc::Refine() {
    // not implemented
    return 0;
}

int Proc::LoadBalance() {
    // not implemented
    return 0;
}


int Proc::BuildGhosts() {
    stat.timers["build_ghosts"].Start();
    // TODO 1
    // нужна структура, чтобы знать, где искать гостов
    ghosts_in = new LinearTree[mpiInfo.comm_size];
    ghosts_out = new LinearTree[mpiInfo.comm_size];

    for (Cell c: mesh.cells) {
        vector<Cell> neighs;
        
        vector<GlobalNumber_t> neigh_ids = c.get_all_possible_neighbours_ids();

        // find their owners
        for (GlobalNumber_t neigh: neigh_ids) {
            int owner = find_owner(neigh);
            if (owner == -1) {
                // тогда ячейки нет
                continue;
            }
            if (owner != mpiInfo.comm_rank) {
                // std::cout << "owner not me\n";
                // ghosts_in[owner].cells.resize(ghosts_in[owner].cells.size());
                // ghosts_out[owner].cells.push_back(c); // todo: удалить повторения
            }
        }
    }

    std::cout << mpiInfo.comm_rank << " ghosts built\n";
    stat.timers["build_ghosts"].Stop();
    return 0;
}

int Proc::find_owner(GlobalNumber_t cell_id) {
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        if (meta.GetMetaOfProc(i).procStart > cell_id) {
            return (i-1);
        }
    }
    if (meta.GetMetaOfProc(mpiInfo.comm_size).procStart <= cell_id) {
        return mpiInfo.comm_size;
    }
    // std::cout << "owner not found for cell_id=" << cell_id << std::endl;
    return -1;
}



void Proc::FillStart(double (*start_func)(double, double)) {
    double x, y;
    for (int i = 0; i < mesh.cells.size(); i++) {
        mesh.cells[i].get_spacial_coords(&x, &y);
        
        mesh.cells[i].temp[0] = start_func(x, y);
        std::cout << mpiInfo.comm_rank << " cell(" << mesh.cells[i].lvl << ", " << mesh.cells[i].i << ", " << mesh.cells[i].j << "): (" << x << ", " << y << ") --> " << mesh.cells[i].temp[0] << "\n";
    }

    std::cout << mpiInfo.comm_rank << " filled start at temp[0]\n";
}
void Proc::MakeStep() {
    stat.timers["step"].Start();
    int temp_l_corr = time_step_n % 2; // чтобы брать значение temp[0] или temp[1]
    int cur_temp_idx = temp_l_corr;
    int next_temp_idx = (temp_l_corr + 1) % 2;
    time_step_n++;

    // todo exchange ghosts

    for (int i = 0; i < mesh.cells.size(); i++) {
        char border_cond_type;
        double (*cond_func)(double, double, double);
        Cell cell = mesh.cells[i];
        stat.timers["get_border_cond"].Start();
        cell.get_border_cond(&border_cond_type, &cond_func);
        stat.timers["get_border_cond"].Stop();
        double x, y;
        cell.get_spacial_coords(&x, &y);

        double new_T = 0;

        if (border_cond_type == -1) {  // внутренняя ячейка, но мб не у нас
            
            vector<double> termo_flows;

            stat.timers["get_possible_neighs"].Start();
            vector<GlobalNumber_t> possible_neigh_ids = cell.get_all_possible_neighbours_ids();
            stat.timers["get_possible_neighs"].Stop();
            for (GlobalNumber_t id: possible_neigh_ids) {
                int owner = find_owner(id);
                if (owner == -1) {
                    continue;
                }
                Cell neigh_cell;

                if (owner == mpiInfo.comm_rank) { // ячейка у меня
                    
                    stat.timers["find_cell"].Start();
                    int status = mesh.FindCell(id, &neigh_cell);
                    stat.timers["find_cell"].Stop();
                    if (status == -1) {
                        continue;
                    }
                } else { // ячейка-призрак
                    // int status = ghosts_in[owner].FindCell(id, &neigh_cell);
                    // if (status == -1) {
                    //     continue;
                    // }
                    continue;
                }

                double neigh_x, neigh_y;
                neigh_cell.get_spacial_coords(&neigh_x, &neigh_y);
                double flow = (neigh_cell.temp[cur_temp_idx] - cell.temp[cur_temp_idx]) / dist(x, y, neigh_x, neigh_y);
                if (cell.lvl != neigh_cell.lvl) {
                    flow *= 0.3 * sqrt(10);
                }
                termo_flows.push_back(flow);
            }

            double flows_sum = 0.0;
            for (double f: termo_flows) {
                flows_sum += f;
            }

            new_T = cell.temp[cur_temp_idx] + flows_sum / cell.get_S() + Area::Q(x, y, tau * time_step_n);

        } else if (border_cond_type == 1) {
            new_T = cond_func(x, y, time_step_n * tau); // ? в центрах ячеек по-другому считается?
        } else if (border_cond_type == 2) {
            std::cout << "border cond type = 2 not implemented\n";
            // not imptemented
        }

        mesh.cells[i].temp[next_temp_idx] = new_T;
    }
    stat.timers["step"].Stop();
}

void Proc::WriteT(string filename) {
    stat.timers["io"].Start();

    vector<char> buf = mesh.GenWriteStruct();
    int len = buf.size();
    int *lens = new int[mpiInfo.comm_size];
    stat.timers["communication"].Start();
    MPI_Allgather(&len, 1, MPI_INT, (void *)lens, 1, MPI_INT, mpiInfo.comm);
    stat.timers["communication"].Stop();
    int offset = 0;
    for (int i = 0; i < mpiInfo.comm_rank; i++) {
        offset += lens[i];
    }

    std::cout << mpiInfo.comm_rank << " offset=" << offset << " wr_len=" << len << std::endl;

    MPI_File fh;
    MPI_File_open( mpiInfo.comm, 
                filename.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_at(fh, offset, &buf[0], len, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    stat.timers["io"].Stop();
}

// vector<char> Proc::get_write_cells_struct() {
//     vector<char> buf;
//     std::cout << mpiInfo.comm_rank << " gen structs for write start\n";
//     for (auto c: mesh.cells) {
//         // int local_sz = sizeof(int) + sizeof(int) + sizeof(int) + sizeof(double); // lvl, i, j, temp[0]
//         // buf.resize(buf.size() + local_sz);
//         char *tmp = (char *)(&c.lvl);
//         for (int i = 0; i < sizeof(int); i++) {
//             buf.push_back(tmp[i]);
//         }
//         tmp = (char *)(&c.i);
//         for (int i = 0; i < sizeof(int); i++) {
//             buf.push_back(tmp[i]);
//         }
//         tmp = (char *)(&c.j);
//         for (int i = 0; i < sizeof(int); i++) {
//             buf.push_back(tmp[i]);
//         }
//         if (c.temp[0] > 1) {
//             std::cout << (c.temp[0] > 5) << std::endl;
//         }
        
//         tmp = (char *)(&c.temp[0]);
//         for (int i = 0; i < sizeof(double); i++) {
//             buf.push_back(tmp[i]);
//         }
//     }
//     std::cout << mpiInfo.comm_rank << " gen structs for write finished\n";
//     return buf;
// }

void Proc::WriteStat(string filename) {
    std::cout << "MpiInfo: " << mpiInfo.toString() << std::endl;
}

int Proc::ISendGhosts() {
    // MPI_ISend()

    return 0;
}
int Proc::IRecvGhosts() {
    // MPI_IRecv()
    // запомнить шутку, по которой ждать потом
    return 0;
}
int Proc::WaitallGhosts() {
    return 0;
}

// void MPI_Allgather_wrapper()
