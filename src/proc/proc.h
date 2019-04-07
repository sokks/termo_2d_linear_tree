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

extern double base_tau;
extern double tau;
extern int    time_steps;

void Init(int gs_x, int ts_n);


/* * * * * * * * * * * СТАТИСТИКА * * * * * * * * * * */

class MpiTimer {
    double s_time;
    double dur = 0;

public:
    void   Start() { s_time = MPI_Wtime(); }
    double Stop() { double cur_dur = MPI_Wtime() - s_time; dur += cur_dur; return cur_dur; }
    double FullDur() { return dur; }
};

struct MpiInfo {
    MPI_Comm comm;
    int comm_rank;
    int comm_size;

    string toString();
};

struct Stat {
    int proc_cells;
    map<string, MpiTimer> timers;

    string toString();
};



/* * * * * * * * * * * РАБОТА ПРОЦЕССА * * * * * * * * * * */

struct MetaInfo {
    GlobalNumber_t procStart;
    GlobalNumber_t procEnd;
    GlobalNumber_t procG;

    string toString();
};

class Proc {
    MpiInfo mpiInfo;
    Stat    stat;

    int time_step_n = 0;

    
    GlobalNumber_t totalG;
    GlobalNumber_t procG;
    GlobalNumber_t offsetG;
    MetaInfo *meta = nullptr;
    LinearTree tree;

    LinearTree *ghosts_in = nullptr;
    LinearTree *ghosts_out = nullptr;

public:
    Proc();
    ~Proc();

    /// MPIInit must be called before start work!
    int MPIInit(int argc, char **argv);
    int MPIFinalize();

    /// InitMesh создает базовую структуру сетки
    int InitMesh();
    int InitMesh(string filename);

    /// MarkToRefine
    void MarkToRefine();

    // Refine измельчает помеченные ячейки
    //
    // Оперирует с линейным деревом, добавляя элементы-ячейки 
    int Refine();

    // LoadBalance балансирует нагрузку (линейное дерево) по процессорам
    int LoadBalance();

    // BuildGhosts строит структуры для обмена границами
    int BuildGhosts();

    // FillStart заполняет начальные заначения температуры
    void FillStart(double (*start_func)(double, double));

    // MakeStep делает один временной шаг
    void MakeStep();

    // I/O
    void WriteT(string filename);
    void WriteStat(string filename);

private:

    int find_owner(GlobalNumber_t cell_id);
    vector<char> get_write_cells_struct();

    int ISendGhosts();
    int IRecvGhosts();
    int WaitallGhosts();

// запрещаем копирование
    Proc(Proc&);
    Proc& operator=(Proc&);
};

