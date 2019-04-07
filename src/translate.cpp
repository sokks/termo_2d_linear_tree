#include <iostream>
#include <fstream>
#include <vector>

using namespace std;


struct WriteCell {
    int lvl, i, j;
    double temp;
};

void copy_buf(char *from_buf, int from_start, int from_end, 
              char *to_buf,   int to_start,   int to_end);
void write_txt_cells(vector<WriteCell> cells);
vector<WriteCell> to_write_cells(vector<char> buf);

ofstream fout;


int main(int argc, char **argv) {
    if (argc < 3) {
        cout << "not enough arguments\n";
    }
    string filename_in  = argv[1];
    string filename_out = argv[2];

    ifstream fin(filename_in, ios::binary);

    streampos fileSize;
    fin.seekg(0, std::ios::end);
    fileSize = fin.tellg();
    std::vector<char> buffer(fileSize);

    fin.seekg(0, std::ios::beg);
    fin.read((char*) &buffer[0], fileSize);

    vector<WriteCell> cells = to_write_cells(buffer);

    fout.open(filename_out);
    write_txt_cells(cells);

    fin.close();
    fout.close();
    return 0;
}

vector<WriteCell> to_write_cells(vector<char> buf) {
    char *p = &buf[0];
    int pos = 0;

    int one_sz = 3 * sizeof(int) + sizeof(double);

    int lvl_offset = 0;
    int i_offset = sizeof(int);
    int j_offset = 2 * sizeof(int);
    int temp_offset = 3 * sizeof(int);

    vector<WriteCell> ret;
    while (pos < buf.size()) {
        WriteCell c;
        c.lvl  = * ((int *)(&p[pos+lvl_offset]));
        c.i    = * ((int *)(&p[pos+i_offset]));
        c.j    = * ((int *)(&p[pos+j_offset]));
        c.temp = * ((double *)(&(p[pos+temp_offset])));
        ret.push_back(c);

        pos += one_sz;
    }
    
    return ret;
}

void write_txt_cells(vector<WriteCell> cells) {
    for (WriteCell c: cells) {
        fout << c.lvl << "," << c.i << "," << c.j << "," << c.temp << endl;
    }
}


