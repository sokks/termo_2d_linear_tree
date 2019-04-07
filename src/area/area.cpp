#include "area.h"

// описание области
double Area::x_start = 0.0;
double Area::y_start = 0.0;
double Area::x_end = 1.0;
double Area::y_end = 1.0;

// double ro = 1190.0, lambda = 0.16, c = 900.0;
double Area::ro = 1, Area::lambda = 0.01, Area::c = 1;
double Area::a = lambda / (ro * c);
// double a = 1;



// начальное распределение температуры и граничные условия
double Area::T0(double x, double y) {
    double r = 0.3;
    double c = 0.5;
    if ( (x-c) * (x-c) + (y-c) * (y-c) < r*r ) { //,0078125 < 0,09
        return 10;
    }
    return 0.0;
}

double Area::BorderCond1(double x, double y, double t) { // down & left
    return 0.0;
}

double Area::BorderCond2(double x, double y, double t) { // right
    return -2 * t;
}

double Area::BorderCond3(double x, double y, double t) { // up
    return 100 * t;
}

void Area::get_border_cond(Border b, char *cond_type, double (**cond_func)(double, double, double)) {

    switch (b) {
        case Border::LEFT: {
            *cond_type = 1;
            *cond_func = &BorderCond1;
            break;
        }
        case Border::RIGHT: {
            // *cond_type = 2;
            // *cond_func = &BorderCond2;
            *cond_type = 1;
            *cond_func = &BorderCond1;
            break;
        }
        case Border::UP: {
            // *cond_type = 2;
            // *cond_func = &BorderCond3;
            *cond_type = 1;
            *cond_func = &BorderCond1;
            break;
        }
        case Border::DOWN: {
            *cond_type = 1;
            *cond_func = &BorderCond1;
            break;
        }
    }
}



double Area::Q(double x, double y, double t) {
    if (t > 0.2) {
        return 0;
    }
    double eps = 0.01;
    if ((fabs(x - 0.45) < eps) && (fabs(y - 0.5) < eps)) {
        return 0.01;
    }
    if ((fabs(x - 0.55) < eps) && (fabs(y - 0.5) < eps)) {
        return 0.01;
    }
    return 0;
}


int Area::max_grid_levels = 2;
double Area::refine_thresholds[] = {10};
