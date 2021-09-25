#include <stdlib.h> // malloc
#include <stdio.h> // printf
#include <iostream> // cout
#include <cmath> // math lib
#include <fstream> // file I/O
#include <iomanip> 
#include<string>
using namespace std;

void visualization_planex(double** gen, int nx, int ny) {
    cout << setprecision(4) << fixed;

    for (int j = ny - 1; j >= 0; j--) {
        for (int i = 0; i <= nx - 1; i++) {
            cout << gen[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n\n";
}



void initualization(double** u, double** u_new, double** v, double** v_new, double** p, double** p_new, double** F, double** G, double** p_rhs, int nx, int ny) {
    for (int i = 0; i <= nx - 1; i++) {
        for (int j = 0; j <= ny - 1; j++) {
            u[i][j] = 0;
            u_new[i][j] = 0;
            v[i][j] = 0;
            v_new[i][j] = 0;
            p[i][j] = 0;
            p_new[i][j] = 0;
            F[i][j] = 0;
            G[i][j] = 0;
            p_rhs[i][j] = 0;

        }
    }
}


void setbound(double** u_new, double** v_new, int nx, int ny) {

    for (int i = 0; i <= nx - 1; i++) {
        u_new[i][0] = -u_new[i][1];
        v_new[i][0] = 0;

        u_new[i][ny - 1] = -u_new[i][ny - 2];
        v_new[i][ny - 2] = 0;


    }

    for (int j = 0; j <= ny - 1; j++) {
        u_new[0][j] = 1;
        v_new[0][j] = -v_new[1][j];

        u_new[nx - 2][j] = u_new[nx - 3][j];
        v_new[nx - 1][j] = v_new[nx - 2][j];
    }

}


void alittlestepforhumankind(double** u, double** u_new, double** v, double** v_new, int nx, int ny) {

    for (int i = 0; i <= nx - 1; i++) {
        for (int j = 0; j <= ny - 1; j++) {

            u[i][j] = u_new[i][j];
            v[i][j] = v_new[i][j];

        }
    }

}

void comp_FG(double** u, double** v, double** F, double** G, double dx, double dy, double donor, double Re, double gx, double gy, double ts, int nx, int ny) {

    ///////////////////////////only from 1 to n-2////////////////////////////
    for (int j = 1; j <= ny - 2; j++) {
        F[0][j] = u[0][j];
        F[nx - 2][j] = u[nx - 2][j];

    }

    for (int i = 1; i <= nx - 2; i++) {

        G[i][0] = v[i][0];
        G[i][ny - 2] = v[i][ny - 2];

    }



    //////////////////////////////n-3 confusing U deee***********8/////////////
    for (int i = 1; i <= nx - 3; i++) {//change
        for (int j = 1; j <= ny - 2; j++) {

            double duu_dx = ((pow(((u[i][j] + u[i + 1][j]) / 2.0), 2) - pow(((u[i - 1][j] + u[i][j]) / 2.0), 2)) / dx) + (donor / dx) * (abs(u[i][j] + u[i + 1][j]) * (u[i][j] - u[i + 1][j]) / 4.0 - abs(u[i - 1][j] + u[i][j]) * (u[i - 1][j] - u[i][j]) / 4.0);

            double duv_dy = (((v[i][j] + v[i + 1][j]) * (u[i][j] + u[i][j + 1]) / 4.0 - (v[i][j - 1] + v[i + 1][j - 1]) * (u[i][j - 1] + u[i][j]) / 4.0) / dy) + (donor / dy) * (abs(v[i][j] + v[i + 1][j]) * (u[i][j] - u[i][j + 1]) / 4.0 - abs(v[i][j - 1] + v[i + 1][j - 1]) * (u[i][j - 1] - u[i][j]) / 4.0);

            double ddu_dxdx = (u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) / pow(dx, 2);
            double ddu_dydy = (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / pow(dy, 2);

            //cout << "f/re = " << (ddu_dxdx+ddu_dydy+ddu_dzdz)/Re  <<"\n"; //- (ts/dx)*(p_new[i+1][j][k] - p_new[i][j][k]) << "\n";
            F[i][j] = u[i][j] + ts * ((ddu_dxdx + ddu_dydy) / Re - duu_dx - duv_dy + gx);

        }
    }


    //---------------------------------------------for G


    for (int j = 1; j <= ny - 3; j++) {//change
        for (int i = 1; i <= nx - 2; i++) {


            double dvv_dy = ((pow(((v[i][j] + v[i][j + 1]) / 2.0), 2) - pow(((v[i][j - 1] + v[i][j]) / 2.0), 2)) / dy) + (donor / dy) * (abs(v[i][j] + v[i][j + 1]) * (v[i][j] - v[i][j + 1]) / 4.0 - abs(v[i][j - 1] + v[i][j]) * (v[i][j - 1] - v[i][j]) / 4.0);

            double duv_dx = (((v[i][j] + v[i + 1][j]) * (u[i][j] + u[i][j + 1]) / 4.0 - (v[i - 1][j] + v[i][j]) * (u[i - 1][j + 1] + u[i - 1][j]) / 4.0) / dx) + (donor / dx) * (abs(u[i][j] + u[i][j + 1]) * (v[i][j] - v[i + 1][j]) / 4.0 - abs(u[i - 1][j + 1] + u[i - 1][j]) * (v[i - 1][j] - v[i][j]) / 4.0);

            double ddv_dxdx = (v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j]) / pow(dx, 2);
            double ddv_dydy = (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]) / pow(dy, 2);

            G[i][j] = v[i][j] + ts * ((ddv_dxdx + ddv_dydy) / Re - dvv_dy - duv_dx + gy);

        }
    }


    //////////////////////////////n-3 confusing U deee***********8/////////////
}


void comp_p_rhs(double** F, double** G, double** p_rhs, double ts, double dx, double dy, int nx, int ny) {
    //cout << "rhs" << "\n";
    for (int i = 1; i <= nx - 2; i++) {
        for (int j = 1; j <= ny - 2; j++) {

            p_rhs[i][j] = ((F[i][j] - F[i - 1][j]) / dx + (G[i][j] - G[i][j - 1]) / dy) / ts;

        }
    }
}
/////////////////////////////cal only newly created p's////////////////////
double find_max(double** mat, int nx, int ny) {
    double maxElement = 0.;
    for (int i = 1; i <= nx - 2; i++) {
        for (int j = 1; j <= ny - 2; j++) {

            if (abs(mat[i][j]) > maxElement) {
                maxElement = abs(mat[i][j]);
            }

        }
    }

    // finally return maxElement
    return maxElement;
}
/////////////////////////////cal only newly created p's////////////////////

void comp_p(double** p, double** p_new, double** p_rhs, double** res, double dx, double dy, int nx, int ny, int it_max, double eps) {
    int it = 1;
    double r = 9999;

    /////////////////update p's periodic bc every iteration/////////////////
    while (it <= it_max && abs(r) >= eps) {
        for (int i = 1; i <= nx - 2; i++) {

            p_new[i][0] = p[i][1];
            p_new[i][ny - 1] = p[i][ny - 2]; // presure strip for north and south bound

        }
        for (int j = 1; j <= ny - 2; j++) {

            p_new[0][j] = p[1][j];
            p_new[nx - 1][j] = -p_new[nx - 2][j]; // pressure strip for east and west bound


        }

        /////////////////update p's periodic bc every iteration/////////////////


        double sum_rsqr = 0;
        //cout << "pressure_lang_up_kob" << "\n";
        //visualization_gen(p_new,nx,ny);
        for (int i = 1; i <= nx - 2; i++) {//////////
            for (int j = 1; j <= ny - 2; j++) {////////
                p_new[i][j] = (-0.7) * p[i][j] + (1.7 / (2 / pow(dx, 2) + 2 / pow(dy, 2))) * ((p[i + 1][j] + p_new[i - 1][j]) / pow(dx, 2) + (p[i][j + 1] + p_new[i][j - 1]) / pow(dy, 2) - p_rhs[i][j]);
            }
        }

        ///////////////////// 2 to n-2 (1 is already updated)/////////////////////////


        //cout << "pressure_lang_up_value" << "\n";
        //visualization_gen(p_new,nx,ny);

        /*
        for ( int i = 1; i<=ny-2; i++){
          for ( int j = 1; j<=nx-2; j++){
        double rsqr = pow(((p_new[i][j+1] - 2*p_new[i][j] + p_new[i][j-1])/pow(dx,2) + (p_new[i-1][j] - 2*p_new[i][j] + p_new[i+1][j])/pow(dy,2) - p_rhs[i][j]),2);
        sum_rsqr = sum_rsqr + rsqr;      F[0][j][k] = u[0][j][k];
        cout << rsqr << " ";
          }
          cout << "\n";
        }
        r = pow((sum_rsqr/ ((nx-2)*(ny-2))),0.5);
        */
        //cout << "res" << "\n";



        for (int i = 1; i <= nx - 2; i++) {
            for (int j = 1; j <= ny - 2; j++) {

                res[i][j] = ((p_new[i + 1][j] - p_new[i][j]) - (p_new[i][j] - p_new[i - 1][j])) / (dx * dx)
                    + ((p_new[i][j + 1] - p_new[i][j]) - (p_new[i][j] - p_new[i][j - 1])) / (dy * dy) - p_rhs[i][j];

            }
        }

        //visualization_gen(res,nx,ny);

        //cout << "\n";
        r = find_max(res, nx, ny);

        it++;
        //cout << "residual = " << r << "\n" << "interation =" << it << "\n";
        //update p
        for (int i = 0; i <= nx - 1; i++) {
            for (int j = 0; j <= ny - 1; j++) {

                p[i][j] = p_new[i][j];

            }
        }
    }
    cout << "residual = " << abs(r) << "\n" << "interation =" << it << "\n";



}


/////////////////////////////// sperate cal (follow the book) ////////////////////////////////////
void finalcomp(double** u_new, double** v_new, double** F, double** G, double** p_new, double ts, double dx, double dy, int nx, int ny)
{
    for (int i = 1; i <= nx - 3; i++) {
        for (int j = 1; j <= ny - 2; j++) {

            u_new[i][j] = F[i][j] - (ts / dx) * (p_new[i + 1][j] - p_new[i][j]);

        }
    }
    for (int j = 1; j <= ny - 3; j++) {//change
        for (int i = 1; i <= nx - 2; i++) {

            v_new[i][j] = G[i][j] - (ts / dy) * (p_new[i][j + 1] - p_new[i][j]);

        }
    }

    /////////////////////////////// sperate cal (follow the book) ////////////////////////////////////
}
void write_result(double** var, double dx, double dy, int nx, int ny, int ts, string x) {
    //var = simulated data (array)
    //nx, ny = size (2d)
    //ts = timestep
    ofstream myfile;
    string file_name = "vtk_" + x + "/TwoD_" + to_string(ts) + ".vtk";
    myfile.open(file_name);
    // Paraview stuff
    //Header
    myfile << "# vtk DataFile Version 2.0\n";
    myfile << "FlowField\n";
    myfile << "ASCII\n";
    // Grid
    myfile << "DATASET STRUCTURED_GRID\n";
    myfile << "DIMENSIONS " << nx << " " << ny << " " << 1 << "\n";
    myfile << "POINTS " << nx * (ny) << " float\n";
    for (int j = 0; j <= ny - 1; j++) {
        for (int i = 0; i <= nx - 1; i++) {
            myfile << i * dx << " " << j * dy << " 0" << endl;
        }
    }
    // Dataset
    myfile << "\n";
    myfile << "POINT_DATA";
    myfile << " " << nx * (ny) << "\n";
    // Point data
    myfile << "\n";
    myfile << "SCALARS " + x + " float 1\n";
    myfile << "LOOKUP_TABLE default\n";
    for (int j = 0; j <= ny - 1; j++) {
        for (int i = 0; i <= nx - 1; i++) {
            myfile << var[i][j] << "\n";
        }
    }
    // end
    myfile.close();
}
int main() {

    int nx = 100;
    int ny = 31;

    double x_size = 10.8;//streamwise
    double y_size = 0.9;//spanwise

    double dx = x_size / (nx - 2);
    double dy = y_size / (ny - 2);

    double donor = 0;
    double gx = 0;
    double gy = 0;

    int t_max = 5000;
    int it_max = 1000; //for poisson iteration
    double eps = 0.01; //for poisson iteration



    //double u_in = 0;
    double Re = 50;//1500;

    //double ts = Re/((1/pow(dx,2))+(1/pow(dy,2)))/2*1.1; // choose ts based on 3.5 with a sf of 1.1
    double ts = 0.007;


    ////////////////////////////////////////////////////CREATE ARRAYS///////////////////////////////////////////

    double** u = (double**)malloc(nx * sizeof(double*));
    for (int i = 0; i <= nx - 1; i++) {
        u[i] = (double*)malloc(ny * sizeof(double));

    }

    double** u_new = (double**)malloc(nx * sizeof(double*));
    for (int i = 0; i <= nx - 1; i++) {
        u_new[i] = (double*)malloc(ny * sizeof(double));

    }

    double** v = (double**)malloc(nx * sizeof(double*));
    for (int i = 0; i <= nx - 1; i++) {
        v[i] = (double*)malloc(ny * sizeof(double));

    }

    double** v_new = (double**)malloc(nx * sizeof(double*));
    for (int i = 0; i <= nx - 1; i++) {
        v_new[i] = (double*)malloc(ny * sizeof(double));

    }


    double** p = (double**)malloc(nx * sizeof(double*));
    for (int i = 0; i <= nx - 1; i++) {
        p[i] = (double*)malloc(ny * sizeof(double));

    }

    double** p_new = (double**)malloc(nx * sizeof(double*));
    for (int i = 0; i <= nx - 1; i++) {
        p_new[i] = (double*)malloc(ny * sizeof(double));

    }

    double** F = (double**)malloc(nx * sizeof(double*));
    for (int i = 0; i <= nx - 1; i++) {
        F[i] = (double*)malloc(ny * sizeof(double));

    }

    double** G = (double**)malloc(nx * sizeof(double*));
    for (int i = 0; i <= nx - 1; i++) {
        G[i] = (double*)malloc(ny * sizeof(double));

    }

    double** p_rhs = (double**)malloc(nx * sizeof(double*));
    for (int i = 0; i <= nx - 1; i++) {
        p_rhs[i] = (double*)malloc(ny * sizeof(double));

    }
    double** res = (double**)malloc(nx * sizeof(double*));
    for (int i = 0; i <= nx - 1; i++) {
        res[i] = (double*)malloc(ny * sizeof(double));
    }

    initualization(u, u_new, v, v_new, p, p_new, F, G, p_rhs, nx, ny);

    setbound(u_new, v_new, nx, ny);

    alittlestepforhumankind(u, u_new, v, v_new, nx, ny);

    //showinipls ( u, v, p, nx, ny);

    // ptr=&u[0][0];



    for (int t = 0; t <= t_max - 1; t++) {
        cout << "t= " << t << "\n";
        comp_FG(u, v, F, G, dx, dy, donor, Re, gx, gy, ts, nx, ny);

        comp_p_rhs(F, G, p_rhs, ts, dx, dy, nx, ny);

        comp_p(p, p_new, p_rhs, res, dx, dy, nx, ny, it_max, eps);

        finalcomp(u_new, v_new, F, G, p_new, ts, dx, dy, nx, ny);

        setbound(u_new, v_new, nx, ny);
        alittlestepforhumankind(u, u_new, v, v_new, nx, ny);
       
        write_result(u_new, dx, dy, nx, ny, t, "u");
        write_result(p_new, dx, dy, nx, ny, t, "p");
        write_result(v_new, dx, dy, nx, ny, t, "v");
        cout << "\n\n";
    }

}
