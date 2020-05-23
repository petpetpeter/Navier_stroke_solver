#include <stdlib.h> // malloc
#include <stdio.h> // printf
#include <iostream> // cout
#include <cmath> // math lib
#include <fstream> // file I/O
#include <iomanip> //forsomething
#include<string>
using namespace std;

void visualization_planex(double*** gen, int nx, int ny, int nz) {// later
    cout << setprecision(4) << fixed;

    for (int k = nz - 1; k >= 0; k--) {
        for (int i = 0; i <= nx - 1; i++) {
            cout << gen[i][ny / 2][k] << " ";
        }
        cout << "\n";
    }
    cout << "\n\n";
}


void initualization(double*** u, double*** u_new, double*** v, double*** v_new, double*** w, double*** w_new, double*** p, double*** p_new, double*** F, double*** G, double*** H, double*** p_rhs, int nx, int ny, int nz) {
    for (int i = 0; i <= nx - 1; i++) {
        for (int j = 0; j <= ny - 1; j++) {
            for (int k = 0; k <= nz - 1; k++) {
                u[i][j][k] = 0;
                u_new[i][j][k] = 0;
                v[i][j][k] = 0;
                v_new[i][j][k] = 0;
                w[i][j][k] = 0;
                w_new[i][j][k] = 0;
                p[i][j][k] = 0;
                p_new[i][j][k] = 0;
                F[i][j][k] = 0;
                G[i][j][k] = 0;
                H[i][j][k] = 0;
                p_rhs[i][j][k] = 0;
            }
        }
    }
}


void setbound(double*** u_new, double*** v_new, double*** w_new, double*** p_new, double*** p, int nx, int ny, int nz) {
    for (int j = 0; j <= ny - 1; j++) {
        for (int k = 1; k <= nz - 2; k++) {

            u_new[0][j][k] = u_new[nx - 3][j][k];
            v_new[0][j][k] = v_new[nx - 3][j][k];
            w_new[0][j][k] = w_new[nx - 3][j][k];
            v_new[1][j][k] = v_new[nx - 2][j][k];
            w_new[1][j][k] = w_new[nx - 2][j][k];
            //p_new[1][j][k] = p_new[nx-2][j][k];// periodic  condition for west bound.

            u_new[nx - 2][j][k] = u_new[1][j][k];
            v_new[nx - 1][j][k] = v_new[2][j][k];
            w_new[nx - 1][j][k] = w_new[2][j][k];// periodic (copy) for east bound

        }
    }

    for (int i = 0; i <= nx - 1; i++) {
        for (int j = 0; j <= ny - 1; j++) {

            u_new[i][j][nz - 1] = -u_new[i][j][nz - 2];
            v_new[i][j][nz - 1] = -v_new[i][j][nz - 2];
            w_new[i][j][nz - 2] = 0;// no slip for north bound

            u_new[i][j][0] = -u_new[i][j][1];
            v_new[i][j][0] = -v_new[i][j][0];
            w_new[i][j][0] = 0;// no slip for south bound

        }
    }

    for (int i = 0; i <= nx - 1; i++) {
        for (int k = 1; k <= nz - 2; k++) {

            u_new[i][0][k] = u_new[i][ny - 3][k];
            u_new[i][1][k] = u_new[i][ny - 2][k];
            v_new[i][0][k] = v_new[i][ny - 3][k];
            w_new[i][0][k] = w_new[i][ny - 3][k];
            w_new[i][1][k] = w_new[i][ny - 2][k];
            //p_new[i][1][k] = p_new[i][ny-2][k];// periodic (copy) for outward bound

            u_new[i][ny - 1][k] = u_new[i][2][k];
            v_new[i][ny - 2][k] = v_new[i][1][k];
            w_new[i][ny - 1][k] = w_new[i][2][k];// periodic (copy) for inward bound
        }
    }
    //************** Fix p_new = p_new******************************//
    for (int i = 1; i <= nx - 2; i++) {
        for (int j = 1; j <= ny - 2; j++) {
            p_new[i][j][0] = p_new[i][j][1];
            p_new[i][j][nz - 1] = p_new[i][j][nz - 2]; // presure strip for north and south bound
        }
    }
    for (int j = 1; j <= ny - 2; j++) {
        for (int k = 1; k <= nz - 2; k++) {
            p_new[0][j][k] = p_new[1][j][k];
            p_new[nx - 1][j][k] = p_new[nx - 2][j][k]; // pressure strip for east and west bound
            p_new[1][j][k] = p_new[nx - 2][j][k]; //stream_wise_periodic
        }
    }
    for (int i = 1; i <= nx - 2; i++) {
        for (int k = 1; k <= nz - 2; k++) {
            p_new[i][0][k] = p_new[i][1][k];
            p_new[i][ny - 1][k] = p_new[i][ny - 2][k];
            p_new[i][1][k] = p_new[i][ny - 2][k]; //span_wise_periodic
        }
    }
    //************** Fix p_new = p_new******************************//
}


//*****************update p also`*********************************//
void alittlestepforhumankind(double*** u, double*** u_new, double*** v, double*** v_new, double*** w, double*** w_new, double*** p, double*** p_new, int nx, int ny, int nz) {

    for (int i = 0; i <= nx - 1; i++) {
        for (int j = 0; j <= ny - 1; j++) {
            for (int k = 0; k <= nz - 1; k++) {
                u[i][j][k] = u_new[i][j][k];
                v[i][j][k] = v_new[i][j][k];
                w[i][j][k] = w_new[i][j][k];
                p[i][j][k] = p_new[i][j][k];

            }
        }
    }
    //*****************update p also`*********************************//
}

void comp_FG(double*** u, double*** v, double*** w, double*** F, double*** G, double*** H, double dx, double dy, double dz, double donor, double Re, double gx, double gy, double gz, double ts, int nx, int ny, int nz) {

    ///////////////////////////only from 1 to n-2////////////////////////////
    for (int j = 1; j <= ny - 2; j++) {
        for (int k = 1; k <= nz - 2; k++) {
            F[0][j][k] = u[0][j][k];
            F[nx - 2][j][k] = u[nx - 2][j][k];
        }
    }

    for (int i = 1; i <= nx - 2; i++) {
        for (int k = 1; k <= nz - 2; k++) {
            G[i][0][k] = v[i][0][k];
            G[i][ny - 2][k] = v[i][ny - 2][k];
        }
    }

    for (int j = 1; j <= ny - 2; j++) {
        for (int i = 1; i <= nx - 2; i++) {
            H[i][j][0] = w[i][j][0];
            H[i][j][nz - 2] = w[i][j][nz - 2];
        }
    }
    ///////////////////////////only from 1 to n-2////////////////////////////
      // --------------------------------------------------------for F

    //////////////////////////////n-3 confusing U deee***********8/////////////
    for (int i = 1; i <= nx - 3; i++) {//change
        for (int j = 1; j <= ny - 2; j++) {
            for (int k = 1; k <= nz - 2; k++) {

                double duu_dx = ((pow(((u[i][j][k] + u[i + 1][j][k]) / 2.0), 2) - pow(((u[i - 1][j][k] + u[i][j][k]) / 2.0), 2)) / dx) + (donor / dx) * (abs(u[i][j][k] + u[i + 1][j][k]) * (u[i][j][k] - u[i + 1][j][k]) / 4.0 - abs(u[i - 1][j][k] + u[i][j][k]) * (u[i - 1][j][k] - u[i][j][k]) / 4.0);

                double duv_dy = (((v[i][j][k] + v[i + 1][j][k]) * (u[i][j][k] + u[i][j + 1][k]) / 4.0 - (v[i][j - 1][k] + v[i + 1][j - 1][k]) * (u[i][j - 1][k] + u[i][j][k]) / 4.0) / dy) + (donor / dy) * (abs(v[i][j][k] + v[i + 1][j][k]) * (u[i][j][k] - u[i][j + 1][k]) / 4.0 - abs(v[i][j - 1][k] + v[i + 1][j - 1][k]) * (u[i][j - 1][k] - u[i][j][k]) / 4.0);

                double duw_dz = (((w[i][j][k] + w[i + 1][j][k]) * (u[i][j][k] + u[i][j][k + 1]) / 4.0 - (w[i][j][k - 1] + w[i + 1][j][k - 1]) * (u[i][j][k - 1] + u[i][j][k]) / 4.0) / dz) + (donor / dz) * (abs(w[i][j][k] + w[i + 1][j][k]) * (u[i][j][k] - u[i][j][k + 1]) / 4.0 - abs(w[i][j][k - 1] + w[i + 1][j][k - 1]) * (u[i][j][k - 1] - u[i][j][k]) / 4.0);

                double ddu_dxdx = (u[i + 1][j][k] - 2.0 * u[i][j][k] + u[i - 1][j][k]) / pow(dx, 2);
                double ddu_dydy = (u[i][j + 1][k] - 2.0 * u[i][j][k] + u[i][j - 1][k]) / pow(dy, 2);
                double ddu_dzdz = (u[i][j][k + 1] - 2.0 * u[i][j][k] + u[i][j][k - 1]) / pow(dz, 2);

                //cout << "f/re = " << (ddu_dxdx+ddu_dydy+ddu_dzdz)/Re  <<"\n"; //- (ts/dx)*(p_new[i+1][j][k] - p_new[i][j][k]) << "\n";
                F[i][j][k] = u[i][j][k] + ts * ((ddu_dxdx + ddu_dydy + ddu_dzdz) / Re - duu_dx - duv_dy - duw_dz + gx);
            }
        }
    }


    //---------------------------------------------for G


    for (int j = 1; j <= ny - 3; j++) {//change
        for (int i = 1; i <= nx - 2; i++) {
            for (int k = 1; k <= nz - 2; k++) {

                double dvv_dy = ((pow(((v[i][j][k] + v[i][j + 1][k]) / 2.0), 2) - pow(((v[i][j - 1][k] + v[i][j][k]) / 2.0), 2)) / dy) + (donor / dy) * (abs(v[i][j][k] + v[i][j + 1][k]) * (v[i][j][k] - v[i][j + 1][k]) / 4.0 - abs(v[i][j - 1][k] + v[i][j][k]) * (v[i][j - 1][k] - v[i][j][k]) / 4.0);

                double duv_dx = (((v[i][j][k] + v[i + 1][j][k]) * (u[i][j][k] + u[i][j + 1][k]) / 4.0 - (v[i - 1][j][k] + v[i][j][k]) * (u[i - 1][j + 1][k] + u[i - 1][j][k]) / 4.0) / dx) + (donor / dx) * (abs(u[i][j][k] + u[i][j + 1][k]) * (v[i][j][k] - v[i + 1][j][k]) / 4.0 - abs(u[i - 1][j + 1][k] + u[i - 1][j][k]) * (v[i - 1][j][k] - v[i][j][k]) / 4.0);

                double dvw_dz = (((w[i][j + 1][k] + w[i][j][k]) * (v[i][j][k + 1] + v[i][j][k]) / 4.0 - (w[i][j + 1][k - 1] + w[i][j][k - 1]) * (v[i][j][k - 1] + v[i][j][k]) / 4.0) / dz) + (donor / dz) * (abs(w[i][j + 1][k] + w[i][j][k]) * (v[i][j][k] - v[i][j][k + 1]) / 4.0 - abs(w[i][j + 1][k - 1] + w[i][j][k - 1]) * (v[i][j][k - 1] - v[i][j][k]) / 4.0);

                double ddv_dxdx = (v[i + 1][j][k] - 2.0 * v[i][j][k] + v[i - 1][j][k]) / pow(dx, 2);
                double ddv_dydy = (v[i][j + 1][k] - 2.0 * v[i][j][k] + v[i][j - 1][k]) / pow(dy, 2);
                double ddv_dzdz = (v[i][j][k + 1] - 2.0 * v[i][j][k] + v[i][j][k - 1]) / pow(dz, 2);


                G[i][j][k] = v[i][j][k] + ts * ((ddv_dxdx + ddv_dydy + ddv_dzdz) / Re - dvv_dy - duv_dx - dvw_dz + gy);
            }
        }
    }

    for (int k = 1; k <= nz - 3; k++) {//change
        for (int i = 1; i <= nx - 2; i++) {
            for (int j = 1; j <= ny - 2; j++) {

                double dww_dz = ((pow(((w[i][j][k] + v[i][j][k + 1]) / 2), 2) - pow(((w[i][j][k - 1] + w[i][j][k]) / 2), 2)) / dz) + (donor / dz) * (abs(w[i][j][k] + w[i][j][k + 1]) * (w[i][j][k] - w[i][j][k + 1]) / 4 - abs(w[i][j][k - 1] + w[i][j][k]) * (w[i][j][k - 1] - w[i][j][k]) / 4);

                double duw_dx = (((u[i][j][k] + u[i][j][k + 1]) * (w[i][j][k] + w[i + 1][j][k]) / 4 - (u[i - 1][j][k] + u[i - 1][j][k + 1]) * (w[i - 1][j][k] + w[i][j][k]) / 4) / dx) + (donor / dx) * (abs(u[i][j][k] + u[i][j][k + 1]) * (w[i][j][k] - w[i + 1][j][k]) / 4 - abs(u[i - 1][j][k] + u[i - 1][j][k + 1]) * (w[i - 1][j][k] - w[i][j][k]) / 4);

                double dvw_dy = (((v[i][j][k] + v[i][j][k + 1]) * (w[i][j][k] + w[i][j + 1][k]) / 4 - (v[i][j - 1][k] + v[i][j - 1][k + 1]) * (w[i][j - 1][k] + w[i][j][k]) / 4) / dy) + (donor / dy) * (abs(v[i][j][k] + v[i][j][k + 1]) * (w[i][j][k] - w[i][j + 1][k]) / 4 - abs(v[i][j - 1][k] + v[i][j - 1][k + 1]) * (w[i][j - 1][k] - w[i][j][k]) / 4);

                double ddw_dxdx = (w[i + 1][j][k] - 2 * w[i][j][k] + w[i - 1][j][k]) / pow(dx, 2);
                double ddw_dydy = (w[i][j + 1][k] - 2 * w[i][j][k] + w[i][j - 1][k]) / pow(dy, 2);
                double ddw_dzdz = (w[i][j][k + 1] - 2 * w[i][j][k] + w[i][j][k - 1]) / pow(dz, 2);


                H[i][j][k] = w[i][j][k] + ts * ((ddw_dxdx + ddw_dydy + ddw_dzdz) / Re - duw_dx - dvw_dy - dww_dz + gz);
            }
        }
    }
    //////////////////////////////n-3 confusing U deee***********8/////////////
}


void comp_p_rhs(double*** F, double*** G, double*** H, double*** p_rhs, double ts, double dx, double dy, double dz, int nx, int ny, int nz) {
    //cout << "rhs" << "\n";
    for (int i = 1; i <= nx - 2; i++) {
        for (int j = 1; j <= ny - 2; j++) {
            for (int k = 1; k <= nz - 2; k++) {
                p_rhs[i][j][k] = ((F[i][j][k] - F[i - 1][j][k]) / dx + (G[i][j][k] - G[i][j - 1][k]) / dy + (H[i][j][k] - H[i][j][k - 1]) / dz) / ts;
            }
        }
    }
}
/////////////////////////////cal only newly created p's////////////////////
double find_max(double*** mat, int nx, int ny, int nz) {
    double maxElement = 0.;
    for (int i = 2; i <= nx - 2; i++) {
        for (int j = 2; j <= ny - 2; j++) {
            for (int k = 2; k <= nz - 2; k++) {
                if (abs(mat[i][j][k]) > maxElement) {
                    maxElement = abs(mat[i][j][k]);
                }
            }
        }
    }

    // finally return maxElement
    return maxElement;
}
/////////////////////////////cal only newly created p's////////////////////

void comp_p(double*** p, double*** p_new, double*** p_rhs, double*** res, double dx, double dy, double dz, int nx, int ny, int nz, int it_max, double eps) {
    int it = 1;
    double r = 9999;




    /////////////////update p's periodic bc every iteration/////////////////
    while (it <= it_max && abs(r) >= eps) {
        for (int i = 1; i <= nx - 2; i++) {
            for (int j = 1; j <= ny - 2; j++) {
                p_new[i][j][0] = p[i][j][1];
                p_new[i][j][nz - 1] = p[i][j][nz - 2]; // presure strip for north and south bound
            }
        }
        for (int j = 1; j <= ny - 2; j++) {
            for (int k = 1; k <= nz - 2; k++) {
                p_new[0][j][k] = p[1][j][k];
                p_new[nx - 1][j][k] = p[nx - 2][j][k]; // pressure strip for east and west bound
                p_new[1][j][k] = p[nx - 2][j][k]; //stream_wise_periodic
            }
        }
        for (int i = 1; i <= nx - 2; i++) {
            for (int k = 1; k <= nz - 2; k++) {
                p_new[i][0][k] = p[i][1][k];
                p_new[i][ny - 1][k] = p[i][ny - 2][k];
                p_new[i][1][k] = p[i][ny - 2][k]; //span_wise_periodic
            }
        }
        /////////////////update p's periodic bc every iteration/////////////////

        ///////////////////// 2 to n-2 (1 is already updated)/////////////////////////
        double sum_rsqr = 0;
        //cout << "pressure_lang_up_kob" << "\n";
        //visualization_gen(p_new,nx,ny);
        for (int i = 2; i <= nx - 2; i++) {//////////
            for (int j = 2; j <= ny - 2; j++) {////////
                for (int k = 1; k <= nz - 2; k++) {
                    p_new[i][j][k] = (-0.7) * p[i][j][k] + (1.7 / (2 / pow(dx, 2) + 2 / pow(dy, 2) + 2 / pow(dz, 2))) * ((p[i + 1][j][k] + p_new[i - 1][j][k]) / pow(dx, 2) + (p[i][j + 1][k] + p_new[i][j - 1][k]) / pow(dy, 2) + (p[i][j][k + 1] + p_new[i][j][k - 1]) / pow(dz, 2) - p_rhs[i][j][k]);
                }
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

        /*/////////////////////////up again before cal res (not ecessary) ///////////////////////////////
        for (int i = 1; i <= nx - 2; i++) {
            for (int j = 1; j <= ny - 2; j++) {
                p_new[i][j][0] = p[i][j][1];
                p_new[i][j][nz - 1] = p[i][j][nz - 2]; // presure strip for north and south bound
            }
        }
        for (int j = 1; j <= ny - 2; j++) {
            for (int k = 1; k <= nz - 2; k++) {
                p_new[0][j][k] = p[1][j][k];
                p_new[nx - 1][j][k] = p[nx - 2][j][k]; // pressure strip for east and west bound
                p_new[1][j][k] = p[nx - 2][j][k]; //stream_wise_periodic
            }
        }
        for (int i = 1; i <= nx - 2; i++) {
            for (int k = 1; k <= nz - 2; k++) {
                p_new[i][0][k] = p[i][1][k];
                p_new[i][ny - 1][k] = p[i][ny - 2][k];
                p_new[i][1][k] = p[i][ny - 2][k]; //span_wise_periodic
            }
        }
        /////////////////////////up again before cal res///////////////////////////////*/


        for (int i = 2; i <= nx - 2; i++) {
            for (int j = 2; j <= ny - 2; j++) {
                for (int k = 2; k <= nz - 2; k++) {
                    res[i][j][k] = ((p_new[i + 1][j][k] - p_new[i][j][k]) - (p_new[i][j][k] - p_new[i - 1][j][k])) / (dx * dx)
                        + ((p_new[i][j + 1][k] - p_new[i][j][k]) - (p_new[i][j][k] - p_new[i][j - 1][k])) / (dy * dy) + (p_new[i][j][k + 1] - 2 * p_new[i][j][k] + p_new[i][j][k - 1]) / pow(dz, 2) - p_rhs[i][j][k];
                }
            }
        }

        //visualization_gen(res,nx,ny);

        //cout << "\n";
        r = find_max(res, nx, ny, nz);

        it++;
        //cout << "residual = " << r << "\n" << "interation =" << it << "\n";
        //update p
        for (int i = 0; i <= nx - 1; i++) {
            for (int j = 0; j <= ny - 1; j++) {
                for (int k = 0; k <= nz - 1; k++) {
                    p[i][j][k] = p_new[i][j][k];
                }
            }
        }
    }
    cout << "residual = " << abs(r) << "\n" << "interation =" << it << "\n";

}

/////////////////////////////// sperate cal (follow the book) ////////////////////////////////////
void finalcomp(double*** u_new, double*** v_new, double*** w_new, double*** F, double*** G, double*** H, double*** p_new, double ts, double dx, double dy, double dz, int nx, int ny, int nz)
{
    for (int i = 1; i <= nx - 3; i++) {
        for (int j = 1; j <= ny - 2; j++) {
            for (int k = 1; k <= nz - 2; k++) {
                u_new[i][j][k] = F[i][j][k] - (ts / dx) * (p_new[i + 1][j][k] - p_new[i][j][k]);
            }
        }
    }
    for (int j = 1; j <= ny - 3; j++) {//change
        for (int i = 1; i <= nx - 2; i++) {
            for (int k = 1; k <= nz - 2; k++) {
                v_new[i][j][k] = G[i][j][k] - (ts / dy) * (p_new[i][j + 1][k] - p_new[i][j][k]);
            }
        }
    }
    for (int k = 1; k <= nz - 3; k++) {//change
        for (int i = 1; i <= nx - 2; i++) {
            for (int j = 1; j <= ny - 2; j++) {
                w_new[i][j][k] = H[i][j][k] - (ts / dz) * (p_new[i][j][k + 1] - p_new[i][j][k]);
            }
        }
    }
    /////////////////////////////// sperate cal (follow the book) ////////////////////////////////////
}
void write_result(double*** gen, double dx, double dy, double dz,int nx, int ny, int nz, int t,string x) {
    //var = simulated data (array)
    //nx, ny = size (2d)
    //ts = timestep
    ofstream myfile;
    string file_name = "vtk_"+ x +"/ThreeD_" + to_string(t) + ".vtk";
    myfile.open(file_name);
    // Paraview stuff
    //Header
    myfile << "# vtk DataFile Version 2.0\n";
    myfile << "FlowField\n";
    myfile << "ASCII\n";
    // Grid
    myfile << "DATASET STRUCTURED_GRID\n";
    myfile << "DIMENSIONS " << nx  << " " << ny << " " << nz << endl;
    myfile << "POINTS " << (nx ) * ny * (nz) << " float\n";
    for (int k = 0; k <= nz - 1; k++) {
        for (int j = 0; j <= ny - 1; j++) {
            for (int i = 0; i <= nx - 1; i++) {
                myfile << i * dx << " " << j *dy << " " << k * dz << endl;
            }
        }
    }
    // Dataset
    myfile << "\n";
    myfile << "POINT_DATA";
    myfile << " " << (nx ) * ny * (nz) << "\n";
    // Point data
    myfile << "\n";
    myfile << "SCALARS "+x+" float 1\n";
    myfile << "LOOKUP_TABLE default\n";
    for (int k = 0; k <= nz - 1; k++) {
        for (int j = 0; j <= ny - 1; j++) {
            for (int i = 0; i <= nx - 1; i++) {
                myfile << gen[i][j][k] << endl;
            }
        }
    }
    // end
    myfile.close();
}



int main() {

    int nx = 11;
    int ny = 11;
    int nz = 21;
    double x_size = 4.0 * 3.1416;//streamwise
    double y_size = 2.0 * 3.1416;//spanwise
    double z_size = 2.0;//chanel thickness
    double dx = x_size / (nx - 3);
    double dy = y_size / (ny - 3);
    double dz = z_size / (nz - 2);
    double donor = 0;
    double gx = 1;
    double gy = 0;
    double gz = 0;
    int t_max = 100000;
    int it_max = 100; //for poisson iteration
    double eps = 0.01; //for poisson iteration



    //double u_in = 0;
    double Re = 50;//1500;

    //double ts = Re/((1/pow(dx,2))+(1/pow(dy,2)))/2*1.1; // choose ts based on 3.5 with a sf of 1.1
    double ts = 0.002;


    ////////////////////////////////////////////////////CREATE ARRAYS///////////////////////////////////////////

    double*** u = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        u[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            u[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }

    double*** u_new = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        u_new[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            u_new[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }

    double*** v = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        v[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            v[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }

    double*** v_new = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        v_new[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            v_new[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }

    double*** w = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        w[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            w[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }

    double*** w_new = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        w_new[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            w_new[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }

    double*** p = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        p[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            p[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }

    double*** p_new = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        p_new[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            p_new[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }

    double*** F = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        F[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            F[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }

    double*** G = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        G[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            G[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }

    double*** H = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        H[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            H[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }

    double*** p_rhs = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        p_rhs[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            p_rhs[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }
    double*** res = (double***)malloc(nx * sizeof(double**));
    for (int i = 0; i <= nx - 1; i++) {
        res[i] = (double**)malloc(ny * sizeof(double*));
        for (int j = 0; j <= ny - 1; j++) {
            res[i][j] = (double*)malloc(nz * sizeof(double));
        }
    }

    initualization(u, u_new, v, v_new, w, w_new, p, p_new, F, G, H, p_rhs, nx, ny, nz);

    setbound(u_new, v_new, w_new, p_new, p, nx, ny, nz);

    alittlestepforhumankind(u, u_new, v, v_new, w, w_new, p, p_new, nx, ny, nz);

    //showinipls ( u, v, p, nx, ny);

    // ptr=&u[0][0];



    for (int t = 0; t <= t_max - 1; t++) {
        cout << "t= " << t << "\n";
        comp_FG(u, v, w, F, G, H, dx, dy, dz, donor, Re, gx, gy, gz, ts, nx, ny, nz);

        comp_p_rhs(F, G, H, p_rhs, ts, dx, dy, dz, nx, ny, nz);
        //cout << "F\n";
        //visualization_planex(F, nx, ny,nz);

        comp_p(p, p_new, p_rhs, res, dx, dy, dz, nx, ny, nz, it_max, eps);

        finalcomp(u_new, v_new, w_new, F, G, H, p_new, ts, dx, dy, dz, nx, ny, nz);

        setbound(u_new, v_new, w_new, p_new, p, nx, ny, nz);
        // write_result(u_new,nx , ny,t);
        alittlestepforhumankind(u, u_new, v, v_new, w, w_new, p, p_new, nx, ny, nz);
        
        //visualization_planex(u_new, nx, ny, nz);
        
        //visualization_planex(p_new, nx, ny,nz);
        write_result(u_new, dx,dy,dz, nx, ny, nz, t, "u");
        write_result(v_new, dx, dy, dz, nx, ny, nz, t, "v");
        write_result(w_new, dx, dy, dz, nx, ny, nz, t, "w");
        write_result(p_new, dx, dy, dz, nx, ny, nz, t, "p");
        //write_result(u_new,nx , ny,t);
        //visualization_v(v_new,nx,ny,t);
        //visualization_p(p_new,nx,ny,t);
    }

}
