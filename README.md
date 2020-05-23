# Navier_stroke_solver
2D-Channel Flow inlet boundary condition and 3D-Channel body-force driven flow with periodic boundary condition
## How to run the code
1. Clone the provided respiratory
```bash
git clone https://github.com/petpetpeter/Navier_stroke_solver.git
```
2. Create output folders "in the respiratory folder" (name format : "variable" + "_vtk")
```bash
cd Navier_stroke_solver
mkdir vtk_u
mkdir vtk_v
mkdir vtk_w
mkdir vtk_p
```
3. Edit parameter in the main function
All parameters should be set approprieately to ensure stability.  Recommended proven combinations are provided as the defualt combinations for both the 2D and 3D files.
```c++
int main() {
    int nx = 11; //number_of_grid
    int ny = 11; 
    int nz = 21; 
    double x_size = 4.0 * 3.1416;//streamwise_length
    double y_size = 2.0 * 3.1416;//spanwise_length
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
    double Re = 50;
    double ts = 0.002;
    //In case of Inlet Boundary Condition
    double u_in = 1.0;
```
4. Compile channel_flow_*.cpp
```bash
g++ -o your_program_name channel_flow_2d.cpp
g++ -o your_program_name channel_flow_3d.cpp
```
5. Run your_program_name. The output files will apear as .vtk files (which can be viewed in Paraview) in the folders created in step 2.
- It may take more than 2 hours to finish 100,000 time steps for the 3D-simulation.
- For a faster 3D simulation, you may change the parameters. A proven combination that provides faster simulation (aprox. 15 min) but lower detail is shown below as an example. 
```c++
int main() {
    int nx = 7; //number_of_grid
    int ny = 7; 
    int nz = 9; 
    double x_size = 4.0 * 3.1416;//streamwise_length
    double y_size = 2.0 * 3.1416;//spanwise_length
    double z_size = 2.0;//chanel thickness
    double dx = x_size / (nx - 3);
    double dy = y_size / (ny - 3);
    double dz = z_size / (nz - 2);
    double donor = 0;
    double gx = 1;
    double gy = 0;
    double gz = 0;
    int t_max = 30000;
    int it_max = 100; //for poisson iteration
    double eps = 0.01; //for poisson iteration
    double Re = 50;
    double ts = 0.005;
    //In case of Inlet Boundary Condition
    double u_in = 1.0;
```
6. (OPTIONAL) Unwanted outputs can be deselected by commenting out "writeresult" functions associated with the unwanted outputs at the bottom of the .cpp files ( all outputs:u,v,w,and p are selected by defualt).
For an example, if p values are unwanted as an ouput:
```c++
 write_result(u_new, dx,dy,dz, nx, ny, nz, t, "u");
 write_result(v_new, dx, dy, dz, nx, ny, nz, t, "v");
 write_result(w_new, dx, dy, dz, nx, ny, nz, t, "w");
// write_result(p_new, dx, dy, dz, nx, ny, nz, t, "p");
