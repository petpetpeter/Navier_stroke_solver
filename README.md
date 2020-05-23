# Navier_stroke_solver
2D-Channel Flow inlet boundary condition and 3D-Channel Flow periodic boundary condition
## How to run a code
1. Clone an respiratory
```bash
git clone https://github.com/petpetpeter/Navier_stroke_solver.git
```
2. Create output folder in respiratory folder (name format : "variable" + "_vtk")
```bash
cd Navier_stroke_solver
mkdir vtk_u
mkdir vtk_v
mkdir vtk_w
mkdir vtk_p
```
3. Edit parameter in the main function
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
5. Run your_program_name and View the result in Paraview
- It may take 2 hours to finish 100,000 time step for 3D-simulation.
- For faster simulation with lower detail, you may change the parameter in main as follow
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
