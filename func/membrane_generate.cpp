#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <tuple>
#include <sstream>
#include <iomanip>
#include <math.h>
using namespace std;



struct Lipid { //lipid structure (description)
    double x, y, z;
    string name;
    double dipole_x,dipole_y,dipole_z,dipole_total;
};

tuple<double, double, double, double> ustawDipol(double px, double dpx, double py, double dpy, double pz, double dpz, double ptot);


vector<Lipid> generateBilayer(double width, double height, double depth, double apl_A, double apl_B, double apl_C, double percent_A, double percent_B, double percent_C,
                              double dipole_A_x, double dipole_error_A_x, double dipole_A_y, double dipole_error_A_y, double dipole_A_z, double dipole_error_A_z, double dipole_A_total,
                              double dipole_B_x, double dipole_error_B_x, double dipole_B_y, double dipole_error_B_y, double dipole_B_z, double dipole_error_B_z, double dipole_B_total,
                              double dipole_C_x, double dipole_error_C_x, double dipole_C_y, double dipole_error_C_y, double dipole_C_z, double dipole_error_C_z, double dipole_C_total) {
/*
This function allows to generate the membrane consisted of appropriate dipoles. For now up tp 3 (three) lipid types are acceptable. The vector will take multiple input arguments
like membrane size (width, height, depth), area per lipid (A,B,C), membrane composition (% of lipids), and lipid dipole moments (X,Y,Z, TOT +std).
Input: double (values)....
Output: bilayer vector with X,Y,Z and tot of lipid's dipole moment /double/, lipid coordinates (X,Y,Z)/double/ and lipid type /string/.
*/
    double total_percent = percent_A + percent_B + percent_C;

    //bilayer composition handling
    if (percent_A==0) {
        apl_A, dipole_A_x,dipole_A_y,dipole_A_z,dipole_error_A_x,dipole_error_A_y,dipole_error_A_z==0;}
    else if (percent_B==0) {
        apl_B,dipole_B_x,dipole_B_y,dipole_B_z,dipole_error_B_x,dipole_error_B_y,dipole_error_B_z==0;}
    else if (percent_C==0) {
        apl_C,dipole_C_x,dipole_C_y,dipole_C_z,dipole_error_C_x,dipole_error_C_y,dipole_error_C_z==0;}

    if (total_percent != 100) {


        throw out_of_range("The lipid percentage incorrectly defined. Execution terminated.");
    }

    else {


    double total_apl=(width *height)*0.90; // Calculate the total APL and include 0.5nm spacing > is about 10%
    int num_lipids_A = static_cast<int>(total_apl/apl_A * percent_A*0.01);
    int num_lipids_B = static_cast<int>(total_apl/apl_B * percent_B*0.01);
    int num_lipids_C = static_cast<int>(total_apl/apl_C * percent_C*0.01);

    double mean_apl,total_lipids_num;
    int lipids_per_line;

    total_lipids_num=num_lipids_A+num_lipids_B+num_lipids_C;
    lipids_per_line= sqrt(total_lipids_num); // calculate the max number of lipids to put in X and Y row


    vector<Lipid> bilayer;
    bilayer.reserve(num_lipids_A + num_lipids_B + num_lipids_C); //memory reserve

    //cout<<"A :"<< num_lipids_A << " B :"<<num_lipids_B<< " C :"<<num_lipids_C<<endl;
    //cout<< "lipids_per_line: "<<lipids_per_line/2<<endl;

    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> dist(0, 1); //for extra fluctuation in z axis (defects introduction)

    vector<string> lipid_types_first;
    lipid_types_first.reserve(num_lipids_A+num_lipids_B+num_lipids_C);

    for (int i = 0; i < num_lipids_A; ++i) {
        lipid_types_first.push_back("A");
    }

    for (int i = 0; i < num_lipids_B; ++i) {
        lipid_types_first.push_back("B");
    }

    for (int i = 0; i < num_lipids_C; ++i) {
        lipid_types_first.push_back("C");
    }

    shuffle(lipid_types_first.begin(), lipid_types_first.end(), gen); //shuffling the lipids in the monolayer to avoid ordering

    vector<string> lipid_types_second = lipid_types_first;
    shuffle(lipid_types_second.begin(), lipid_types_second.end(), gen);

    int lipid_index = 0;

    // First monolayer at z = 2 +- std
    for (double y = 0; y <= lipids_per_line/2; y += 0.5) { //here the 5A spacing between dipoles is introduced
        for (double x = 0; x <= lipids_per_line/2; x += 0.5) {
            double z = abs(depth/2.0) + dist(gen) * 0.01; // Randomize z-coordinate with small deviation (membrane fluctuation)

            double dipole_x, dipole_y, dipole_z,dipole_total;
            if (lipid_types_first[lipid_index] == "A") { //generate the dipole moments (X,Y,Z, total) for appropriate type of lipid.

                tie(dipole_x, dipole_y, dipole_z, dipole_total)=ustawDipol(dipole_A_x, dipole_error_A_x, dipole_A_y, dipole_error_A_y, dipole_A_z, dipole_error_A_z, dipole_A_total);
            }
            else if (lipid_types_first[lipid_index] == "B") {

                tie(dipole_x, dipole_y, dipole_z, dipole_total)=ustawDipol(dipole_B_x, dipole_error_B_x, dipole_B_y, dipole_error_B_y, dipole_B_z, dipole_error_B_z, dipole_B_total);
            }
            else if (lipid_types_first[lipid_index] == "C") {

                tie(dipole_x, dipole_y, dipole_z, dipole_total) = ustawDipol(dipole_C_x, dipole_error_C_x, dipole_C_y, dipole_error_C_y, dipole_C_z, dipole_error_C_z, dipole_C_total);
            }

            bilayer.push_back({x, y, z, lipid_types_first[lipid_index], dipole_x, dipole_y, dipole_z, dipole_total}); //that push back all necessary lipid information for further calculations
            lipid_index++;
        }
    }

    lipid_index = 0; // Reset lipid index for second monolayer

    // Second monolayer at z = -2 +- std
    for (double y = 0; y <= lipids_per_line/2; y += 0.5) {
        for (double x = 0; x <= lipids_per_line/2; x += 0.5) {
            double z = -abs(depth/2.0)+ dist(gen) * 0.01; // Randomize z-coordinate with small deviation

            double dipole_x, dipole_y, dipole_z,dipole_total;
            if (lipid_types_second[lipid_index] == "A") {

                tie(dipole_x, dipole_y, dipole_z, dipole_total)=ustawDipol(dipole_A_x, dipole_error_A_x, dipole_A_y, dipole_error_A_y, dipole_A_z, dipole_error_A_z, dipole_A_total);
            }
            else if (lipid_types_second[lipid_index] == "B") {

                tie(dipole_x, dipole_y, dipole_z, dipole_total)=ustawDipol(dipole_B_x, dipole_error_B_x, dipole_B_y, dipole_error_B_y, dipole_B_z, dipole_error_B_z, dipole_B_total);
            }
            else if (lipid_types_second[lipid_index] == "C") {

                tie(dipole_x, dipole_y, dipole_z, dipole_total)=ustawDipol(dipole_C_x, dipole_error_C_x, dipole_C_y, dipole_error_C_y, dipole_C_z, dipole_error_C_z, dipole_C_total);
            }

            bilayer.push_back({x, y, z, lipid_types_second[lipid_index], dipole_x, dipole_y, dipole_z, dipole_total});
            lipid_index++;
        }
    }

    return bilayer;
}}


double normal( double mean,  double  std)
{
   double pii=3.1415927;
   double r_max=RAND_MAX+1;
   return std*sqrt(-2*log((rand()+1)/r_max))*sin(2*pii*rand()/r_max)+mean; //losowanie z rozkladu normalnego
}

/*double normal(double mean, double stdDev) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(mean, stdDev);
    return dist(gen);
}


 void ustawDipol(int ii, double px, double dpx, double py, double dpy, double pz, double dpz, double ptot)
{
    line1:
    d[ii][0] = normal(px, dpx);
    d[ii][1] = normal(py, dpy);
    d[ii][2] = normal(pz, dpz);
    d[ii][3] = d[ii][0] * d[ii][0] + d[ii][1] * d[ii][1] + d[ii][2] * d[ii][2];
    d[ii][3] = sqrt(d[ii][3]);

    if (abs(d[ii][3]) >= ptot)
        goto line1;

    ofstream outputFile("dipoles.txt");
    if (outputFile.is_open()) {
        for (int i = 0; i < ii; i++) {
            outputFile << d[i][0] << " " << d[i][1] << " " << d[i][2] << " " << d[i][3] << endl;
        }
        outputFile.close();
}

}
*/

tuple<double, double, double, double> ustawDipol(double px, double dpx, double py, double dpy, double pz, double dpz, double ptot) {
    /*
    That function allows to calculate the dipole moments to each lipid type. Dipoles are drawn from mean+std with normal distribution function.
    Input: Dipoles X,Y,Z + std, tot /double/
    Output: dipoles for X Y Z /double/
    */

    double dipole_x, dipole_y, dipole_z, dipole_total;
    line1:
    double qc=3.33564*1e-30;
    dipole_x = normal(px, dpx);
    dipole_y = normal(py, dpy);
    dipole_z = normal(pz, dpz);
    dipole_total = dipole_x * dipole_x + dipole_y * dipole_y + dipole_z * dipole_z;
    dipole_total = sqrt(dipole_total);

    if (abs(dipole_total) >= ptot)
        goto line1;

    return make_tuple(dipole_x, dipole_y, dipole_z, dipole_total);
}



void printLipidCoordinates(const vector<Lipid>& bilayer) {
    for (const auto& lipid : bilayer) {
        cout << "x: " << lipid.x << ", y: " << lipid.y << ", z: " << lipid.z << ", type: " << lipid.name<<endl;

    }
}

void saveLipidCoordinatesToCSV(const vector<Lipid>& bilayer, const string& filename) {
    ofstream outputFile(filename);
    if (outputFile.is_open()) {
        for (const auto& lipid : bilayer) {
            outputFile <<lipid.dipole_x<< ","  <<lipid.dipole_y<< ","  <<lipid.dipole_z<< ","  <<lipid.dipole_total<< ","  <<lipid.x << ","  << lipid.y << "," << lipid.z <<  "," << lipid.name << endl;
        }
        outputFile.close();
        cout << "Pomyslnie zapisano koordynaty lipidow do pliku " << filename << endl;
    } else {
        cout << "Blad: Nie mozna otworzyc pliku " << filename << endl;
    }
}



int main() {
    double width = 12.0; //nm
    double height = 12.0;
    double depth = 4.0;
    double apl_A = 0.65;
    double apl_B = 0.72;
    double apl_C = 0.72;
    double percent_A = 50.0; // Procentowy udzia³ lipida A
    double percent_B = 50.0; // Procentowy udzia³ lipida B
    double percent_C = 0.0; // Procentowy udzia³ lipida C

    double dipole_A_x = 0.14;
    double dipole_error_A_x = 10.17;
    double dipole_A_y = 0.29;
    double dipole_error_A_y = 10.09;
    double dipole_A_z = -35.08;
    double dipole_error_A_z = 8.29;
    double dipole_A_total=39.19;

    double dipole_B_x = 0.34;
    double dipole_error_B_x = 11.29;
    double dipole_B_y = -0.37;
    double dipole_error_B_y = 11.39;
    double dipole_B_z = 1.65;
    double dipole_error_B_z = 8.44;
    double dipole_B_total=18.27;

    //double dipole_B_x,dipole_B_y,dipole_B_z,dipole_error_B_x,dipole_error_B_y,dipole_error_B_z;
    double dipole_C_x,dipole_C_y,dipole_C_z,dipole_error_C_x,dipole_error_C_y,dipole_error_C_z;

    vector<Lipid> bilayer = generateBilayer(width, height, depth, apl_A, apl_B, apl_C, percent_A, percent_B, percent_C,
                                            dipole_A_x,dipole_error_A_x,dipole_A_y,dipole_error_A_y,dipole_A_z,dipole_error_A_z,dipole_A_total,
                                            dipole_B_x,dipole_error_B_x,dipole_B_y,dipole_error_B_y,dipole_B_z,dipole_error_B_z,dipole_B_total,
                                            dipole_C_x,dipole_error_C_x,dipole_C_y,dipole_error_C_y,dipole_C_z,dipole_error_C_z,50.0);

    printLipidCoordinates(bilayer);

    string filename = "lipid_coordinates.csv";
    saveLipidCoordinatesToCSV(bilayer, filename);



    return 0;
}
