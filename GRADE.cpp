//
//                  ***GRADE_v1.00***
//  MyFunctions.hpp
//
//  Created by Farbod Mahmoudinobar and Cristiano L. Dias on 12/6/18.
//  Copyright © 2018 Farbod Mahmoudinobar and Cristiano L. Dias.
//
//  This file is part of GRADE.
//
//  GRADE is a free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  GRADE is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with GRADE.  If not, see <https://www.gnu.org/licenses/>.

#include <bits/stdc++.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <unistd.h>
#include<chrono>
#include "omp.h"
#include "MyFunctions.hpp"

using namespace std;

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::vector<std::pair<T1, T2>>& vec) {
    os << "[ ";
    for (const auto& pair : vec) {
        os << "(" << pair.first << ", " << pair.second << ") ";
    }
    os << "]";
    return os;
}

//overload to be able to print out maps
template<typename key, typename value>
ostream& operator<<(ostream& os, const std::map<key,value>& M)
{
    os << "{";
    for(auto pair : M)
        os << pair.first << ":" << pair.second << ",  ";
    os << "}";
    return os;

}

//overload ostream to be able to cout vectors
template <typename T>
ostream& operator<<( ostream& os, const vector<T>& vec)
{
    for ( auto element : vec)
    {os << element<< "  ";}
    return os;
}    

bool contains(const std::vector<std::string>& container, std::string element)
{
    auto it = std::find(container.begin(), container.end(), element);
    return it != container.end();
}

//template concat function
template <typename T> 
std::vector<T> concat(std::vector<T>& a , std::vector<T>& b)
{
    std::vector<T> c;
  for(  auto e: a){
      c.push_back(e);
    }
    for( auto e: b)
    {
        c.push_back( e);
    }
    return c;
}


// Main function ----------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    double const Version = 1.00;

    int cage_512_count;
    int cage_62512_count;
    int cage_64512_count;
    //test contains()method
    //vector<int> ar2 = {9,8,10};
    //vector<int> ar3 = concat(ar1, ar2);
    //cout << "concat test: " << ar3 << endl;
    //if( contains(arr, "vegana"))
            //cout << "Haan bhai hai" << endl;
    //else cout << "Nahi hai" << endl;

    // output program information.
    cout << std::fixed;
    cout << std::setprecision(2);
    cout << "\n\n\t\t** GRADE - VERSION " << Version << " **\n\n";

    string inputFilename, s1, s2;
    string outputFilename, rawFilename, outputFilename2;
    int DT = 1, FR = 1, THETA = 45;
    int in_theta = 0, in_fr = 0, in_dt = 0, in_r = 0, in_F4 = 0, in_d1 = 0, in_d2 = 0, in_s1 = 0, in_s2 = 0; // This parameter is for whether theta is given as input parameter (1) or taken as default value (0).
    double HBOND_DIST = 0.35;                                                                                // Command line input parameter for Hbond_distance cutoff, taken from '-r' flag.
    double delta_p = 0.18, delta_h = 0.26;                                                                   // Command line input parameter for delta constraints for pentagon and hexagons.
    string F4 = "no";

    // Parsing command line parameters (input taken from "-i" and output taken from "-o")
    std::string arg;
    std::string inName, outName;

    //-------------------------------------------------------------------------------------

    int argi = 1;

    while (argi < argc)
    {
        const char *args = argv[argi++];

        switch (*args)
        {
        case '-':
            if (strcmp(args, "-i") == 0)
            {
                inputFilename = argv[argi];
            }
            else if (strcmp(args, "-o") == 0)
            {
                outputFilename = argv[argi];
            }
            else if (strcmp(args, "-dt") == 0)
            {
                DT = atoi(argv[argi]);
                in_dt = 1;
            }
            else if (strcmp(args, "-fr") == 0)
            {
                FR = atoi(argv[argi]);
                in_fr = 1;
            }
            else if (strcmp(args, "-theta") == 0)
            {
                THETA = atoi(argv[argi]);
                in_theta = 1;
            }
            else if (strcmp(args, "-r") == 0)
            {
                HBOND_DIST = atof(argv[argi]);
                in_r = 1;
            }
            else if (strcmp(args, "-f4") == 0 || strcmp(args, "-F4") == 0)
            {
                in_F4 = 1;
                F4 = "YES";
                if (argi < argc)
                    F4 = argv[argi];
                if (F4 == "no" || F4 == "No" || F4 == "NO")
                {
                    in_F4 = 0;
                    F4 = "NO";
                }
            }
            else if (strcmp(args, "-d1") == 0)
            {
                delta_p = atof(argv[argi]);
                in_d1 = 1;
            }
            else if (strcmp(args, "-d2") == 0)
            {
                delta_h = atof(argv[argi]);
                in_d2 = 1;
            }

            else if (strcmp(args, "-s1") == 0)
            {
                s1 = argv[argi];
                in_s1 = 1;
            }
            else if (strcmp(args, "-s2") == 0)
            {
                s2 = argv[argi];
                in_s2 = 1;
            }

            else
                cout << "Skipped unknown option(s): '" << args << " " << argv[argi] << "'"
                     << "\n\n";
            break;
        }
    }

    if (inputFilename.empty())
    {
        cout << "**No Input Provided!**"
             << "\n\n";
        print_usage();
        cout << "**No Input Provided!**"
             << "\n\n";
        return 0;
    }
    else
    {
        print_usage();
    }

    if (outputFilename.empty())
    {
        outputFilename = inputFilename;
    }

    cout << "Command line:\n";
    cout << argv[0] << " -i " << inputFilename << " -o " << outputFilename;
    if (in_dt == 1)
    {
        cout << " -dt " << DT;
    }
    if (in_fr == 1)
    {
        cout << " -fr " << FR;
    }
    if (in_theta == 1)
    {
        cout << " -theta " << THETA;
    }
    if (in_r == 1)
    {
        cout << " -r " << HBOND_DIST;
    }
    if (in_F4 == 1)
    {
        cout << " -f4 "
             << "YES";
    }
    if (in_d1 == 1)
    {
        cout << " -d1 " << delta_p;
    }
    if (in_d2 == 1)
    {
        cout << " -d2 " << delta_h;
    }
    if (in_s1 == 1)
    {
        cout << " -s1 " << s1;
    }
    if (in_s2 == 1)
    {
        cout << " -s1 " << s2;
    }

    //cout << "\n\n";
    //-------------------------------------------------------------------------------------

    vector<vector<int>> My_neigh;
    int Natoms, count_solute = 0, count_solute2 = 0;
    vector<int> Nneigh;
    vector<vector<double>> atom_positions;
    vector<vector<double>> solvent_positions;
    vector<vector<double>> solute_positions;
    double F4_value = 0; // Value of F4 order parameter for each frame.
    vector<double> current_F4_line;
    vector<vector<double>> time_vs_F4; // This variable holds Time/frame_counter in first column and value of F4 order parameter in second column.
    vector<std::string> molstrings;

    string line = "NONE", str1, str2, str3;
    int int1;
    double x, y, z;
    double boxX = 0.0, boxY = 0.0, boxZ = 0.0;
    vector<double> temp_vect;
    //vector<int> temp_vect2;
    int lineNumber = 0;
    int firstSOL = 0;
    int frameCounter = 0;
    int topSolute = 0;
    size_t methane_512 = 0, methane_62512 = 0, methane_64512 = 0;
    string time;
    Natoms = 0;
    string solute1_norm = "AAA";
    string solute2_norm = "BBB";
    string solute1,solute2,solute3,solute4;
    map<string,int>mp,map_count,real_map;

    ofstream outFile;
    ofstream outFile2;
    size_t found_ext = outputFilename.find_last_of(".");
    rawFilename = outputFilename.substr(0, found_ext);
    outputFilename = rawFilename + ".xvg";
    outputFilename2 = rawFilename + "_cageF4_info.dat";
    string temp1, temp2, solute1_atom_finder, solute2_atom_finder;
    temp1 = rawFilename + "_cage62512.gro";
    temp2 = rawFilename + "_cage512.gro";
    remove(temp1.c_str());
    remove(temp2.c_str());

    // Create the header for outputfile.
    outFile.open(outputFilename, ofstream::app);
    outFile2.open(outputFilename2, ofstream::app);
    outFile << "# ------------------------------------------------------------------------------------------------------- " << endl;
    outFile << "#|(Frame) Time(ps)\t\t|cage\t|filled_cage\t|cage\t|filled_cage\t|cage\t|filled_cage\t|" << endl;
    outFile << "#|\t\t\t\t|5¹²\t|5¹²\t\t|6²5¹²\t|6²5¹²\t\t|6⁴5¹²\t|6⁴5¹²\t\t|" << endl;
    outFile << "# ------------------------------------------------------------------------------------------------------- " << endl;
    outFile2 << "Frame " << std::string(10, ' ') << "Time(ps)" << std::string(35-14, ' ') << "5¹²cage" << " 5¹²cage-filled" << " 6²5¹²cage" << " 6²5¹²cage-filled" << " 6⁴5¹²cage" << " 6⁴5¹²cage-filled " << "F4_value " << std::endl;

    ifstream fileIN;
    fileIN.open(inputFilename);
    //permute( inputFilename);
    //fileIN.open(permutedFile);

    //ofstream log("oxygen_atoms", ios::out);
    // Error Check
    //if (fileIN.fail())
    //{
        //cerr << "Error Reading File"
             //<< "\n";
        //exit(1);
    //}

    ofstream outdata;
    outdata.open(rawFilename+"_f4_coords.xyz");
    outdata.close();
    //This file is not needed
    //ofstream outFile_F4;
    //if(in_F4 == 1)                          //If F4 flag option is on, open a file for F4 as a function of time.
    //{
    //    remove("F4.xvg");       //Remove any existing F4.xvg file and create a new one. 
    //    outFile_F4.open("F4.xvg", ofstream::app);
    //    outFile_F4 << "# --------------------------------------- \n" ;
    //    outFile_F4 << "#|Frame\t|F4\t\t|Time(ps)\t|" << endl;
    //    outFile_F4 << "# --------------------------------------- \n" ;

    //}

    // ... //

    // Start reading the input file.
    int iter = 0;

    while (!fileIN.eof())
    {
        cout << "\nStarted reading input file" << endl;
        getline(fileIN, line);
        lineNumber++;
        // Following if statement makes sure the "Natoms" always has the correct number, even if the code skips the first frame due to FR being greater than 1. (Added in v1.18)
        if (lineNumber == 2 || (lineNumber == (2 + 3 * frameCounter + frameCounter * Natoms) && !fileIN.eof()))
        {
            cout << "Reading 2nd line of .gro file for the no.of atoms" << endl;

            istringstream streamA(line);

            streamA >> Natoms;

            cout << " Number of atoms = " << Natoms << endl;
        }
        size_t found = 0;
        size_t found_time = 0;
        if (line.find("t="))
        {
            found_time = line.find("t=");
            //cout << "Found time " << std::endl;
        }
        if (lineNumber == 1 || (lineNumber == (1 + 3 * frameCounter + frameCounter * Natoms) && !fileIN.eof())) // This is to find the first line of gro file.
        {
            frameCounter++;

            if (((frameCounter - 1) % FR) != 0)
                continue; // If frameCounter-1 is not a multiple of FR, continue to next frame.
                          //"-1" is to ignore the frame with t=0.000. This ensures that FR=20 reads frames with t=0,20,40,... .
            if (found_time != string::npos) // If found_time is not null, do the following.
            {

                time = line.substr(found_time + 3);
                outFile << "(" << frameCounter << "\t) " << time << "\t\t| ";
                //cout <<"NO core dumped here " << std::endl;
                //if(to_string(frameCounter).length() < 5 && time.length()< 16)
                outFile2 << frameCounter << std::string(5 - to_string(frameCounter).length(),' ') << time << std::string(35 - time.length(), ' ');
                //else outFile2 << frameCounter << " " << time << " " ;
                //cout << "line 393" << std::endl;
            }
            else
            {   outFile << frameCounter << "\t\t\t| ";
                outFile2 << frameCounter << std::string(4, ' ') << std::string(16, ' ');
            }

            cout << " frame#: " << frameCounter << ", ";
            if (found_time != string::npos)
                cout << line.substr(found_time) << " ps\n";
            else
                //cout << "\n";

            if (found == string::npos)
                cout << line << endl;

            getline(fileIN, line); // Read 2nd line (number of atoms)
            lineNumber++;

            istringstream streamA(line);

            streamA >> Natoms;

            //vector<int> temp_vec = {0, 0, 0};

            int count_solvent = 0;
            vector<pair<std::string,vector<double>>> solute_atoms;
            vector<Molecule> molecules;
            vector<std::string> arr_str2;
            vector<std::string> solute_names;
            count_solute = 0;
            temp_vect = {0, 0, 0};
            atom_positions.clear();
            molstrings.clear();
            //molecules.clear();
            atom_positions.push_back(temp_vect);

            //cout << "line 396" << endl;
            int it_count = 0;
            //while (atom_positions.size() <= Natoms && atom_positions.size() < 9999)
            //cout << "Natoms = " << Natoms << std::endl;
            while (atom_positions.size() <= Natoms)
            {

                temp_vect.clear();

                getline(fileIN, line);
                lineNumber++;

                istringstream streamA(line);

                //cout << "No core dumped here " << std::endl;
                if(atom_positions.size() <= 9999)
                {
                     streamA >> str1 >> str2 >> int1 >> x >> y >> z;

                     Point coordinates(x,y,z);
                     Atom atom1(str2, coordinates);
                     if(!contains(molstrings, str1)){
                         molstrings.push_back(str1);
                         Molecule mol1(str1);
                         mol1.append_atom(atom1);
                         molecules.push_back(mol1);
                     }
                     else{
                         Molecule mol1;
                         mol1 = molecules.back();
                         mol1.append_atom(atom1);
                         molecules.back() = mol1;
                     } 
                             //cout << "str1 = " << str1 << " str2 = " << str2 << " x = " << x << " ";
                }
                if( atom_positions.size() > 9999)
                {
                     streamA >> str1 >> str2 >> x >> y >> z;

                     Point coordinates(x,y,z);
                     Atom atom1(str2, coordinates);
                     if(!contains(molstrings, str1)){
                         molstrings.push_back(str1);
                         Molecule mol1(str1);
                         mol1.append_atom(atom1);
                         molecules.push_back(mol1);
                     }
                     else{
                         Molecule mol1;
                         mol1 = molecules.back();
                         mol1.append_atom(atom1);
                         molecules.back() = mol1;
                     } 
                }

                temp_vect.push_back(x);
                temp_vect.push_back(y);
                temp_vect.push_back(z);

                atom_positions.push_back(temp_vect);

                if (line.find("wat") != string::npos || line.find( "SOL") != string::npos) // If line includes "SOL", then add to number of count_solvent.          Else, add to number of count_solute.
                {
                    count_solvent++;
                    //cout << " count_solvent = " << count_solvent << endl;
                    //cout << " lineNumber = " << lineNumber << endl;
                    solvent_positions.push_back(temp_vect);
                    //cout << "atom_positions " << atom_positions[atom_positions.size() -1]<< endl;
                    //cout << "Solvent position = " << temp_vect << endl;
                    
                    Molecule mol1;
                    mol1 = molecules.back();
                    mol1.is_solute = false;
                    molecules.back() = mol1;
                    if (count_solvent == 1)
                    {
                        firstSOL = lineNumber;
                    }
                    if (line.find("OW") != string::npos)
                    {
                        //log << line;
                        //log <<"\n";
                    }
                }
                else
                {
                    // if (count_solute == 0)
                    // {
                    //     solute1_norm = line.substr(5, 7); // Get the name of first solute in system.
                    //     size_t found_space = solute1_norm.find_first_of(" ");
                    //     solute1_norm = solute1_norm.substr(0, found_space);
                    // }

                    // if (line.substr(5, 7) != solute1_norm)
                    // {

                    //     solute2 = line.substr(5, 7); // Get the name of second solute in system.
                    //     size_t found_space = solute2.find_first_of(" ");
                    //     solute2 = solute2.substr(0, found_space);

                    //     count_solute2++;
                    // }
                    arr_str2.push_back(str2);
                    solute_positions.push_back(temp_vect);
                    solute1_norm = line.substr(5, 7); // Get the name of first solute in system.
                    //cout << " solute1_norm = " << solute1_norm << endl;
                    
                    pair<string, vector<double>> p1;
                    //cout << " lineNumber " << lineNumber << " has solute atom!" << endl; // to comment out later
                    //cout << line << endl;
                    size_t found_space = solute1_norm.find_first_of(" ");
                    solute1_norm = solute1_norm.substr(0, found_space);

                    p1.first = solute1_norm;
                    p1.second = temp_vect;
                    if(!contains( solute_names, solute1_norm)){
                        //cout << " New solute encountered: " << solute1_norm << endl;
                        //cout << " already encountered this solute name " << solute1_norm << endl;
                        solute_names.push_back(solute1_norm);
                        //cout << " lineNumber = " << lineNumber << endl;
                    }
                    else {
                    //    //cout << " New solute encountered: " << solute1_norm << endl;
                          //cout << " already encountered this solute name " << solute1_norm << endl;
                    //    //cout << " lineNumber = " << lineNumber << endl;
                    //    //cout << solute1_norm << " added to container" << endl;
                    }
                    //cout << "solute1_norm "<< solute1_norm << endl;
                    if(mp.find(solute1_norm)==mp.end()){
                       mp[solute1_norm]=count_solute+count_solvent;
                    }
                    map_count[solute1_norm]++;
                    //cout << "map_count[solute1_norm] = " << map_count[solute1_norm] <<endl;
                    count_solute++;
                    solute_atoms.push_back( p1);
                    //molecules.push_back(mol);
                }
                //cout << " count_solvent = " << count_solvent << endl;
	
		//cout << "count_solute = " << count_solute << std::endl;
                it_count++;
            }
            //atom_positions.clear();
            //cout << "atom_positions before concat = " << atom_positions.size( ) << endl;
            atom_positions.erase( atom_positions.begin( )+1,atom_positions.end()) ;

            auto v1 = concat(atom_positions, solute_positions);
            //cout << " v1 = " << v1.size( ) << endl;
            atom_positions = concat( v1, solvent_positions);
            //atom_positions = concat(solute_positions, solvent_positions);
            v1.clear();
            solvent_positions.clear();
            //cout << "atom_positions " << atom_positions.size() << endl;

            while (atom_positions.size() <= Natoms && atom_positions.size() > 9999)
            {
                //This loop is never executed

                temp_vect.clear();
                cout << " THIS LOOP IS BEING EXECTU" << endl;

                getline(fileIN, line);
                lineNumber++;
                cout << " Executing the next while loop " << endl;

                istringstream streamA(line);

                streamA >> str1 >> str2 >> x >> y >> z;

                temp_vect.push_back(x);
                temp_vect.push_back(y);
                temp_vect.push_back(z);

                //atom_positions.push_back(temp_vect);

                if (line.find("wat") != string::npos) // If line includes "SOL"
                {
                    count_solvent++;
                    //cout << "count_solvent = " << count_solvent << endl;
                }
                else
                {
                    

                    // if (line.substr(5, 7) != solute1)
                    // {
                    //     // solute2 = line.substr(5, 7);
                    //     // size_t found_space = solute2.find_first_of(" ");
                    //     // solute2 = solute2.substr(0, found_space);
                    //     //count_solute2++;
                    // }
                        solute2_norm = line.substr(5, 7);
                        size_t found_space = solute2_norm.find_first_of(" ");
                        solute2_norm = solute2_norm.substr(0, found_space);
                        if(mp.find(solute2_norm)==mp.end()){
                            mp[solute2_norm]=count_solute+count_solvent;
                        }
                        map_count[solute2_norm]++;
                        count_solute++;
                        pair<string, vector<double>> p2;
                        p2.first = solute2_norm;
                        p2.second = temp_vect;

                        solute_atoms.push_back(p2);
                }
                cout << "count_solvent = " << count_solvent << endl;
            }

            getline(fileIN, line); // get the box size from last line of frame
            // End of reading each frame.

            string box_size_xyz;
            box_size_xyz = line;
            lineNumber++;

            istringstream streamB(line);
            streamB >> boxX >> boxY >> boxZ;

            if (frameCounter == 1)
            {
                cout << "frameCounter = 1" << endl;
                topSolute = firstSOL - 3;
                //cout << " firstSOL = " << firstSOL << endl;
                //cout << "topSolute = " << topSolute << endl;
                if (in_s1 == 0)
                {
                    cout << "topSolute = " << topSolute << endl;
                    //cout << "solute1: " << solute1 << " " << topSolute << " atoms ";
                    //cout <<  "solute1: " << topSolute << " atoms"; 
                    cout << "count of solvent atoms = " << count_solvent << endl;
                    cout << " count of solute atoms = " << count_solute << endl;
                    int assign = 0;
                    cout << " solute_names = " << solute_names << endl;
                    for(auto solute : solute_names){
                        real_map[solute] = assign;
                        cout << " count of " << solute << " atoms = " << map_count[solute] << endl;
                        assign += map_count[solute];
                        //real_map[solute] = map_count[solute];
                    }
                    cout << " 1. real_map = " << real_map << endl;
                    //cout << "solute1: " << map_count[solute1_norm] << endl;
                }
                else
                {
                    //cout << "topSolute = " << topSolute;
                    cout << "solute1: " << s1 << " " << topSolute << " atoms ";
                }
                // if( (in_s1==0 && in_s2==1) || (in_s1==1 && in_s2==1) || (in_s1==0 && in_s2==0))
                {

                    if (strcmp(solute1.c_str(), solute2.c_str()) != 0)
                    {
                        if (in_s2 == 0)
                            cout << ", solute2: " << solute2 << " " << count_solute2 - topSolute << " atoms\n";
                        else
                            cout << ", solute2: " << s2 << " " << count_solute2 - topSolute << " atoms\n";
                    }
                    else
                        cout << "\n";
                }
            }
          //real_map=mp;
          //int val=INT_MAX;
          //string now;
          //for(auto i:mp){
                //if(val>i.second){
                    //val=i.second;
                    //now=i.first;
                //}
            //}
            //solute1=now;
            //mp.erase(now);
            //val=INT_MAX;
            //if(mp.size()>0){
                //for(auto i:mp){
                    //if(val>i.second){
                        //val=i.second;
                        //now=i.first;
                //}
                //}
                //solute2=now;
                //mp.erase(now);
                //if(mp.size()>0){
                    //val=INT_MAX;
                    //for(auto i:mp){
                    //if(val<i.second){
                    //val=i.second;
                    //now=i.first;
                //}
                //}
                //solute3=now;
                //mp.erase(now);
                //if(mp.size()>0){
                    //val=INT_MAX;
                    //for(auto i:mp){
                    //if(val<i.second){
                    //val=i.second;
                    //now=i.first;
                //}
                //}
                //solute4=now;
                //mp.erase(now);
                //}
            //}
            //}
            
            // Start of calculations for each frame

            cout << "count_solute = " << count_solute << endl;
            topSolute =  count_solute;
            cout << "count_solvent = "<< count_solvent << endl;
            //cout << " running calc_Distance( )" << endl;
            //cout <<"My_neigh = " << My_neigh<<endl;
            calc_Distance(count_solvent, count_solute, My_neigh, atom_positions, boxX, boxY, boxZ, Nneigh, Natoms, time, HBOND_DIST);
            //cout << " Nneigh = " << Nneigh.size()<< endl;

            vector<vector<int>> ring5, ring6, ring5_temp, ring6_temp;
            vector<vector<int>> My_neigh_ring6, My_neigh_ring5, My_neigh_ring6_ring5;

            //cout << " Nneigh = " << Nneigh << endl;
            //auto start = std::chrono::high_resolution_clock.now( );
            //cout << " My_neigh after calc_Distance()= " << My_neigh.size() << endl;
            ring_Finder(count_solute, Natoms, Nneigh, My_neigh, ring5_temp, ring6_temp, topSolute, count_solvent, atom_positions, boxX, boxY, boxZ, HBOND_DIST, delta_p, delta_h);
            //cout << "ring5_tempafter running ring_Finder = " << ring5_temp.size() << endl;
            //auto stop = std::chrono::_V2::high_resolution_clock.now();
            //auto duration = std::chrono::duration_cast<std::chrono::microseconds>( start - stop);
            //cout << "ring_Finder runtime: " << duration.count() << endl;
            //cout << " ring_Finder ran succesfully" << endl;

            if (ring5_temp.size() > 0)
            {
                coplanar_Points(ring5_temp, atom_positions, time, ring5, THETA); // Find the 5-rings which form a plane and get rid of the rest.
            }

            //cout << "removing duplicate 5-rings" << endl;
            //cout << "ring5 before removing duplicates " << ring5.size() << endl;
            int count_ring5 = remove_duplicates_map_rings(ring5); // Remove duplicate lines from 5-rings.

            cout << "ring[5]: " << count_ring5 << "\n";

            if (ring6_temp.size() > 0)
                coplanar_Points(ring6_temp, atom_positions, time, ring6, THETA);

            int count_ring6 = remove_duplicates_map_rings(ring6);

            cout << "ring[6]: " << count_ring6 << "\n";

            vector<unsigned long int> N_ring5_neigh; // Vector of ring-neighbours of all rings.

            //cout << " Running find_shared_edges_ring5 " << endl; // to comment out later
            find_shared_edges_ring5(count_ring5, ring5, My_neigh_ring5, N_ring5_neigh);

            vector<unsigned long int> N_ring6_neigh;

            vector<unsigned long int> N_ring6_ring5_neigh;

            //cout << " Running find_shared_edges_ring6_ring5 " << endl; // to comment out later
            find_shared_edges_ring6_ring5(count_ring6, count_ring5, ring5, ring6, My_neigh_ring6_ring5, N_ring6_ring5_neigh);

            if (count_ring5 == 0 && count_ring6 == 0)
            {
                cout << "\n**NO RINGS FOUND!**\n\n";
            }

            vector<vector<int>> cup512;
            //cout << " Running cup_512_Finder; " << endl; // comment out later
            cup_512_Finder(ring5, count_ring5, N_ring5_neigh, My_neigh_ring5, cup512);
            int count_512_cups = 0;
            if (cup512.size() != 0)
            {
                count_512_cups = remove_duplicates_map(cup512);
            } // Remove duplicate lines from 512 cups.

            if (count_512_cups == 0)
            {
                cout << "\n**NO 5⁶ CUPS FOUND**\n\n";
            }
            else
            {
                cout << "# 5⁶ \t cup: " << count_512_cups << "\n";
            }

            vector<vector<int>> cup62512;
            //cout << " Running cup_62512_Finder; " << endl;
            cup_62512_Finder(ring6, count_ring6, My_neigh_ring6_ring5, N_ring6_ring5_neigh, ring5, N_ring5_neigh, My_neigh_ring5, cup62512);
            int count_62512_cups = 0;
            if (cup62512.size() != 0)
            {
                count_62512_cups = remove_duplicates_map(cup62512);
            }

            if (count_62512_cups == 0)
            {
                cout << "\n**NO 6¹5⁶ CUPS FOUND**\n\n";
            }
            else
            {
                cout << "# 6¹5⁶ \t cup: " << count_62512_cups << "\n";
            }

            vector<vector<int>> cage_512, cage_62512, cage_64512;
            vector<vector<int>> cage_512_rings, cage_62512_rings, cage_64512_rings;

            cage_512_count = 0;

            if (count_512_cups > 0)
            {

                //cout << "count_512_cups = " << count_512_cups << endl;
                cage_512_count = cage_Finder(cup512, count_512_cups, My_neigh_ring5, cage_512, cage_512_rings, time);

                //cout << "cage_512_count = " << cage_512_count << endl;
                cage_512_count = remove_duplicates_map(cage_512_rings);

                cout << "real_map = " << real_map << endl;
                //cout << "No core dumped here either" << endl;
                //cout << "solute1 = " << solute1 << endl;
                //cout << "solute2 = " << solute2 << endl;
                //cout << "solute3 = " << solute3 << endl;
                //cout << "methane_512 = " << methane_512 << endl;
                //solute1 = "CO2";
                //solute2 = "METH";
                //solute3 = "Met";
                // Write output gro file only if frameCounter is a multiple of DT. if -dt is not provided, every frame will be written.
                if (cage_512_count > 0 && (frameCounter % DT == 0))
                    print_vmd_cage64512_frings(cup512, cage_512_count, cage_512_rings, ring5, ring6, atom_positions, time, rawFilename, box_size_xyz, solute_positions, methane_512, topSolute, solute1, solute2,solute3, solute4, frameCounter,map_count,real_map,arr_str2, solute_atoms, molecules);
            }

            cout << " Ab dekhte hai" << endl;
            cout << "# 5¹²\tcage: " << cage_512_count << "\n";

            cage_62512_count = 0;
            if (count_62512_cups > 0)
            {
                cage_62512_count = cage_Finder(cup62512, count_62512_cups, My_neigh_ring5, cage_62512, cage_62512_rings, time);

                cage_62512_count = remove_duplicates_map(cage_62512_rings);

                // Write output gro file only if frameCounter is a multiple of DT. if -dt is not provided, every frame will be written.
                if (cage_62512_count > 0 && (frameCounter % DT) == 0)
                    print_vmd_cage64512_frings(cup62512, cage_62512_count, cage_62512_rings, ring5, ring6, atom_positions, time, rawFilename, box_size_xyz, solute_positions, methane_62512, topSolute, solute1, solute2,solute3, solute4, frameCounter,map_count,real_map,arr_str2, solute_atoms, molecules);
            }
            //cout << "No core dumped until here " << endl;
            cout << "# 6²5¹²\tcage: " << cage_62512_count << "\n";
            /**********************************************/
            /*Finding Cage 64512*/
            cage_64512_count = 0;

            cage_64512_count = cage_Finder_64512(cup62512, count_62512_cups, cage_64512_rings);

            //cout << " Now showing the output of method cage_Finder_64512() " << endl; // comment out later
            //cout << " solute_atoms array = " << solute_atoms << endl;
            cout << "# 6⁴5¹²\tcage: " << cage_64512_count << "\n\n";

            //cout << "arr_str2 = " << arr_str2 << endl;
            print_vmd_cage64512_frings(cup62512, cage_64512_count, cage_64512_rings, ring5, ring6, atom_positions, time, rawFilename, box_size_xyz, solute_positions, methane_64512, topSolute, solute1, solute2,solute3, solute4, frameCounter,map_count,real_map, arr_str2, solute_atoms, molecules);

            /**********************************************/

            outFile << cage_512_count << "\t| " << methane_512 << "\t\t| " << cage_62512_count << "\t| " << methane_62512 << "\t\t| " << cage_64512_count << "\t| " << methane_64512 << endl;
            //cout << "methane_512 = " << methane_512 << std::endl;   
            outFile2 << std::string(5,' ')<< cage_512_count << std::string(8-to_string(cage_512_count).length(), ' ') << methane_512 << std::string(15-to_string(methane_512).length(), ' ') << cage_62512_count << std::string(10-to_string(cage_62512_count).length(), ' ') << methane_62512 << std::string(18-to_string(methane_62512).length(), ' ') << cage_64512_count << std::string(10-to_string(cage_64512_count).length(), ' ') << methane_64512 << std::string(18-to_string( methane_64512).length(), ' ');
                ;

            if (in_F4 == 1) // If F4 flag is provided, calculate F4.
            {
				
				// ... //

                F4_value = calc_F4(count_solvent, count_solute, My_neigh, atom_positions, boxX, boxY, boxZ, Nneigh, Natoms, topSolute, time, rawFilename, HBOND_DIST);
                outFile2 << F4_value << std::endl;
                
                // ... //

                //These conditonals are not required
                //if(F4_value > 0) outFile_F4 << frameCounter << "\t| " << F4_value << "\t| " << time << endl;
                //else outFile_F4 << frameCounter << "\t|" << F4_value << "\t| " << time << endl;        //The if-else condition takes care of the extra space needed for "-" sign and alligns the output(lazy way!).
                
                // ... //
            }
            //atom_positions.clear();
            arr_str2.clear();
        }

        atom_positions.clear();
        //arr_str2.clear();
        //solute_atoms.clear();
        solute_positions.clear();
        iter++;
        cout << " Finished " << iter << " number of iterations " << endl;
        cout <<"\n";
        //cout << "iter = " << iter << endl;
    } // End of reading the input file.


    // Reading all the *cage* files
    cout << "Reading cage files for caged solute atoms info " << endl;


    //fstream OFILE(rawFilename+"_cageF4_info.dat", std::ios::app);
    fstream OFILE(rawFilename+"_guest_info.dat", std::ios::app);
    OFILE << "time" << std::string(5,' ');
    for(auto it : real_map)
        OFILE << "512"+it.first+"X" << std::string(5,' ');
    for(auto it : real_map)
        OFILE << "62512"+it.first+"X" << std::string(5,' ');
    for(auto it : real_map)
        OFILE << "64512"+it.first+"X" << std::string(5,' ');

    for(auto it : real_map)
        OFILE << it.first + "X" << std::string(5, ' ');
    //OFILE.close();
    //OFILE << std::string(5, ' ') << "512_total";
    //OFILE << std::string(5, ' ') << "62512_sum";
    //OFILE << std::string(5, ' ') << "64512_sum";
    OFILE << "\n";

    int last = 1;
    for( int f = 1; f <= frameCounter; f++)
    {
        fstream concat_file512(inputFilename.substr(0,inputFilename.length()-4) + "_cage512_concat.gro", std::ios::app);
        fstream concat_file62512(inputFilename.substr(0,inputFilename.length()-4) + "_cage62512_concat.gro", std::ios::app);
        fstream concat_file64512(inputFilename.substr(0,inputFilename.length()-4) + "_cage64512_concat.gro", std::ios::app);

        //fstream OFILE("caged.info", std::ios::app);
        OFILE << f << std::string(9-std::to_string(f).length(), ' ');
        //OFILE.close();
        cout << "frame#: " << f << endl;
        //if(cage_512_count)
        //
            std::string cage_file = inputFilename.substr(0,inputFilename.length()-4) + "_cage512-" + std::to_string(f) + ".gro";
            std::map<std::string, int> caged_map = func(real_map, cage_file);
            std::map<std::string, int> caged_512_map = func(real_map, cage_file);
            
            cout << "caged512_map = " << caged_map << endl;

            //Write to file caged.info

            std::string cage_type = "512";
            for(auto it : caged_map)
            {
                int int1 = cage_type.size() + it.first.size() + 5 - std::to_string(it.second).size();
                if(cage_512_count)
                        OFILE << it.second << std::string(int1, ' ');
                else OFILE << "0" << std::string(int1, ' ');
            }

            // concat all the cage512-*.gro files into one file
        //if(cage_512_count)
        {
            ifstream cage512_file_stream(cage_file);
            std::string line;
            while(std::getline(cage512_file_stream, line))
                concat_file512 << line << endl;
        }
        //if(cage_62512_count)
        //
            //std::string cage_file = inputFilename.substr(0,inputFilename.length()-4) + "_cage62512-" + std::to_string(f) + ".gro";
            cage_file = inputFilename.substr(0,inputFilename.length()-4) + "_cage62512-" + std::to_string(f) + ".gro";
            std::map<std::string, int> caged_62512_map = func(real_map, cage_file);
            caged_map = func( real_map, cage_file);
            cout << "caged62512_map = " << caged_map << std::endl;
            
            //Write to caged.info 
            //std::string cage_type = "62512";
            cage_type = "62512";
            for(auto it : caged_62512_map)
            {
                int int1 = cage_type.size() + it.first.size() + 5 - std::to_string(it.second).size();
                if(cage_62512_count)
                        OFILE << it.second << std::string(int1, ' ');
                else OFILE << "0" << std::string(int1, ' ');
                //OFILE << it.second << std::string(int1, ' ');
            }

        //if(cage_62512_count)
        {
            // concat all the cage62512-*.gro files into one file
            ifstream cage62512_file_stream(cage_file);
            std::string line;
            while(std::getline(cage62512_file_stream, line))
                concat_file62512 << line << endl;
        }
        //if(cage_64512_count)
        //
            //std::string cage_file = inputFilename.substr(0,inputFilename.length()-4) + "_cage64512-" + std::to_string(f) + ".gro";
            cage_file = inputFilename.substr(0,inputFilename.length()-4) + "_cage64512-" + std::to_string(f) + ".gro";
            std::map<std::string, int> caged_64512_map = func(real_map, cage_file);
            caged_map = func(real_map, cage_file);
            cout << "caged64512_map = " << caged_map << endl;

            //Write to caged.info 
            //std::string cage_type = "64512";
            cage_type = "64512";
            for(auto it : caged_64512_map)
            {
                int int1 = cage_type.size() + it.first.size() + 5 - std::to_string(it.second).size();
                if(cage_64512_count)
                        OFILE << it.second << std::string(int1, ' ');
                else OFILE << "0" << std::string(int1, ' ');
                //OFILE << it.second << std::string(int1, ' ');
            }

        //if(cage_64512_count)
        {
            // concat all the cage64512-*.gro files into one file
            ifstream cage64512_file_stream(cage_file);
            std::string line;
            while(std::getline(cage64512_file_stream, line))
                concat_file64512 << line << endl;
        }
        
        
        //fstream OFILE("caged.info", std::ios::app);
        //OFILE << std::endl;
        //OFILE.close();
        for(auto it = caged_512_map.begin(); it != caged_512_map.end(); it++)
        {
            std::string key = it->first;
            int tot_count = caged_512_map[key] + caged_62512_map[key] + caged_64512_map[key];
            cout << "1. " << caged_512_map[key];
            cout << "1. " << caged_62512_map[key];
            cout << "1. " << caged_64512_map[key];
            OFILE << tot_count << std::string(key.length()+5-to_string(tot_count).length(), ' ');
            
        }
        OFILE <<"\n";

        concat_file512.close();
        concat_file62512.close();
        concat_file64512.close();

        cout << "\n";

        last = f;

    }
    
    std::map<int, double> f4_vs_z = calc_time_avg_f4_z("binning_z.xvg", last);
    //for(auto it1: f4_vs_z)
        //it1.second = it1.second/frameCounter;
    cout << f4_vs_z;

    ofstream obj("time_avergaed_f4.xvg");
    for(auto it1 : f4_vs_z)
        obj <<  it1.first << " " << it1.second << std::endl;
    obj.close();
    OFILE.close();

    fileIN.close();
    outFile.close();
    outFile2.close();


    return 0;

} // Close main function
