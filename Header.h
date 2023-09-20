//
//  Header.h
//  Haykcpp
//
//  Created by Hayk Gargaloyan on 4/8/23.
//



// NEXT STEP parametrize infinitely many solutions?

#ifndef Header_h
#define Header_h

#include <iostream>
#include <cctype>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>

using namespace std;


const int MAX_VARS = 20;
const int MAX_EQS = 25;
const double EPSILON = 1E-14;

/*
 
 3v+4w+5x-6y+z=10
 v-w+x+y-z=12
 4v-3w+6x-5y+12z=9
 5v+4w+12x-y+3z=3
 
 5v+4w+12x-y=3
 y=69
 4w-6y=40
 3w+4y=39

 */

// option to show steps or not?

// rename class

// Write validate function (???)

// initialize from main, using addEq function

// function to erase from console???

// let user decide decimal precision and setwidth?

class Logger        // rename??
{
public:
    virtual void Write(string s) = 0;
    // virtual void Input(string& s) = 0; let them handle input in main or do through logger?
};



class LinearSystem
{
    public:
        LinearSystem(Logger* pLog);
        //void Initialize();
        void AddEquation(string eq);
        void FinishAndSolve();  //solving same matrix twice? should be ok
        void SetShowWork(bool show_wrk);
        void Showsolution() const { if(solved) pLogger->Write(solution_ss.str()); }     //Write solution again if desired
            
    private:
        int numEq = 0;      // number of rows
        int numVars = 0;    // number of columns - 1
        string vars;        // stores all the variables
        int equationsSolved = 0;
        Logger* pLogger;
        double constants[MAX_EQS];      // number after = sign
        double augMatrix[MAX_EQS][MAX_VARS + 1] = {0}; // actual dimensions numVar + 1 by numEq
        int solutionsArr[MAX_EQS];
        bool solved = false;
        bool showWork = false;
        stringstream workRecord;
        stringstream solution_ss;

        bool ParseNum(size_t& pos, string& num, const string& s) const;
        bool EquationToMatrix(string& s);
        bool ProcessEnd(string& s, size_t pos);
        void CleanUp();
        void ResetConstants();
    
        void R_plus_xR(int r1, int r2, double multiplier);
        void R_divide(int row, double divisor);
        void SwapRows(int r1, int r2);
    
        string GetMatrix() const;     // make public? don't want it being called until after initialization
        
        bool IsLeading(int row, int pos) const;
        int FindLead(int row) const;
    
        void Solve();
        void SolveWide();
        void SolveNonWide();
        void REFWide();
        void REFSquareOrTall();
        
        int Solutions(int row) const;
        void DealWithInfinity(int row);
    
        void DisplayWork() const;
        void SetNoSolutions() { solution_ss.str(""); solution_ss.clear(); solution_ss << "No Solutions\n"; }

};


#endif /* Header_h */
