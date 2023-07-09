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
        void Showsolution() const { if(solved) pLogger->Write(solution_ss.str()); }
            
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
    
        void ResetConstants();
};


///////////////////////////////
/// INITIALIZING    ////
/////////////////////////////

LinearSystem::LinearSystem(Logger* pLog)
{
    pLogger = pLog;
    if(pLog == nullptr)
        exit(1);    //force to have some kind of input/output method
}

void LinearSystem::AddEquation(string eq)
{
    if(eq.size() != 0) // only allow adding equations if matrix has not been solved yet, do nothing for empty string
    {
        if(numEq >= MAX_EQS - 1)
            pLogger->Write("maximum number of equations reached");
        else if(EquationToMatrix(eq))
        {
            numEq++;
            numVars = (int) vars.size();
            //pLogger->Write("enter additional equations or press enter to quit: \n");
        }
        else
        {
            CleanUp();  //undo any changes made in EquationToMatrix
            pLogger->Write( "Please enter a valid linear equation" );
        }
    }
}

void LinearSystem::FinishAndSolve()
{
    for(int i = 0; i < numEq; i++)
        augMatrix[i][numVars] = constants[i];
    
    if(numEq != 0)
    {
        workRecord << "Augmented Matrix: \n";
        workRecord << GetMatrix();
        Solve();
    }
}

bool LinearSystem::ParseNum(size_t& pos, string& num, const string& s) const
{
    size_t startPos = pos;
    while(pos < s.size() && (isdigit(s.at(pos)) || s.at(pos) == '.') ) //put number into a string, iterating until the end of the number
    {
        if( isdigit(s.at(pos)) )
            num += s.at(pos);
        else if (s.at(pos) == '.' && num.find('.') == string::npos) //if it's a decimal and there haven't been any decimals yet
            num += s.at(pos);
        else
            return false;   // 2 decimals in one number
            
        pos++;    //look at next char
    }
    return (pos != startPos);   // if there was no number to parse, pos will be unchanged, otherwise if it succesfully parses, it will return true as pos will not be the same anymore
}

bool LinearSystem::EquationToMatrix(string& s)
{
    size_t eqPos = s.find('=');
    if (eqPos == string::npos)     // no = in equation
        return false;

    size_t i;
    int multiplier = 1; //positive or negative
    bool adjacent = false;  //keep track of whether or not variables are adjacent to each other or not, ex: xy is not allowed but x+y is
    
    for(i = 0; i < eqPos; i++)
    {
        if(! ( isdigit(s.at(i)) || s.at(i) == '.' || isalpha(s.at(i)) || isspace(s.at(i)) || s.at(i) == '+' || s.at(i) == '-' ) ) //anything other than the characters allowed
            return false;
        else if(s.at(i) == '-') //next number will be negative
        {
            adjacent = false;
            multiplier *= -1;
        }
        else if(s.at(i) == '+')
            adjacent = false;
        else if ( isalpha(s.at(i)) )    //variable without coefficient
        {
            if(adjacent)
                return false;
            if (vars.find(s.at(i)) == string::npos )
                vars += s.at(i);
            augMatrix[numEq][vars.find(s.at(i))] += multiplier;
            multiplier = 1; //resets multiplier
            adjacent = true;
        }
        else if(s.at(i) != ' ')    //it's a digit or decimal
        {
            if(adjacent)    //disallow cases like 3x2y, need 3x+2y or 3x-2y
                return false;

            string num;
            if ( ParseNum(i, num, s) /*valid number*/ && isalpha(s.at(i)))   //need var immediately after number, 9 x not allowed, must be 9x
            {
                if (vars.find(s.at(i)) == string::npos )   //variable letter that is not already logged in the vars string
                    vars += s.at(i);
                augMatrix[numEq][vars.find(s.at(i))] += multiplier * stod(num);  //adds to whatever coefficient this variable already has in this equation
                multiplier = 1; //resets multiplier
                adjacent = true;
            }
            else
                return false;   //no variable after number or not a valid number
        }
    }
    
    if(vars.size() == 0 || vars.size() > MAX_VARS)    // no variables or too many variables
        return false;
    
    return (ProcessEnd(s, eqPos));
}

bool LinearSystem::ProcessEnd(string& s,size_t pos)
{
    if(pos == s.size() - 1) // nothing after equal sign
        return false;
    
    pos++;  //right after equal sign
    
    int multiplier = 1;
    
    while( pos < s.size() - 1 && (isspace(s.at(pos)) || s.at(pos) == '-')  )       //iterate until first non-space, non '-' after the '='
    {
        if(s.at(pos) == '-')
            multiplier *= -1;   //if there is a leading negative or negatives, resolve that issue, allow ax + by ... = ----5, = 5
        pos++;
    }
    
    if(!(isdigit(s.at(pos)) || s.at(pos) == '.') )  //first non-space, non '-' must be a digit or decimal
        return false;
    
    string num;
    
    if(ParseNum(pos, num, s))   // pos is now the index after the end of the number
    {
        while( s.at(s.size()-1) == ' ' )
            s.pop_back();
        
        if( pos == s.size() )
        {
            constants[numEq] = (stod(num) * multiplier);    //add this number to the constants array
            return true;
        }
    }
    
    return false;
}

void LinearSystem::CleanUp()  //procedure for returning false from EquationToMatrix, reset the row of augmatrix
{
    for(size_t i = 0; i <= vars.size(); i++)
        augMatrix[numEq][i] = 0;
    
    while(vars.size() > numVars)    //if any new variables in faulty eq, remove them
        vars.pop_back();
}

string LinearSystem::GetMatrix() const
{
    stringstream ss;
    for(int i = 0; i < numEq; i++)
    {
        for(int j = 0; j < numVars + 1; j++)
        {
            ss << left << setw(10) << setprecision(3) << fixed << augMatrix[i][j] << ' ';
        }
        ss << endl;
    }
    return ss.str() + "\n\n";
}

//////////////////////
/// WORK    ////
/////////////////////

void LinearSystem::SetShowWork(bool show_wrk)
{
    showWork = show_wrk;
    DisplayWork();
}

void LinearSystem::DisplayWork() const
{
    if(showWork && solved)
        pLogger->Write(workRecord.str());
}

////////////////////////////////////
/////   TRAVERSAL   ////////
//////////////////////////////////

int LinearSystem::FindLead(int row) const
{
    for(int i = 0; i < numVars; i++)
    {
        if(abs(augMatrix[row][i]) > EPSILON)
            return i;
    }
    return -1;
}

bool LinearSystem::IsLeading(int row, int col) const
{
    for(int i = 0; i < col; i++)
    {
        if(abs(augMatrix[row][i]) > EPSILON)
            return false;
    }
    return abs(augMatrix[row][col]) > EPSILON;
}

////////////////////////////////////
/////   OPERATIONS   //////
/////////////////////////////////

void LinearSystem::R_plus_xR(int r1, int r2, double multiplier)
{
    if(multiplier != 0)
    {
        workRecord << "R" << r1 + 1 << " + " << multiplier << "*R" << r2 + 1 << ":\n\n";
        
        for(int i = 0; i <= numVars; i++)
        {
            augMatrix[r1][i] += augMatrix[r2][i] * multiplier;
            if(abs(augMatrix[r1][i]) < EPSILON) //avoid negative zero
                augMatrix[r1][i] = 0.0;
        }
        workRecord << GetMatrix();
    }
}

void LinearSystem::R_divide(int row, double divisor)
{
    if(divisor == 0)
    {
        cout << "cannot divide by zero, something's wack" << endl;
        exit(1);
    }
    if(divisor != 1)
    {
        workRecord << "Divide R" << row + 1 << " by " << divisor << ":\n\n";

        for(int i = 0; i <= numVars; i++)
        {
            augMatrix[row][i] /= divisor;
            if(abs(augMatrix[row][i]) < EPSILON) //avoid negative zero
                augMatrix[row][i] = 0.0;
        }
        workRecord << GetMatrix();
    }
}

void LinearSystem::SwapRows(int r1, int r2)   // a nice way with pointers?
{
    workRecord << "Swap row " << r1 + 1 << " and " << r2 + 1 << ":\n\n";

    for(int i = 0; i <= numVars; i++)
    {
        double hold = augMatrix[r1][i];
        augMatrix[r1][i] = augMatrix[r2][i];
        augMatrix[r2][i] = hold;
    }
    workRecord << GetMatrix();
}

//////////////////////////////////////////////////
/////   SOLUTION HANDLING   ////////
/////////////////////////////////////////////////

int LinearSystem::Solutions(int row) const    // return 0 for no solutions, 1 for normal case, and 2 for infinite solutions
{
    int lead = FindLead(row);
    if( lead == -1) //row of 0 coefficients
    {
        if(abs(augMatrix[row][numVars]) < EPSILON && row < numVars)  // constant is 0 and this row corresponds to one of the variables
            return 2;   //infinitely many solutions
        else if (abs(augMatrix[row][numVars]) > EPSILON)  //constant is nonzero, 0x = !0, no solutions
            return 0;
    }
    else
    {
        int i;
        for(i = lead + 1; i < numVars; i++) // lead guaranteed to be non-negative and less than numVars
        {
            if(abs(augMatrix[row][i]) > EPSILON)
                return 2;   // free variable
        }
    }
    return 1;   //nothing wrong, 1 solution in that row
}

void LinearSystem::DealWithInfinity(int row)  // at this point matrix should be upper triangular, infinite solutions only allowed in rows that are supposed to have pivot variables
{
    if(abs(augMatrix[row][row]) > EPSILON)
    {
        solution_ss << vars[row] << " = " << augMatrix[row][numVars];
        for(int col = row + 1; col < numVars; col++)
        {
            if(abs(augMatrix[row][col]) > EPSILON)
                solution_ss << " + " << -1 * augMatrix[row][col] << vars[col];
        }
        solution_ss << endl;
    }
    else
        solution_ss << vars[row] << " is a free variable." << endl;
}

//////////////////////////////////////////////////
/////   ROW ECHELON FORM   ///////
////////////////////////////////////////////////

void LinearSystem::REFWide()  //place wide matrix into Row Echelon Form
{
    int rowsComplete = 0;
    int col = 0;
    int flipRow;
    
    while(rowsComplete < numEq)
    {
        flipRow = -1;   //just a value to initialize, can't be confused for a row or for numEq
        
        if(abs(augMatrix[rowsComplete][col]) < EPSILON && col < numVars)    // if 0 coefficient in column col
        {
            for(flipRow = rowsComplete + 1; flipRow < numEq; flipRow++) // all preceeding rows have pivot variable with only 0s below, so start iterating at rowsComplete + 1
            {
                if(IsLeading(flipRow, col))
                    break;
            }
            if (flipRow == numEq)
                col++;  //no leading digit for this var anywhere, look at next variable
            else
                SwapRows(rowsComplete, flipRow);    // swap rows so we can get a leading 1
        }
        
        if(col == numVars)
            return;  //used all variables
        
        if(flipRow != numEq)    //no flip was necessary or succesful flip performed to get leading digit
        {
            R_divide(rowsComplete, augMatrix[rowsComplete][col]);  //divide row to get leading 1, just ensured no division by 0
            
            for(int i = rowsComplete + 1; i < numEq; i++)
                R_plus_xR(i, rowsComplete, -1 * augMatrix[i][col]);    //ensures all subsequent rows have 0 in that column
            
            rowsComplete++;
            col++;
        }
    }
}

void LinearSystem::REFSquareOrTall()  //Place non-wide matrix into Row Echelon Form
{
    int varsComplete = 0;
    while(varsComplete < numVars)
    {
        int flipRow = -1;   // junk value to initialize
        if( abs(augMatrix[varsComplete][varsComplete]) < EPSILON )
        {
            // traverse rows until there is a leading digit in the desired column
            bool noFlipFound = true;
            for(flipRow = 1;flipRow < varsComplete; flipRow++) // no sense flipping 0th row, first var entered will be col 1 so guaranteed leading digit in row 0
            {
                if(IsLeading(flipRow, varsComplete))
                {
                    noFlipFound = false;
                    break;
                }
            }
            
            if(noFlipFound)
            {
                for(flipRow = varsComplete + 1; flipRow < numEq; flipRow++)
                {
                    if(IsLeading(flipRow, varsComplete))
                        break;
                }
            }
            if (flipRow == numEq) //column of all 0s below, don't need to worry about this variable anymore as it's free
                varsComplete++;
            else
                SwapRows(varsComplete, flipRow);
        }
        
        if(flipRow != numEq)    // no swapping needed or succesful swapping
        {
            R_divide(varsComplete, augMatrix[varsComplete][varsComplete]);  //divide row to get leading 1, just ensured no division by 0
        
            for(int i = varsComplete + 1; i < numEq; i++)
                R_plus_xR(i, varsComplete, -1 * augMatrix[i][varsComplete]);    //ensures all subsequent rows have 0 in that column
            
            varsComplete++;
        }
    }
}

///////////////////////////////
/////   SOLVING   ///////
////////////////////////////

void LinearSystem::SolveWide()
{
    REFWide();
    // left with upper triangular-ish matrix
    
    int leads[MAX_EQS];
    leads[numEq - 1] = FindLead(numEq - 1);

    for(int i = numEq - 1; i > 0; i--)  // make sure all leading digits are the only nonzeros in their columns
    {
        leads[i - 1] = FindLead(i - 1);     // find lead in preceeding row, lead at index 0 will be useful
        
        if(leads[i] != -1)  // not row of all 0s
        {
            for(int j = 0; j < i; j++)
                R_plus_xR(j, i, -1 * augMatrix[j][leads[i]]);
        }
    }

    string freeVars;
    for(int row = 0; row < numEq; row++)
    {
        if(leads[row] == -1 && abs(augMatrix[row][numVars]) > EPSILON)     // 0 = !0
        {

            SetNoSolutions();
            return;
        }
        
        solution_ss << vars[ leads[row] ] << " = " << augMatrix[row][numVars];
        for(int c = leads[row] + 1; c < numVars; c++)
        {
            if(abs(augMatrix[row][c]) > EPSILON)
            {
                solution_ss << " + " << -1 * augMatrix[row][c] << vars[c];
                
                string new_free_var;
                new_free_var.push_back(vars[c]);
                new_free_var += " is a free variable.\n";
                
                if(freeVars.find((new_free_var)) == string::npos)
                    freeVars += new_free_var;
            }
        }
        solution_ss << endl;
    }
    
    solution_ss << freeVars;
}

void LinearSystem::SolveNonWide()
{
    REFSquareOrTall();
    // left with upper triangular matrix
    for(int i = 0; i < numEq; i++)
    {
        for(int j = i + 1; j < numVars; j++)
        {
            if(IsLeading(j, j))
                R_plus_xR(i, j, -1 * augMatrix[i][j]);  //take care of upper triangle to put into RREF
        }
    }
    
    for(int i = 0; i < numVars; i++)  //store solutions, avoid jumping to conclusions before confirming that there are solutions
    {
        solutionsArr[i] = Solutions(i);
        if(solutionsArr[i] == 0)
        {
            SetNoSolutions();
            return;
        }
        
        switch(solutionsArr[i]) //now that there are guaranteed to be 1 or more solutions, deal with each case accordingly
        {
            case 1:
                solution_ss << vars[i] << " = " << augMatrix[i][numVars] << endl;
                break;
            case 2:
                DealWithInfinity(i);
        }
    }
    
    for(int j = numVars; j < numEq; j++)
    {
        if(Solutions(j) == 0)
        {
            SetNoSolutions();
            return;
        }
    }
}

void LinearSystem::ResetConstants()
{
    for(int i = 0; i < numEq; i++)
    {
        constants[i] = augMatrix[i][numVars];
        augMatrix[i][numVars] = 0;
    }
}

void LinearSystem::Solve()
{
    if(numEq >= numVars)    // not a wide Matrix
        SolveNonWide();
    
    else    // less rows than vars, either no solution or infinitely many solutions, wide Matrix
        SolveWide();
    
    solved = true;
    
    DisplayWork();
    pLogger->Write(solution_ss.str());
    
    // prepare for solving system again
    equationsSolved = numEq;
    workRecord.str("");
    workRecord.clear();
    solution_ss.str("");
    solution_ss.clear();
    
    ResetConstants();
}

#endif /* Header_h */
