/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 M. J. Rincón
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Runtime expression parser for symbolic regression expressions.
    
    Uses ExprTK library for fast expression evaluation.
    
    To enable ExprTK:
    1. Download ExprTK: git clone https://github.com/ArashPartow/exprtk.git external/exprtk
    2. Define EXPRTK_AVAILABLE in Make/options or add -DEXPRTK_AVAILABLE to EXE_INC
    3. Add -Iexternal/exprtk to EXE_INC in Make/options

\*---------------------------------------------------------------------------*/

#include "expressionParser.H"
#include "dictionary.H"
#include "IOstreams.H"
#include "dimensionedScalar.H"
#include "wordList.H"
#include <cctype>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

expressionParser::expressionParser(const std::string& expr)
:
    expressionStr_(expr),
#ifdef USE_EXPRTK
    expression_(new exprtk::expression<double>()),
    parser_(new exprtk::parser<double>()),
    symbolTable_(new exprtk::symbol_table<double>()),
#endif
    compiled_(false)
{
    if (!expressionStr_.empty())
    {
        compileExpression();
    }
}


expressionParser::expressionParser(const dictionary& dict)
:
    expressionStr_(""),
#ifdef USE_EXPRTK
    expression_(new exprtk::expression<double>()),
    parser_(new exprtk::parser<double>()),
    symbolTable_(new exprtk::symbol_table<double>()),
#endif
    compiled_(false)
{
    read(dict);
}


expressionParser::~expressionParser()
{}


void expressionParser::setExpression(const std::string& expr)
{
    expressionStr_ = expr;
    compiled_ = false;
    compileExpression();
}


void expressionParser::registerFieldVariable
(
    const std::string& name,
    const volScalarField& field
)
{
    fieldVariables_[name] = &field;
}


void expressionParser::registerScalarVariable(const std::string& name, double* value)
{
    scalarVariables_[name] = value;
}


void expressionParser::registerConstant
(
    const std::string& name,
    const dimensionedScalar& value
)
{
    registerConstant(name, value.value());
}


void expressionParser::registerConstant(const std::string& name, double value)
{
    constants_[name] = value;
}


bool expressionParser::compileExpression()
{
    if (expressionStr_.empty())
    {
        compiled_ = false;
        return false;
    }
    
    // Basic validation: check for empty expression
    if (expressionStr_.find_first_not_of(" \t\n\r") == std::string::npos)
    {
        compiled_ = false;
        return false;
    }
    
#ifdef USE_EXPRTK
    // Clear previous state
    symbolTable_->clear();
    expression_->release();
    
    // Register constants
    for (const auto& constPair : constants_)
    {
        symbolTable_->add_constant(constPair.first, constPair.second);
    }
    
    // Register scalar variables (for single-value evaluation)
    for (const auto& varPair : scalarVariables_)
    {
        symbolTable_->add_variable(varPair.first, *varPair.second);
    }
    
    // Register placeholder variables for field variables that will be set during evaluation
    // This allows compilation even if variables haven't been registered yet
    // We'll use a static map to store placeholder values (persists across calls)
    static std::map<std::string, double> placeholderValues;
    
    // Check which variables are referenced in the expression but not yet registered
    // Common invariant names: I1, I2, I3, I4, I5
    // Register ALL common variables that appear in the expression as placeholders
    // This allows compilation even if they haven't been registered as fields yet
    std::vector<std::string> commonVars = {"I1", "I2", "I3", "I4", "I5"};
    for (const auto& varName : commonVars)
    {
        // Check if variable appears in expression (using word boundary check)
        // Look for varName as a standalone word (not part of another word like "I10")
        std::string pattern = varName;
        // Simple check: look for varName followed by non-alphanumeric or end of string
        bool foundInExpression = false;
        size_t pos = 0;
        while ((pos = expressionStr_.find(varName, pos)) != std::string::npos)
        {
            // Check if it's a standalone word
            bool isStart = (pos == 0) || !std::isalnum(expressionStr_[pos - 1]);
            bool isEnd = (pos + varName.length() == expressionStr_.length()) ||
                         !std::isalnum(expressionStr_[pos + varName.length()]);
            if (isStart && isEnd)
            {
                foundInExpression = true;
                break;
            }
            pos += varName.length();
        }
        
        if (foundInExpression)
        {
            // Check if not already registered as field or scalar variable
            bool notRegistered = (fieldVariables_.find(varName) == fieldVariables_.end() &&
                                  scalarVariables_.find(varName) == scalarVariables_.end());
            
            if (notRegistered)
            {
                // Initialize placeholder value if needed
                if (placeholderValues.find(varName) == placeholderValues.end())
                {
                    placeholderValues[varName] = 0.0;
                }
                // Add to symbol table (ExprTK will use this placeholder during compilation)
                // The actual values will be set during evaluateField() when fields are registered
                // Note: ExprTK's add_variable takes a reference, so the value must persist
                // Using static map ensures the value persists across function calls
                double& varRef = placeholderValues[varName];
                symbolTable_->add_variable(varName, varRef);
            }
        }
    }
    
    // Also register any field variables that are already registered (for compilation)
    // These need placeholder values too since we can't use field pointers during compilation
    for (const auto& fieldPair : fieldVariables_)
    {
        // Check if not already added above
        bool alreadyAdded = false;
        for (const auto& varName : commonVars)
        {
            if (varName == fieldPair.first)
            {
                alreadyAdded = true;
                break;
            }
        }
        
        if (!alreadyAdded)
        {
            // Initialize placeholder value if needed
            if (placeholderValues.find(fieldPair.first) == placeholderValues.end())
            {
                placeholderValues[fieldPair.first] = 0.0;
            }
            // Add to symbol table
            symbolTable_->add_variable(fieldPair.first, placeholderValues[fieldPair.first]);
        }
    }
    
    // Set symbol table in expression
    expression_->register_symbol_table(*symbolTable_);
    
    // Compile expression
    if (!parser_->compile(expressionStr_, *expression_))
    {
        WarningInFunction
            << "Failed to compile expression: " << expressionStr_ << nl
            << "Parser error: " << parser_->error() << endl;
        compiled_ = false;
        return false;
    }
    
    compiled_ = true;
    return true;
#else
    // Without ExprTK, just mark as compiled if non-empty
    WarningInFunction
        << "ExprTK not available. Expression will not be evaluated." << nl
        << "To enable ExprTK:" << nl
        << "  1. Download: git clone https://github.com/ArashPartow/exprtk.git external/exprtk" << nl
        << "  2. Add -DEXPRTK_AVAILABLE to EXE_INC in Make/options" << nl
        << "  3. Add -Iexternal/exprtk to EXE_INC in Make/options" << endl;
    compiled_ = false;
    return false;
#endif
}


double expressionParser::evaluate() const
{
    if (!compiled_)
    {
        WarningInFunction
            << "Expression not compiled. Returning 0.0."
            << endl;
        return 0.0;
    }
    
#ifdef USE_EXPRTK
    if (expression_)
    {
        return expression_->value();
    }
    else
    {
        WarningInFunction
            << "Expression object not valid. Returning 0.0."
            << endl;
        return 0.0;
    }
#else
    WarningInFunction
        << "ExprTK not available. Cannot evaluate expression."
        << endl;
    return 0.0;
#endif
}


void expressionParser::evaluateField(volScalarField& result) const
{
    if (!compiled_)
    {
        WarningInFunction
            << "Expression not compiled. Setting field to zero."
            << endl;
        result = dimensionedScalar("zero", result.dimensions(), 0.0);
        return;
    }
    
#ifdef USE_EXPRTK
    if (!expression_ || !symbolTable_)
    {
        WarningInFunction
            << "Expression object not valid. Setting field to zero."
            << endl;
        result = dimensionedScalar("zero", result.dimensions(), 0.0);
        return;
    }
    
    // Check if field variables are registered
    if (fieldVariables_.empty())
    {
        WarningInFunction
            << "No field variables registered. Cannot evaluate expression."
            << " Make sure to call registerFieldVariable() before evaluateField()."
            << endl;
        result = dimensionedScalar("zero", result.dimensions(), 0.0);
        return;
    }
    
    // Create symbol table with variable references for efficient field evaluation
    exprtk::symbol_table<double> cellSymbolTable;
    
    // Register constants (same for all cells)
    for (const auto& constPair : constants_)
    {
        cellSymbolTable.add_constant(constPair.first, constPair.second);
    }
    
    // Create variable storage for field variables (one per field variable)
    // Store values in a way that ExprTK can reference
    std::vector<double> fieldVarStorage(fieldVariables_.size());
    std::map<std::string, double*> fieldVarPtrs;
    
    label varIdx = 0;
    for (const auto& fieldPair : fieldVariables_)
    {
        fieldVarPtrs[fieldPair.first] = &(fieldVarStorage[varIdx]);
        // Register variable by reference - ExprTK will use this pointer
        cellSymbolTable.add_variable(fieldPair.first, fieldVarStorage[varIdx]);
        varIdx++;
    }
    
    // Create and compile expression once
    exprtk::expression<double> cellExpression;
    cellExpression.register_symbol_table(cellSymbolTable);
    
    exprtk::parser<double> cellParser;
    if (!cellParser.compile(expressionStr_, cellExpression))
    {
        WarningInFunction
            << "Failed to compile expression for field evaluation: " << expressionStr_ << nl
            << "Parser error: " << cellParser.error() << endl;
        result = dimensionedScalar("zero", result.dimensions(), 0.0);
        return;
    }
    
    // Evaluate expression for each cell
    const label nCells = result.size();
    
    // Debug: check if we have field variables registered
    if (fieldVariables_.empty())
    {
        WarningInFunction
            << "No field variables registered for expression evaluation!" << nl
            << "Expression: " << expressionStr_ << nl
            << "Make sure registerFieldVariable() is called before evaluateField()." << endl;
        result = dimensionedScalar("zero", result.dimensions(), 0.0);
        return;
    }
    
    // Debug: print first evaluation for verification
    static label evalCallCount = 0;
    evalCallCount++;
    if (evalCallCount == 1)
    {
        Info<< "[EXPRESSION_PARSER] Evaluating expression: " << expressionStr_ << nl
            << "    Registered field variables: ";
        for (const auto& fieldPair : fieldVariables_)
        {
            Info<< fieldPair.first << " ";
        }
        Info<< nl << "    Number of cells: " << nCells << endl;
    }
    
    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        // Update field variable values for this cell
        for (const auto& fieldPair : fieldVariables_)
        {
            const volScalarField* field = fieldPair.second;
            if (field && cellI < field->size() && fieldVarPtrs.find(fieldPair.first) != fieldVarPtrs.end())
            {
                // Update value through pointer - ExprTK sees this immediately
                *(fieldVarPtrs[fieldPair.first]) = (*field)[cellI];
            }
            else
            {
                // Set to zero if field is invalid or cell index out of range
                if (fieldVarPtrs.find(fieldPair.first) != fieldVarPtrs.end())
                {
                    *(fieldVarPtrs[fieldPair.first]) = 0.0;
                }
            }
        }
        
        // Evaluate expression for this cell
        result[cellI] = cellExpression.value();
        
        // Debug: print first cell evaluation
        if (evalCallCount == 1 && cellI == 0)
        {
            Info<< "[EXPRESSION_PARSER] First cell evaluation:" << nl
                << "    Cell 0 result: " << result[cellI] << nl
                << "    Variable values: ";
            for (const auto& fieldPair : fieldVariables_)
            {
                if (fieldVarPtrs.find(fieldPair.first) != fieldVarPtrs.end())
                {
                    Info<< fieldPair.first << "=" << *(fieldVarPtrs[fieldPair.first]) << " ";
                }
            }
            Info<< endl;
        }
    }
    
    // Correct boundary conditions
    result.correctBoundaryConditions();
#else
    WarningInFunction
        << "ExprTK not available. Cannot evaluate expression field."
        << " Setting field to zero." << endl;
    result = dimensionedScalar("zero", result.dimensions(), 0.0);
#endif
}


void expressionParser::read(const dictionary& dict)
{
    expressionStr_ = dict.lookupOrDefault<string>("expression", "");
    
    // Remove any trailing quotes, semicolons, or whitespace from expression string
    // This handles cases where the dictionary might have trailing characters
    while (!expressionStr_.empty() && 
           (expressionStr_.back() == '"' || expressionStr_.back() == '\'' || 
            expressionStr_.back() == ';' || expressionStr_.back() == ' ' || 
            expressionStr_.back() == '\t'))
    {
        expressionStr_.pop_back();
    }
    // Also remove leading quotes if present
    if (!expressionStr_.empty() && (expressionStr_[0] == '"' || expressionStr_[0] == '\''))
    {
        expressionStr_ = expressionStr_.substr(1);
    }
    
    // Read normalization constants from variables sub-dictionary
    if (dict.found("variables"))
    {
        const dictionary& varsDict = dict.subDict("variables");
        wordList keys = varsDict.toc();
        forAll(keys, i)
        {
            word key = keys[i];
            scalar value = varsDict.lookupOrDefault<scalar>(key, 0.0);
            registerConstant(key, value);
        }
    }
    
    if (!expressionStr_.empty())
    {
        compileExpression();
    }
    else
    {
        compiled_ = false;
    }
}


void expressionParser::write(Ostream& os) const
{
    os.writeKeyword("expression") << expressionStr_ << token::END_STATEMENT
        << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
