(* *****************************************************************************
  genAlg Author         - Alan D. Freed
  genAlg License        - GNU Lesser General Public License, vs 3 or later
  genAlg Copyright      - (c) Alan D. Freed 2014
  genAlg Pascal Version - 1.3.3
--------------------------------------------------------------------------------
  Use a genetic algorithm for the purpose of parameter estimation.
--------------------------------------------------------------------------------
  This application uses three external libraries, all written by the author:
     'Arrays' provides a collection of array type and procedures to manage them
     'rng'    provides a random number generator and statistical functions.
     'GenAlg' provides a genetic algorithm for use in parameter estimation.
--------------------------------------------------------------------------------
   This is one of five parts to the template.  It produces the program to run.
--------------------------------------------------------------------------------
   References
      Goldberg, D.E., Genetic Algorithms in Search, Optimization, and Machine
         Learning, Addison-Wesley, Boston, 1989.
      Goldberg, D.E., The Design of Innovation: Lessons learned from and for
         competent genetic algorithms.  In: Genetic algorithms and evolutionary
         computation, Vol. 7, Klewer, Boston, 2002.
      Johnson, N.L., "Systems of Frequency Curves Generated by Methods of
         Translation," Biometrika, Vol. 36 (1949), 149-176.       
***************************************************************************** *)

// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

program gaDriver;

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    Classes, SysUtils, gaSolver;
    
  procedure Run;
    var
      reportName : String;
  begin
    // enter the name of the text file that the report is to be written to
    reportName := '...';
    // run the genetic algorithm
    RunGeneticAlgorithm(reportName)
    
    (* the following functions are also available to be called
    
    GetEliteData (var parameters : TRVector; var genome : TSVector);

    GetModelData (selectExperiment : Integer;             // 1..dimE
                  var xExperiment  : TRMatrix;            // output
                  var yExperiment  : TRMatrix;            // output
                  var xModel       : TRMatrix;            // input
                  var yModel       : TRMatrix);           // output         *)
  end;
  
begin
  Run
end.
