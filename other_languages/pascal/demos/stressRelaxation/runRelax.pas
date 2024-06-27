
// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

program RunRelax;

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    Classes, SysUtils, gaRelaxSolver;
    
  procedure Run;
    var
      reportName : String;
  begin
    reportName := 'fitRelaxData';
    RunGeneticAlgorithm(reportName)
  end;
  
begin
  Run
end.
