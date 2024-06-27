program demo1;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Interfaces, // this includes the LCL widgetset
  Forms, tachartlazaruspkg, demo1form,
  { you can add units after this }
  gaSolver;

{$R *.res}

  const
    reportName = 'linear_test_case';

begin
  RunGeneticAlgorithm(reportName);
  RequireDerivedFormResource := True;
  Application.Initialize;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.

