(* *****************************************************************************
   Author         - Alan D. Freed, Ph.D
   License        - GNU Lesser General Public License, vs. 3 or later
   Copyright      - (c) Alan D. Freed 2014-2016
   Pascal Version - 1.3.3
--------------------------------------------------------------------------------
   GenAlg is a free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the Free
   Software Foundation, either version 3 of the License, or (at your option)
   any later version.

   GenAlg is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
   more details.

   You should have received a copy of the GUN Lesser General Public License
   along with GenAlg.  If not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------------------
   Port began on - April    11, 2014
   Last modified - November 21, 2015
--------------------------------------------------------------------------------
   Version 1.2 - 1.3
   Adopted 'Arrays' for creating and managing vector arrays.  
   This library calls librarys 'Arrays' and 'rng' written by the author.
--------------------------------------------------------------------------------
   Writes out a report for the user, post analysis.
***************************************************************************** *)

// global directives
{$MODE OBJFPC}      // allows objects, classes, interfaces, exception handling
{$H+}               // use ANSI strings
{$IFDEF UNIX}
  {$CALLING CDECL}  // uses the GCC calling convention  - for Unix libraries
  {$PIC ON}         // create Position Independent Code - for Unix libraries
{$ENDIF}

unit gaWriteReport;

interface

  uses
    {$IFDEF UNIX}{$IFDEF UseCThreads}
      cthreads,
    {$ENDIF}{$ENDIF}
    gaColonies; // the genetic algorithm's class for managing the population

  procedure ReportHeader (c : TColony; reportFileName : String);
  { writes out the header to the report file }

  procedure ReportBody (c : TColony; reportFileName : String);
  { writes data from each generation to the report file }

  procedure ReportFooter (c : TColony; reportFileName : String);
  { writes out the footer to the report file }

implementation

  uses
    Classes,      // for constructing objects
    SysUtils,     // the system's file handling routines
    gaCore,       // genetic algorithm's basic type definitions
    gaCreatures;  // the genetic algorithm's class for individuals

  procedure ReportHeader (c : TColony; reportFileName : String);
    var
      i, j, k, l,
      populationSize : Integer;
      report         : TextFile;
  begin
    AssignFile(report, reportFileName);
    {$I+} // use exceptions
    try
      ReWrite(report);   // create the file
      WriteLn(report, '');
      Write(report, '-------------------------------------');
      WriteLn(report, '---------------------------------------');
      Write(report, '    An optimization of ', c.dimR, ' random ');
      if c.dimR = 1 then
        Write(report, 'variable extracted from ')
      else
        Write(report, 'variables extracted from ');
      if c.dimE = 1 then
        WriteLn(report, c.dimE, ' experiment ')
      else
        WriteLn(report, c.dimE, ' experiments ');
      Write(report, '    comprising of ', c.dimS, ' total data points ');
      Write(report, 'fit with a model in ', c.dimP);
      if c.dimP = 1 then
        WriteLn(report, ' parameter')
      else
        WriteLn(report, ' parameters');
      Write(report, '    using a genetic algorithm with a population of size ');
      c.PopulationData(populationSize, i, j, k, l);
      WriteLn(report, populationSize, '.');
      Write(report, '-------------------------------------');
      WriteLn(report, '---------------------------------------');
      Write(report, '    |    elite    | -----------------');
      WriteLn(report, '- population fitness ------------------');
      Write(report, ' #  |   fitness   | --- mean --- | - ');
      WriteLn(report, 'std dev - | - skewness - | - kurtosis -');
      CloseFile(report)
    except
      on E: EInOutError do begin
        Writeln('File handling error occurred. Details: '
                + E.ClassName + '/' + E.Message);
      end;
    end
  end;

  procedure ReportBody (c : TColony; reportFileName : String);
    var
      elite       : TCreature;
      generation  : Integer;
      i, j, k, l  : Integer;
      kurt, mu,
      sigma, skew : Real;
      report      : TextFile;
      s           : String;
  begin
    c.PopulationData(i, generation, j, k, l);
    if generation < 10 then
      s := '  '
    else if generation < 100 then
      s := ' '
    else
      s := '';
    s := s + IntToStr(generation);
    s := s + '   ';
    elite := c.BestCreature;
    s := s + FloatToStrF(elite.GetFitness, ffExponent, 7, 1);
    s := s + '      ';
    c.FitnessMoments(mu, sigma, skew, kurt);
    s := s + FloatToStrF(mu, ffExponent, 3, 1);
    s := s + '       ';
    s := s + FloatToStrF(sigma, ffExponent, 3, 1);
    if skew < 0.0 then
      s := s + '       '
    else
      s := s + '        ';
    s := s + FloatToStrF(skew, ffExponent, 3, 1);
    if kurt < 0.0 then
      s := s + '       '
    else
      s := s + '        ';
    s := s + FloatToStrF(kurt, ffExponent, 3, 1);
    {$I+}
    try
      AssignFile(report, reportFileName);
      Append(report);
      WriteLn(report, s);
      CloseFile(report)
    except
      on E: EInOutError do begin
        Writeln('File handling error occurred. Details: '
                + E.ClassName + '/' + E.Message);
      end;
    end;
    // clean up
    elite := nil
  end;

  procedure ReportFooter (c : TColony; reportFileName : String);
    var
      creatures, i  : Integer;
      crossovers    : Integer;
      distribution  : Integer;
      generations   : Integer;
      immigrants    : Integer;
      mutations, j  : Integer;
      delta, gamma  : Real;
      kurtosis      : Real;
      lambda, mu    : Real;
      lowerBound    : TRVector;
      skewness      : Real;
      sigma, xi     : Real;
      elite         : TCreature;
      parameters    : TRVector;
      report        : TextFile;
      s             : String;
      upperBound    : TRVector;
      varyParameter : TBVector;
  begin
    lowerBound := nil;
    parameters := nil;
    upperBound := nil;
    {$I+}
    try
      AssignFile(report, reportFileName);
      Append(report);
      Write(report, '-  -  -  -  -  -  -  -  -  -  -  -  -');
      WriteLn(report, '  -  -  -  -  -  -  -  -  -  -  -  -  -');
      s := '    Statistics from this optimization run include:';
      WriteLn(report, s);
      c.PopulationData(creatures, generations, immigrants, crossovers, mutations);
      s := '        number of generations = ';
      s := s + IntToStr(generations);
      WriteLn(report, s);
      s := '        number of creatures   = ';
      s := s + IntToStr(creatures);
      WriteLn(report, s);
      s := '        number of immigrants  = ';
      s := s + IntToStr(immigrants);
      WriteLn(report, s);
      s := '        number of crossovers  = ';
      s := s + IntToStr(crossovers);
      WriteLn(report, s);
      s := '        number of mutations   = ';
      s := s + IntToStr(mutations);
      WriteLn(report, s);
      Write(report, '-  -  -  -  -  -  -  -  -  -  -  -  -');
      WriteLn(report, '  -  -  -  -  -  -  -  -  -  -  -  -  -');
      s := '    Fitness statistics from the last generation are:';
      WriteLn(report, s);
      Write(report, '    | Normal Distribution || --------');
      WriteLn(report, '----- Johnson   Distribution ----------');
      Write(report, '    |   mean   |  std dev || S? |  gamma ');
      WriteLn(report, '  |  delta   |    xi    |  lambda');
      c.FitnessMoments(mu, sigma, skewness, kurtosis);
      c.FitnessStatistics(distribution, gamma, delta, xi, lambda);
      if mu < 0.0 then
        s := '     '
      else
        s := '      ';
      s := s + FloatToStrF(mu, ffExponent, 3, 2);
      s := s + '   ';
      s := s + FloatToStrF(sigma, ffExponent, 3, 2);
      if distribution = SB then
        s := s + '    SB'
      else if distribution = SL then
        s := s + '    SL'
      else if distribution = SN then
        s := s + '    SN'
      else if distribution = SU then
        s := s + '    SU'
      else
        s := s + '    ST';
      if gamma < 0.0 then
        s := s + '  '
      else
        s := s + '   ';
      s := s + FloatToStrF(gamma, ffExponent, 3, 2);
      s := s + '   ';
      s := s + FloatToStrF(delta, ffExponent, 3, 2);
      if xi < 0.0 then
        s := s + '  '
      else
        s := s + '   ';
      s := s + FloatToStrF(xi, ffExponent, 3, 2);
      s := s + '   ';
      s := s + FloatToStrF(lambda, ffExponent, 3, 2);
      WriteLn(report, s);
      Write(report, '-  -  -  -  -  -  -  -  -  -  -  -  -');
      WriteLn(report, '  -  -  -  -  -  -  -  -  -  -  -  -  -');
      s := '    Parameters from the elite creature are:';
      WriteLn(report, s);
      s := ' #       lower      optimum      upper      name of parameter';
      WriteLn(report, s);
      elite := c.BestCreature;
      elite.GetAllParameters(parameters);
      varyParameter := c.VaryParameters;
      c.GetBoundaries(lowerBound, upperBound);
      j := 0;
      for i := 1 to c.dimP do begin
        if i < 10 then
          s := ' '
        else
          s := '';
        s := s + IntToStr(i) + '  ';
        if not varyParameter[i] then 
          s := s + '            '
        else begin
          Inc(j);
          if lowerBound[j] >= 0 then
            s := s + '   '
          else
            s := s + '  ';
          s := s + FloatToStrF(lowerBound[j], ffExponent, 4, 2)
        end;
        if parameters[i] >= 0 then
          s := s + '   '
        else
          s := s + '  ';
        s := s + FloatToStrF(parameters[i], ffExponent, 4, 2);
        if not varyParameter[i] then
          s :=  s + '            '
        else begin
          if upperBound[j] >= 0 then
            s := s + '   '
          else
            s := s + '  ';
          s := s + FloatToStrF(upperBound[j], ffExponent, 4, 2)
        end;
        s := s + '    ';
        s := s + c.names[i];
        WriteLn(report, s)
      end;
      Write(report, '-------------------------------------');
      WriteLn(report, '---------------------------------------');
      WriteLn(report, '');
      CloseFile(report)
    except
      on E: EInOutError do begin
        Writeln('File handling error occurred. Details: '+E.ClassName+'/'+E.Message);
      end;
    end
  end;

end.
