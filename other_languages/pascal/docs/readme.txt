   This is a readme file for GenAlg.
   
   Author         - Alan D. Freed, Ph.D
   Email          - afreed AT tamu.edu
                  - alan.d.freed AT gmail.com
   Home page      - none
   License        - GNU Lesser General Public License, vs. 3 or later
   Copyright      - (c) Alan D. Freed 2014-2015  
   Pascal Version - 1.3

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
   
   Why write GenAlg in Pascal?  Well, mostly because it is a language that
   thinks the way I think.  But the following quote sums this up nicely:
   "Pascal is a language with a remarkably clear syntax, good control struc-
   tures, and some extraordinarily user-friendly compiler implementations. ...
   Because the language is so strongly structured, Pascal programs often
   execute correctly the first time, while FORTRAN and C programs almost
   never do!"  W.H. Press et al., "Numerical Recipes in Pascal", Cambridge 
   University Press, Cambridge, 1989.

   My original Zonnon implementation of this genetic algorithm drew heavily
   upon the Pascal code SGA.  SGA stands for Simple Genetic Algorithm. 
   SGA was written by David Goldberg and is documented in his first book:
      Goldberg, D.E., "Genetic Algorithms in Search, Optimization &
                       Machine Learning", Addison-Wesley, Boston, 1989.
   David Goldberg's original SGA code carries the following copywrite notice:
      { A Simple Genetic Algorithm - SGA - v1.0 }
      { (c)   David Edward Goldberg  1986       }
      {       All Rights Reserved               }
   GenAlg is a complete rewrite of Goldberg's original SGA code, broadly
   enhancing its capabilities, its statistical assessments, its data structure
   and so on.  Nevertheless, many of his ideas are still present.  Formulae
   used to determine the population size, the number of contestants for
   tournament play, and the number of generations needed for convergence came
   from his second book:
      Goldberg, D.E., "The Design of Innovation: Lessons from and for
                       Competent Genetic Algorithms," Kluwer, Boston, 2002.  

--------------------------------------------------------------------------------    
